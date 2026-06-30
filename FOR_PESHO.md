# minSHmap — за Пешо

Минимална, четима реализация на sketch-based read mapper (в стила на minSH).
Един файл [minshmap.py](minshmap.py); [minshmap.cpp](minshmap.cpp) е същият
алгоритъм, бит-в-бит същия изход, за реален мащаб.

## Алгоритъмът (три стъпки)

1. **Sketch — canonical (w, k)-minimizers.** Не пишем свой скеч; идва от чужда
   библиотека [`minimizer-iter`](https://crates.io/crates/minimizer-iter)
   (rust-seq / Igor Martayan). За всеки прозорец от `w` k-мера се избира k-мерът с
   най-малък хеш. *Canonical* = рийдът и неговият reverse-complement дават едни и
   същи minimizer-и, така че един индекс обслужва двете вериги. (`w` е нечётно.)

2. **Index.** Минимайзерите на референцията се обръщат в речник
   `hash → [(сегмент, позиция, верига), …]`.

3. **Map.** За всеки рийд: скеч → seed-ове, подредени най-рядък пръв → разпръскване
   в припокриващи се прозорци → оценка на containment с **seed-евристиката**
   `sh = 1 − (използвани − съвпадения) / m` (горна граница на постижимото; `sh < θ`
   ⇒ прозорецът се отрязва). Най-добрият прозорец → ред PAF, плюс качество на
   мапването `mapq ∈ {0, 60}` по **φ-free** правилото (виж по-долу).

Това е цялата идея — ~150 реда, четат се ред по ред. Същата seed-евристика като в
minSH, само че там ограничава edit distance за подравняване, а тук — containment за
намиране на локус.

## Скечът е от библиотека (и за Python, и за C++)

Твоята забележка беше да ползваме чужд minimizer, не copy-paste. Направено:

- **Python** вика библиотеката директно през тънка PyO3 обвивка
  [minimizer_ext/](minimizer_ext/); в [minshmap.py](minshmap.py) остава само:
  ```python
  from minimizer_ext import canonical_minimizers
  def minimizers(seq, k, w):
      for pos, h, strand in canonical_minimizers(seq, k, w):
          yield pos, h, 1 if strand else 0
  ```
- **C++** вика **същия** crate през C ABI (`mz_compute`/`mz_free`), затова двете
  реализации дават **byte-identical** изход (проверено: py == cpp == cpp `-j4`).

## Билд и пускане

```sh
# еднократно: minimizer-ите (иска Rust toolchain + maturin)
cd minimizer_ext
maturin build --release -i <python> && pip install target/wheels/*.whl   # за Python
cargo build --release --no-default-features                              # статична либ за C++

# Python
python minshmap.py ref.fa reads.fa -k 15 -w 11 -t 0.9

# C++
g++ -O3 -std=c++17 -march=native -pthread -I ../shmap/ext/unordered_dense/include \
    -o minshmap minshmap.cpp -L minimizer_ext/target/release \
    -l:libminimizer_ext.a -lws2_32 -luserenv -lbcrypt -lntdll
./minshmap ref.fa reads.fa -k 15 -w 11 -t 0.9
```

`-w` е нечётно (canonical). `-j N` в C++ е само базов read-level паралелизъм
(изходът е същият като `-j 1`), не оптимизация.

## Премахване на параметъра φ (Def. 5)

В статията Def. 5 („алтернативни мапвания") използваше свободен параметър **φ**
(максимално допустимо припокриване) плюс неудобната клауза „различна посока".
Премахнахме φ, като мерим припокриването в **оригиналните T-координати** на
референцията (покрития интервал `g(s)`, точно както е изходът на Sec. 4.6): тогава
reverse-complement образът на най-доброто мапване **съвпада** с него (припокриване
= 1) и автоматично отпада — без клаузата за посока и без параметър.

**Дефиниция (φ-free, Def. 5′).** *Алтернатива* = всеки кандидат, чийто референтен
интервал е **дизюнктен** с този на най-доброто. `mapq = 60`, ако всяка дизюнктна
алтернатива е по-слаба с повече от `δ` (default 0.15); иначе `mapq = 0`. Самият
**placement (позицията) е непроменен** — добавя се само mapq.

Двете програмки (py и cpp) смятат това идентично. Пълни дефиниции/формули, лема и
емпирична валидация (вкл. тест с дублиран блок, който коректно сваля mapq→0):
[NOTE_phi_elimination.md](NOTE_phi_elimination.md).

## Последни резултати

chr21 (T2T-CHM13v2.0, 2000 reads/датасет, `k=15 w=11`); `shmap` е референтният
baseline (2026-06-21, негови собствени настройки):

| датасет | py r/s | cpp r/s | shmap r/s | py mapped | cpp mapped | shmap mapped |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| hifi | 7.0 | 395.6 | 271.0 | 15 | 15 | 17 |
| ont | 1.7 | 154.3 | 69.0 | 32 | 32 | 32 |
| clr | 12.6 | 430.4 | 1142.0 | 2 | 2 | 2 |

Synthetic (8000 къси reads, `k=15 w=11 t=0.5`): py 14 499 r/s, cpp 50 084 r/s, и
двете мапват 3385, placement precision 100 %.

py и cpp мапват **един и същ брой** reads на всеки датасет — алгоритмично
еквивалентни (еднакъв binary-search lookup), byte-identical PAF. C++ изпреварва
`shmap` на hifi (395.6 vs 271 r/s) и ont (154.3 vs 69 r/s); shmap води на късите
clr. Чистият Python е бавен на ont (1.7 r/s) заради ~35 kb reads, но дава **същите**
плейсменти като C++ (не е алгоритмична грешка). Пълни данни:
[realworld/results_rw/realworld.md](realworld/results_rw/realworld.md).

---

Подробното описание на оригиналния shmap: [ORIGINAL_SHMAP.md](ORIGINAL_SHMAP.md).
Бенчмаркове (само C++ и Python minSHmap): [realworld/](realworld/).

