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
   ⇒ прозорецът се отрязва). Най-добрият прозорец → ред PAF.

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

---

Подробното описание на оригиналния shmap: [ORIGINAL_SHMAP.md](ORIGINAL_SHMAP.md).
Бенчмаркове (само C++ и Python minSHmap): [realworld/](realworld/).

