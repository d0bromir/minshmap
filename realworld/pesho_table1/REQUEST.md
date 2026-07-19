# Какво иска Пешо (Table 1 от main_1.pdf, benchmark-а на shmap)

## Контекст
Показаната таблица е **Table 1 — "Evaluation comparison of long-read mapping
tools"** от статията на Пешо (`../main_1.pdf`). Тя сравнява long-read mappers
(minimap2, winnowmap2, blend, mapquik, **map-shmap**) върху synthetic и real HiFi
reads, мапнати към човешкия геном (Chromosome Y или All chromosomes).

## Колоните в таблицата
| Колона | Значение |
|---|---|
| **Mapped Q60** | брой рийдове, мапнати като **правилни** с mapq = 60 |
| **Q<60 or missed** | процент рийдове или немапнати, или мапнати с mapq < 60 |
| **Wrong Q60** | брой **грешни** мапвания въпреки mapq = 60 (n/a за real reads) |
| **Index sec** | време (сек) за строене на индекса |
| **Map sec** | време (сек) за мапване на рийдовете |
| **Memory GB** | пикова памет |

## Двете конкретни задачи от Пешо
1. **Направи колоните, които изискват "правилен алайнмънт":**
   `Mapped Q60`, `Q<60 or missed`, `Wrong Q60`.
   Тези три се нуждаят от дефиниция за "правилно мапване".
2. **Раздели `Map time` от `Index time`** (в момента ги имаме заедно като
   `index_s` и `map_s`, но не са свързани с правилността → трябва в един ред от
   таблицата, точно както в статията).

## Дефиниция на "правилно мапване" (задължителна, от Пешо)
> It considers correct a mapping which is overlapping with the ground truth more
> than 10% of the united length of the mapping and the ground truth.

Т.е. **IoU-подобно** правило: мапването е **правилно**, ако
`overlap(mapping, truth) / |mapping ∪ truth| > 0.10`
(припокриване спрямо обединената дължина на мапването и истината, > 10%).

## Целта
Пресъздай **цялата Table 1** доколкото можем с нашите инструменти, за сравнение.
- Задължително минимум срещу **minimap2** (другите тулове са прекалено бавни или
  нестабилни, така че minimap2 е достатъчният baseline).
- `map-shmap` = нашата програма (shmap / minSHmap).

## Изход
Ред за всеки (mapper × dataset) с колоните:
`Mapped Q60 | Q<60 or missed % | Wrong Q60 | Index sec | Map sec | Memory GB`.

## Симулатор за HiFi (основните данни) — от `shmap/Makefile` + main_1.pdf
**PBSIM3** (`~/libs/pbsim3/src/pbsim`). Точната команда от `gen_reads` таргета:

```
pbsim --strategy wgs --method sample --sample HG002.fq.tmp \
      --genome <ref.fa> --depth <DEPTH> --prefix <PREFIX> --no-fastq 1
```

- `--method sample` = **семплира error/quality профила от РЕАЛНИ HG002 HiFi
  рийдове** (`HG002.fq`, първите 400000 реда), не фиксиран error-модел → затова
  синтетичните рийдове имитират разпределението на дължина и грешки на реалния
  HiFi dataset (точно както пише в main_1.pdf: „PBSIM3 -method sample").
- PBSIM3 произвежда `.maf` файлове с **ground-truth alignment** → превръщат се
  във FASTA чрез `paftools.js pbsim2fq <ref.fai> *.maf`, който слага истината в
  header-а: `>name!chr!start!end!strand` (напр. `>S1_1!chr21!12345!29456!+`).
- После `seqkit shuffle` разбърква, и по избор liftOver (chain
  `hg002v1.1_to_CHM13v2.0.chain`) когато рийдовете са симулирани от HG002, а
  референцията е CHM13v2.0 (`READSIM_REFNAME != REFNAME`).
- Coverage: chrY-симулация `DEPTH=2..10x`, All-chromosomes `DEPTH=1x` (виж
  Makefile редове 452-459 и таблицата в статията).

**Извод:** за да пресъздам колоните за „правилен алайнмънт", ми трябват
симулирани рийдове С ground-truth координати. Точният път е PBSIM3 →
`paftools.js pbsim2fq` → header `name!chr!start!end!strand`, който после
сравнявам с PAF-а по IoU > 10%.

## Типове данни за редовете (точно от Table 1 на `main_1.pdf`, Sec. 5.1)

Table 1 има **четири dataset-а** (= четири типа данни за редовете), всеки от които
дава по един блок редове (по един ред за всеки mapper). Референцията винаги е
**CHM13v2.0** (`data_rw/hs1.fa`); `chrY` се извлича от нея. Всеки тип живее в свой
подкаталог под `pesho_table1/data/<тип>/reads.fa`, а генериращият скрипт **спазва**
неговите параметри (референция, дължина, покритие, реален/симулиран).

| # | Тип (dataset в статията) | Референция | Дължина | Покритие | Truth? | Как се генерира |
|---|---|---|---|---|---|---|
| 1 | **Chromosome Y, simulated 10kbp, 10×** | `chrY` (от CHM13v2.0) | ~10 kbp | 10× | да | PBSIM3 `--method sample` от реален HG002 HiFi профил (`hifi_sample.fastq`, mean ~12.8 kb ≈ „10kbp") |
| 2 | **All chromosomes, simulated 10kbp, 1×** | цял CHM13v2.0 | ~10 kbp | 1× | да | PBSIM3 `--method sample` (същия профил), но `--genome hs1.fa --depth 1` |
| 3 | **Chromosome Y, simulated 24kbp, 10×** | `chrY` | ~24 kbp | 10× | да | PBSIM3 `--method errhmm --errhmm ERRHMM-SEQUEL.model --length-mean 24000 --length-sd 4000 --accuracy-mean 0.99` (нашият HiFi семпъл стига само ~15 kb, затова 24 kb се постига чрез length-контролиран PacBio модел) |
| 4 | **All chromosomes, real 24kbp от HG002, 1.6×** | цял CHM13v2.0 | реални | 1.6× | **не** | реални HG002 HiFi рийдове (subset `hifi_sample.fastq` → FASTA); **няма ground truth** → `Wrong Q60 = n/a`, а `Mapped Q60` = броят рийдове с mapq = 60 (без проверка за правилност) |

**Правила, които генераторите СПАЗВАТ (по тип):**
- **Референция**: типове 1 и 3 → `data/_ref/chrY.fa`; типове 2 и 4 → целия
  `hs1.fa` (CHM13v2.0). Ред от chrY НЕ се мери срещу цял геном и обратно.
- **Дължина**: типове 1–2 наследяват дължината от реалния HiFi семпъл (`--method
  sample`, ~13 kb); тип 3 налага mean 24 kb чрез `--length-mean`; тип 4 е с реалната
  дължина на изтеглените рийдове.
- **Покритие (depth)**: chrY → `--depth 10`; all-chromosomes симулация → `--depth 1`;
  реалните рийдове дават ~1.6× чрез subset-а.
- **Truth**: симулираните (1–3) носят `name!chr!start!end!strand` в header-а →
  скорват се по IoU > 10 %. Реалните (4) нямат truth → `Wrong Q60 = n/a`.

**Известни отклонения (честно):** нашият HiFi семпъл е от CCS-15kb библиотека
(mean ~12.8 kb), а не HiFi-revio-24kb. Затова тип 3 („24kbp simulated") ползва
length-контролиран PacBio модел вместо `--method sample`, а тип 4 („real 24kbp")
реално е ~13 kb subset. Осите, които са научно важни за таблицата (референция,
покритие, реални-срещу-симулирани, наличие на ground truth), се възпроизвеждат точно.

