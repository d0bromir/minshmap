# minimizer_ext

A tiny [PyO3](https://pyo3.rs) wrapper that exposes the canonical (w, k)-minimizers
of the [`minimizer-iter`](https://crates.io/crates/minimizer-iter) crate
(rust-seq / Igor Martayan) to Python, so `minshmap.py` can use a real library for its
sketch instead of a hand-rolled one.

It exports a single function:

```python
from minimizer_ext import canonical_minimizers
canonical_minimizers(seq: str, k: int, w: int) -> list[tuple[pos, value, strand]]
```

- `k` — minimizer (k-mer) length.
- `w` — window width in k-mers; **must be odd** (the canonical scheme needs an odd
  window to break forward/reverse-complement ties).
- returns, in position order, one `(position, value, strand)` per emitted minimizer;
  `strand` is `True` when the reverse-complement k-mer was the canonical one.

All the sketch logic (rolling hash, sliding-window minimum, canonicalization) lives in
`minimizer-iter`; this crate only marshals its output into Python tuples.

## Build / install

Requires a Rust toolchain (`rustup`) and [`maturin`](https://github.com/PyO3/maturin)
(`pip install maturin`):

```sh
cd minimizer_ext
maturin build --release -i <python>     # e.g. -i C:\Python314\python.exe
pip install --user target/wheels/minimizer_ext-*.whl
```

During development you can instead run `maturin develop --release` from inside an
activated virtualenv.
