//! `minimizer-iter` (rust-seq / Igor Martayan) exposed to BOTH languages so the
//! Python tool (minshmap.py) and the C++ tool (minshmap.cpp) compute the *same*
//! minimizers from the *same* library -- no copy-paste, no divergence.
//!  * Python sees `canonical_minimizers` (PyO3 cdylib, wheel built by maturin).
//!  * C++ sees `mz_compute` / `mz_free` (a C ABI, linked from the static lib).
//! "Canonical" means a sequence and its reverse complement select the same
//! minimizers; `w` must be ODD (the canonical scheme uses an odd window to break
//! the forward/reverse tie at the centre).
use minimizer_iter::MinimizerBuilder;

/// Core: canonical (w, k)-minimizers as (position, value, strand). Empty if seq < k;
/// None if w is invalid (0/even) so each binding reports it its own way.
fn compute(seq: &[u8], k: usize, w: usize) -> Option<Vec<(usize, u64, bool)>> {
    if w == 0 || w % 2 == 0 {
        return None;
    }
    if seq.len() < k {
        return Some(Vec::new());
    }
    Some(
        MinimizerBuilder::<u64>::new()
            .canonical()
            .minimizer_size(k)
            .width(w as u16)
            .iter(seq)
            .map(|(value, pos, strand)| (pos, value, strand))
            .collect(),
    )
}

// ---- Python binding (default features only; maturin builds the cdylib wheel) ----
#[cfg(feature = "python")]
mod py {
    use super::compute;
    use pyo3::exceptions::PyValueError;
    use pyo3::prelude::*;

    /// Canonical (w, k)-minimizers of `seq`: one `(position, value, strand)` per
    /// minimizer in position order (`w` must be odd).
    #[pyfunction]
    fn canonical_minimizers(seq: &str, k: usize, w: usize) -> PyResult<Vec<(usize, u64, bool)>> {
        compute(seq.as_bytes(), k, w).ok_or_else(|| {
            PyValueError::new_err("canonical minimizers require an odd window width w")
        })
    }

    #[pymodule]
    fn minimizer_ext(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add_function(wrap_pyfunction!(canonical_minimizers, m)?)?;
        Ok(())
    }
}

// ---- C ABI (always present; C++ links the static lib, --no-default-features) ----
/// One minimizer for C++: 0-based k-mer start, 64-bit value, strand (0 fw / 1 rc).
#[repr(C)]
pub struct Mz {
    pub pos: u64,
    pub val: u64,
    pub strand: u8,
}

/// Canonical (w, k)-minimizers of `seq[0..len]`; writes count to `out_len`, returns a
/// heap array to read then pass to `mz_free`. Null on invalid w; count 0 when seq < k.
/// # Safety: `seq` must point to `len` bytes; `out_len` must be writable.
#[no_mangle]
pub unsafe extern "C" fn mz_compute(seq: *const u8, len: usize, k: usize, w: usize, out_len: *mut usize) -> *mut Mz {
    let bytes = if len == 0 { &[][..] } else { std::slice::from_raw_parts(seq, len) };
    match compute(bytes, k, w) {
        None => { *out_len = 0; std::ptr::null_mut() }
        Some(v) => {
            let mut boxed: Box<[Mz]> = v.into_iter().map(|(pos, val, s)| Mz { pos: pos as u64, val, strand: s as u8 }).collect();
            *out_len = boxed.len();
            let ptr = boxed.as_mut_ptr();
            std::mem::forget(boxed);
            ptr
        }
    }
}

/// Free an array from `mz_compute`. # Safety: ptr/len must come from `mz_compute`.
#[no_mangle]
pub unsafe extern "C" fn mz_free(ptr: *mut Mz, len: usize) {
    if !ptr.is_null() {
        drop(Box::from_raw(std::slice::from_raw_parts_mut(ptr, len)));
    }
}
