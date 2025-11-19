use num_complex::Complex;
use std::f64::consts::PI;
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

fn bit_reverse_permute(a: &mut [Complex<f64>]) {
    let n = a.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            a.swap(i, j);
        }
    }
}

#[pyfunction]
pub fn fft(mut a: Vec<Complex<f64>>, inverse: bool) -> PyResult<Vec<Complex<f64>>> {
    let n = a.len();
    if n == 0 || !n.is_power_of_two() {
        return Err(PyValueError::new_err("length must be a non-zero power of two"));
    }

    // Reorder array by bit-reversed indices (in-place)
    bit_reverse_permute(&mut a);

    // Iterative Cooley-Tukey
    let mut len = 2usize;
    while len <= n {
        // Note: forward FFT uses exp(-2πi/len), inverse uses +2πi/len
        let angle = 2.0 * PI / (len as f64) * if inverse { 1.0 } else { -1.0 };
        let wlen = Complex::new(angle.cos(), angle.sin());

        let mut i = 0usize;
        while i < n {
            let mut w = Complex::new(1.0, 0.0);
            for j in 0..(len / 2) {
                let u = a[i + j];
                let v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w = w * wlen;
            }
            i += len;
        }

        len <<= 1;
    }

    // If inverse, scale by 1/n
    if inverse {
        let inv_n = 1.0 / (n as f64);
        for x in a.iter_mut() {
            x.re *= inv_n;
            x.im *= inv_n;
        }
    }

    Ok(a)
}