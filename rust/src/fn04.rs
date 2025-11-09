use num_complex::Complex;
use std::f64::consts::PI;
use pyo3::prelude::*;


#[pyfunction]
pub fn dft(signal: Vec<Complex<f64>>) -> PyResult<Vec<Complex<f64>>> {
    let m = signal.len();
    let mut out: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); m];
    for k in 0..m {
        let mut fk = Complex::new(0.0, 0.0);
        for n in 0..m {
            let mk = Complex::new(0.0, - 2.0 * PI * (k as f64) * (n as f64) / (m as f64)).exp();
            fk += mk * signal[n]
        }
        out[k] = fk
    }
    Ok(out)
}


#[pyfunction]
pub fn idft(signal: Vec<Complex<f64>>) -> PyResult<Vec<Complex<f64>>> {
    let m = signal.len();
    let mut out: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); m];
    for k in 0..m {
        let mut fk = Complex::new(0.0, 0.0);
        for n in 0..m {
            let mk = Complex::new(0.0, - 2.0 * PI * (k as f64) * (n as f64) / (m as f64)).exp();
            fk += mk * signal[n]
        }
        out[k] = fk/m as f64
    }
    Ok(out)
}


#[pyfunction]
pub fn gauss(sig: f64, n: usize, split: bool) -> PyResult<Vec<f64>> {
    let mut out: Vec<f64> = vec![0.0; (2*n) as usize];
    for i in -(n as isize)..n as isize {
        out[((n) as isize + i) as usize] = (1.0/(sig*(2.0*PI).sqrt()))*(-0.5*((i as f64)/sig).powf(2.0)).exp()
    }

    if split {
        let v = out.clone();
        for i in 0..2*n {
            out[i] = v[(n + i)%(2*n)]
        }
    }

    Ok(out)
}