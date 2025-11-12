use num_complex::Complex;
use std::f64::consts::PI;
use pyo3::prelude::*;


fn elementwise_mul(a: &[Complex<f64>], b: &[Complex<f64>]) -> Vec<Complex<f64>> {
    assert_eq!(a.len(), b.len(), "Arrays must have the same length");
    a.iter().zip(b.iter()).map(|(x, y)| x * y).collect()
    
}


fn gaussian(sig: f64, n: usize) -> Vec<Complex<f64>> {
    let mut out: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); (n) as usize];
    for i in -((n/2) as isize)..(n/2) as isize {
        out[((n/2) as isize + i) as usize] = Complex::new((1.0/(sig*(2.0*PI).sqrt()))*(-0.5*((i as f64)/sig).powf(2.0)).exp(), 0.0)
    }
    out
}

fn circ_convolve(arr: Vec<Complex<f64>>, filter: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    // let mut out: Vec<f64> = vec![0.0; arr.len()];
    let mut fil = filter.clone();
    for i in 0..fil.len() {
        fil[i] = filter[(filter.len() + 1)%(filter.len())]
    }
    let outy = idft_loc(elementwise_mul(&dft_loc(arr), &dft_loc(filter)));

    outy
}


fn dft_loc(signal: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
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
    out
}


fn idft_loc(signal: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
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
    out
}


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


#[pyfunction]
pub fn filter(signal: Vec<Complex<f64>>, freq: f64) -> PyResult<Vec<Complex<f64>>> {
    // let mut out: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); signal.len()];
    let sig = 1.0 / freq;
    let gausy = gaussian(sig, signal.len());
    let out = circ_convolve(signal, gausy);

    Ok(out)
}