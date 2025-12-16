use pyo3::prelude::*;
mod fn03;
mod fn04;
mod fn05;
mod fn06;
mod fn08;
mod wav;

use crate::fn03::{
    sum_as_string, echo_copy, power_method, jacobi_eigen, lanczos, qr_eigen,
    lanczos_smallest, transposey, q_4, matsum, mul, harm, dodatna,
};

use crate::fn04::{
    dft, idft, gauss, filter
};

use crate::fn05::{
    fft
};

use crate::fn06::{
    euler, analyt, heun, rk2a, rku4, rk45, rkf, pc4
};

use crate::fn08::{shooter};

use crate::wav::branje;


/// A Python module implemented in Rust.
#[pymodule]
fn rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(echo_copy, m)?)?;
    m.add_function(wrap_pyfunction!(power_method, m)?)?;
    m.add_function(wrap_pyfunction!(jacobi_eigen, m)?)?;
    m.add_function(wrap_pyfunction!(lanczos, m)?)?;
    m.add_function(wrap_pyfunction!(qr_eigen, m)?)?;
    m.add_function(wrap_pyfunction!(lanczos_smallest, m)?)?;
    m.add_function(wrap_pyfunction!(transposey, m)?)?;
    m.add_function(wrap_pyfunction!(q_4, m)?)?;
    m.add_function(wrap_pyfunction!(matsum, m)?)?;
    m.add_function(wrap_pyfunction!(mul, m)?)?;
    m.add_function(wrap_pyfunction!(harm, m)?)?;
    m.add_function(wrap_pyfunction!(dodatna, m)?)?;
    m.add_function(wrap_pyfunction!(dft, m)?)?;
    m.add_function(wrap_pyfunction!(gauss, m)?)?;
    m.add_function(wrap_pyfunction!(idft, m)?)?;
    m.add_function(wrap_pyfunction!(branje, m)?)?;
    m.add_function(wrap_pyfunction!(filter, m)?)?;
    m.add_function(wrap_pyfunction!(fft, m)?)?;
    m.add_function(wrap_pyfunction!(euler, m)?)?;
    m.add_function(wrap_pyfunction!(analyt, m)?)?;
    m.add_function(wrap_pyfunction!(heun, m)?)?;
    m.add_function(wrap_pyfunction!(rk2a, m)?)?;
    m.add_function(wrap_pyfunction!(rku4, m)?)?;
    m.add_function(wrap_pyfunction!(rk45, m)?)?;
    m.add_function(wrap_pyfunction!(rkf, m)?)?;
    m.add_function(wrap_pyfunction!(pc4, m)?)?;
    m.add_function(wrap_pyfunction!(shooter, m)?)?;
    Ok(())
}
