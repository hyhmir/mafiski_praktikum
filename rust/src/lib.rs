use pyo3::prelude::*;
mod fn03;
mod fn04;

use crate::fn03::{
    sum_as_string, echo_copy, power_method, jacobi_eigen, lanczos, qr_eigen,
    lanczos_smallest, transposey, q_4, matsum, mul, harm, dodatna,
};

use crate::fn04::{
    dft, idft, gauss
};


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
    Ok(())
}
