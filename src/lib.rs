use pyo3::prelude::*;
use numpy::PyReadonlyArray2;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn diag(matrix: PyReadonlyArray2<f64>) -> PyResult<Vec<f64>> {
    let mat = matrix.as_array();
    Ok(vec![1.0, 2.0, 3.0])
}

/// A Python module implemented in Rust.
#[pymodule]
fn rjast_modul(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(diag, m)?)?;
    Ok(())
}
