use pyo3::prelude::*;


fn matmul(a: &Vec<Vec<f64>>, b: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let n = a.len();
    let m = b[0].len();
    let mut result = vec![vec![0.0; m]; n];
    for i in 0..n {
        for j in 0..m {
            for k in 0..b.len() {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    result
}

fn mat_vec_mul(a: &Vec<Vec<f64>>, v: &[f64]) -> Vec<f64> {
    let n = a.len();
    let mut r = vec![0.0; n];
    for i in 0..n {
        r[i] = doty(&a[i], v);
    }
    r
}

fn identity(n: usize) -> Vec<Vec<f64>> {
    let mut m = vec![vec![0.0; n]; n];
    for i in 0..n { m[i][i] = 1.0; }
    m
}

fn transpose(a: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let n = a.len();
    let m = a[0].len();
    let mut t = vec![vec![0.0; n]; m];
    for i in 0..n {
        for j in 0..m {
            t[j][i] = a[i][j];
        }
    }
    t
}

fn norm(v: &Vec<f64>) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn dot(a: &Vec<f64>, b: &Vec<f64>) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn doty(u: &[f64], v: &[f64]) -> f64 {
    u.iter().zip(v.iter()).map(|(a, b)| a * b).sum()
}

fn householder_tridiagonal(mut a: Vec<Vec<f64>>) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut q = vec![vec![0.0; n]; n];
    for i in 0..n {
        q[i][i] = 1.0;
    }

    for k in 0..(n - 2) {
        let mut x = vec![0.0; n - k - 1];
        for i in 0..(n - k - 1) {
            x[i] = a[k + 1 + i][k];
        }

        let alpha = -x[0].signum() * norm(&x);
        x[0] -= alpha;
        let norm_x = norm(&x);
        if norm_x.abs() < 1e-12 {
            continue;
        }
        for xi in &mut x {
            *xi /= norm_x;
        }

        // Build Householder matrix H = I - 2vv^T
        let mut H = vec![vec![0.0; n]; n];
        for i in 0..n {
            H[i][i] = 1.0;
        }

        for i in 0..(n - k - 1) {
            for j in 0..(n - k - 1) {
                H[k + 1 + i][k + 1 + j] -= 2.0 * x[i] * x[j];
            }
        }

        // Apply H to A and accumulate Q
        a = matmul(&matmul(&transpose(&H), &a), &H);
        q = matmul(&q, &H);
    }

    (a, q)
}


fn qr_decompose(a: &Vec<Vec<f64>>) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut q = vec![vec![0.0; n]; n];
    let mut r = vec![vec![0.0; n]; n];

    let mut a_cols: Vec<Vec<f64>> = (0..n).map(|j| (0..n).map(|i| a[i][j]).collect()).collect();

    for j in 0..n {
        let mut v = a_cols[j].clone();

        for i in 0..j {
            r[i][j] = dot(&q.iter().map(|row| row[i]).collect(), &a_cols[j]);
            for k in 0..n {
                v[k] -= r[i][j] * q[k][i];
            }
        }

        r[j][j] = norm(&v);
        for k in 0..n {
            q[k][j] = v[k] / r[j][j];
        }
    }

    (q, r)
}

fn jacobi(mut a: Vec<Vec<f64>>, tol: f64, max_sweeps: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut v = identity(n);

    let off_norm = |m: &Vec<Vec<f64>>| -> f64 {
        let mut s = 0.0;
        for i in 0..n {
            for j in (i+1)..n {
                s += 2.0 * m[i][j] * m[i][j];
            }
        }
        s.sqrt()
    };

    for _sweep in 0..max_sweeps {
        let mut max_off = 0.0;
        for p in 0..n {
            for q in (p+1)..n {
                let apq = a[p][q].abs();
                if apq > max_off { max_off = apq; }
                if apq <= tol { continue; }

                let app = a[p][p];
                let aqq = a[q][q];
                let tau = (aqq - app) / (2.0 * a[p][q]);
                let t = if tau >= 0.0 {
                    1.0 / (tau + (1.0 + tau*tau).sqrt())
                } else {
                    -1.0 / (-tau + (1.0 + tau*tau).sqrt())
                };
                let c = 1.0 / (1.0 + t*t).sqrt();
                let s = t * c;

                let app_new = c*c*app - 2.0*s*c*a[p][q] + s*s*aqq;
                let aqq_new = s*s*app + 2.0*s*c*a[p][q] + c*c*aqq;
                a[p][p] = app_new;
                a[q][q] = aqq_new;
                a[p][q] = 0.0;
                a[q][p] = 0.0;

                for k in 0..n {
                    if k != p && k != q {
                        let akp = a[k][p];
                        let akq = a[k][q];
                        a[k][p] = c*akp - s*akq;
                        a[p][k] = a[k][p];
                        a[k][q] = s*akp + c*akq;
                        a[q][k] = a[k][q];
                    }
                }

                for k in 0..n {
                    let vkp = v[k][p];
                    let vkq = v[k][q];
                    v[k][p] = c*vkp - s*vkq;
                    v[k][q] = s*vkp + c*vkq;
                }
            }
        }
        if max_off < tol || off_norm(&a) < tol { break; }
    }

    let eigvals = (0..n).map(|i| a[i][i]).collect::<Vec<f64>>();
    (eigvals, v)
}

fn gen_q(n: usize) -> Vec<Vec<f64>> {
    let mut result = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            // cast to signed before subtracting to avoid underflow
            if (i as isize - j as isize).abs() == 1 {
                result[i][j] = 0.5 * (( (i + j + 1) as f64 ).sqrt());
            }
        }
    }
    result
}

fn gen_q_2(n: usize) -> Vec<Vec<f64>> {
    let mut result = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                result[i][j] = 0.5 * ((2 * j + 1) as f64);
            } else if i + 2 == j {
                // here j >= 2 whenever this branch is true, but use f64 to compute safely
                let jf = j as f64;
                result[i][j] = 0.5 * ( (jf * (jf - 1.0)).sqrt() );
            } else if i == j + 2 {
                let jf = j as f64;
                result[i][j] = 0.5 * ( ((jf + 2.0) * (jf + 1.0)).sqrt() );
            }
        }
    }
    result
}

fn gen_q_4(n: usize) -> Vec<Vec<f64>> {
    let mut result = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            let jf = j as f64;
            if i == j {
                result[i][j] = 0.75 * (2.0 * jf * jf + 2.0 * jf + 1.0);
            } else if i + 2 == j {
                // compute polynomial in f64 (use powi)
                let num = 2.0 * jf.powi(3) - 3.0 * jf.powi(2) + 1.0 * jf;
                let den = (jf * (jf - 1.0)).sqrt();
                // protect against division by zero just in case (shouldn't happen for valid j)
                result[i][j] = if den != 0.0 { 0.5 * (num / den) } else { 0.0 };
            } else if i == j + 2 {
                result[i][j] = 0.5 * (((jf + 2.0) * (jf + 1.0)).sqrt() * (2.0 * jf + 3.0));
            } else if i + 4 == j {
                let num = jf.powi(4) - 6.0 * jf.powi(3) + 11.0 * jf.powi(2) - 6.0 * jf;
                let den = (jf * (jf - 1.0) * (jf - 2.0) * (jf - 3.0)).sqrt();
                result[i][j] = if den != 0.0 { 0.25 * (num / den) } else { 0.0 };
            } else if i == j + 4 {
                result[i][j] = 0.25 * (((jf + 4.0) * (jf + 3.0) * (jf + 2.0) * (jf + 1.0)).sqrt());
            }
        }
    }
    result
}



#[pyfunction]
fn q_4(n: usize, tajp: String) -> PyResult<Vec<Vec<f64>>> {
    if tajp == "1" {
        let q = gen_q(n);
        Ok(matmul(&matmul(&q, &q), &matmul(&q, &q)))
    }
    else if tajp == "2" {
        let qq = gen_q_2(n);
        Ok(matmul(&qq, &qq))
    }
    else if tajp == "4" {
        Ok(gen_q_4(n))
    }
    else {Ok(vec![vec![0.0; n]; n])}
}

#[pyfunction]
fn matsum(a: Vec<Vec<f64>>, b: Vec<Vec<f64>>) -> PyResult<Vec<Vec<f64>>> {
    let n = b[0].len();
    let mut m = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            m[i][j] = a[i][j] + b[i][j]
        }
    }
    Ok(m)
}

#[pyfunction]
fn mul(a: f64, b: Vec<Vec<f64>>) -> PyResult<Vec<Vec<f64>>> {
    let n = b[0].len();
    let mut m = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            m[i][j] = a * b[i][j]
        }
    }
    Ok(m)
}

#[pyfunction]
fn harm(n: usize) -> PyResult<Vec<Vec<f64>>> {
    let mut m = vec![vec![0.0; n]; n];
    for i in 0..n {
        m[i][i] = i as f64 + 0.5
    }
    Ok(m)
}

#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}


#[pyfunction]
fn echo_copy(matrix: Vec<Vec<f64>>) -> PyResult<Vec<Vec<f64>>> {

    Ok(matrix)
}


#[pyfunction]
fn power_method(matrix: Vec<Vec<f64>>, max_iters: usize, tol: f64) -> PyResult<(f64, Vec<f64>)> {
    let n = matrix.len();
    let mut b_k = vec![1.0; n];

    for _ in 0..max_iters {
        // Multiply A * b_k
        let mut b_k1 = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                b_k1[i] += matrix[i][j] * b_k[j];
            }
        }

        let norm_bk1 = norm(&b_k1);
        if norm_bk1 < 1e-12 {
            break;
        }

        for i in 0..n {
            b_k1[i] /= norm_bk1;
        }

        // Convergence check
        let diff = b_k
            .iter()
            .zip(&b_k1)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);

        b_k = b_k1;

        if diff < tol {
            break;
        }
    }

    // Approximate eigenvalue λ ≈ (b_kᵀ A b_k) / (b_kᵀ b_k)
    let mut ab = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            ab[i] += matrix[i][j] * b_k[j];
        }
    }

    let lambda = dot(&b_k, &ab) / dot(&b_k, &b_k);
    Ok((lambda, b_k))
}


#[pyfunction]
fn jacobi_eigen(matrix: Vec<Vec<f64>>, max_iters: usize, tol: f64) -> PyResult<(Vec<f64>, Vec<Vec<f64>>)> {
    let n = matrix.len();
    let mut a = matrix.clone();
    let mut v = vec![vec![0.0; n]; n];

    // initialize v as identity
    for i in 0..n {
        v[i][i] = 1.0;
    }

    for _ in 0..max_iters {
        // find largest off-diagonal element
        let mut p = 0;
        let mut q = 1;
        let mut max_val = 0.0;

        for i in 0..n {
            for j in (i + 1)..n {
                if a[i][j].abs() > max_val {
                    max_val = a[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }

        if max_val < tol {
            break;
        }

        let theta = 0.5 * ((2.0 * a[p][q]) / (a[q][q] - a[p][p])).atan();
        let cos = theta.cos();
        let sin = theta.sin();

        // rotate A
        for i in 0..n {
            if i != p && i != q {
                let aip = a[i][p];
                let aiq = a[i][q];
                a[i][p] = cos * aip - sin * aiq;
                a[p][i] = a[i][p];
                a[i][q] = sin * aip + cos * aiq;
                a[q][i] = a[i][q];
            }
        }

        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];

        a[p][p] = cos * cos * app - 2.0 * sin * cos * apq + sin * sin * aqq;
        a[q][q] = sin * sin * app + 2.0 * sin * cos * apq + cos * cos * aqq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;

        // update eigenvector matrix
        for i in 0..n {
            let vip = v[i][p];
            let viq = v[i][q];
            v[i][p] = cos * vip - sin * viq;
            v[i][q] = sin * vip + cos * viq;
        }
    }

    let eigvals = (0..n).map(|i| a[i][i]).collect::<Vec<f64>>();
    Ok((eigvals, v))
}

#[pyfunction]
fn lanczos(matrix: Vec<Vec<f64>>, k: usize) -> PyResult<(Vec<f64>, Vec<Vec<f64>>)> {
    let n = matrix.len();
    let mut v_prev = vec![0.0; n];
    let mut v = vec![1.0; n];
    let norm_v = norm(&v);
    for i in 0..n {
        v[i] /= norm_v;
    }

    let mut alphas = Vec::new();
    let mut betas = Vec::new();
    let mut vs = vec![v.clone()];

    for _ in 0..k {
        let mut w = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                w[i] += matrix[i][j] * v[j];
            }
        }

        let alpha = dot(&v, &w);
        alphas.push(alpha);

        for i in 0..n {
            w[i] -= alpha * v[i] + if betas.len() > 0 { betas.last().unwrap() * v_prev[i] } else { 0.0 };
        }

        let beta = norm(&w);
        if beta < 1e-12 {
            break;
        }
        betas.push(beta);

        v_prev = v.clone();
        for i in 0..n {
            v[i] = w[i] / beta;
        }
        vs.push(v.clone());
    }

    // Build the tridiagonal matrix T
    let m = alphas.len();
    let mut t = vec![vec![0.0; m]; m];
    for i in 0..m {
        t[i][i] = alphas[i];
        if i + 1 < m {
            t[i][i + 1] = betas[i];
            t[i + 1][i] = betas[i];
        }
    }

    // Return T and basis vectors (vs)
    Ok((alphas, t))
}

#[pyfunction]
fn qr_eigen(matrix: Vec<Vec<f64>>, max_iters: usize, tol: f64) -> PyResult<(Vec<f64>, Vec<Vec<f64>>)> {
    let (mut t, mut q_total) = householder_tridiagonal(matrix);
    let n = t.len();

    for _ in 0..max_iters {
        let (q, r) = qr_decompose(&t);
        t = matmul(&r, &q);
        q_total = matmul(&q_total, &q);

        // Convergence check: off-diagonal elements small
        let mut off_diag_sum = 0.0;
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    off_diag_sum += t[i][j].abs();
                }
            }
        }
        if off_diag_sum < tol {
            break;
        }
    }

    let eigvals = (0..n).map(|i| t[i][i]).collect::<Vec<f64>>();
    Ok((eigvals, q_total))
}

#[pyfunction]
fn lanczos_smallest(
    matrix: Vec<Vec<f64>>,
    k: usize,
    m_opt: Option<usize>,       // Lanczos basis size (m >= k). Default: min(n, max(2k, 20))
    _max_iters: Option<usize>,   // overall iteration cap (unused here, but kept)
) -> PyResult<(Vec<f64>, Vec<Vec<f64>>)> {
    // if !is_square(&matrix) {
    //     return Err(PyValueError::new_err("matrix must be square"));
    // }
    let n = matrix.len();
    if n == 0 {
        return Ok((Vec::new(), Vec::new()));
    }
    // if k == 0 || k > n {
    //     return Err(PyValueError::new_err("k must satisfy 1 <= k <= n"));
    // }

    // let sym_err = approx_symmetry_error(&matrix);
    // if sym_err > 1e-8 {
    //     return Err(PyValueError::new_err(format!("matrix is not symmetric (avg abs diff ~ {})", sym_err)));
    // }

    // choose m
    let m = m_opt.unwrap_or_else(|| {
        let guess = (2*k).min(n);
        std::cmp::max(guess, 20.min(n))
    }).min(n);

    // Lanczos storage
    let mut v_space: Vec<Vec<f64>> = Vec::with_capacity(m);
    let mut alpha: Vec<f64> = vec![0.0; m];
    let mut beta: Vec<f64>  = vec![0.0; m+1];

    // initial vector (random-ish deterministic)
    let mut v = vec![1.0f64; n];
    let vnorm = norm(&v);
    // if vnorm == 0.0 { return Err(PyValueError::new_err("initial vector has zero norm")); }
    for x in &mut v { *x /= vnorm; }
    let mut v_prev = vec![0.0f64; n];

    for j in 0..m {
        // w = A v
        let mut w = mat_vec_mul(&matrix, &v);
        let alphaj = doty(&v, &w);
        alpha[j] = alphaj;
        for i in 0..n {
            w[i] -= alphaj * v[i] + beta[j] * v_prev[i];
        }

        // full reorthogonalization against previous v's
        for prev in &v_space {
            let coeff = doty(&w, prev);
            for i in 0..n {
                w[i] -= coeff * prev[i];
            }
        }

        let betaj1 = norm(&w);
        v_space.push(v.clone());
        if betaj1 < 1e-16 {
            // lucky breakdown
            beta[j+1] = 0.0;
            break;
        }
        beta[j+1] = betaj1;
        v_prev = v;
        v = w.iter().map(|x| x / betaj1).collect();
    }

    let actually_m = v_space.len();
    // if actually_m == 0 {
    //     return Err(PyValueError::new_err("Lanczos built zero basis"));
    // }

    // build small tridiagonal T (actually_m x actually_m)
    let mut t = vec![vec![0.0; actually_m]; actually_m];
    for i in 0..actually_m {
        t[i][i] = alpha[i];
        if i + 1 < actually_m {
            t[i][i+1] = beta[i+1];
            t[i+1][i] = beta[i+1];
        }
    }

    // diagonalize T with Jacobi (T is small)
    let (tev, tevvecs) = jacobi(t, 1e-12, 1000);

    // tev holds eigenvalues (unsorted). We want k smallest -> sort ascending by eigenvalue.
    let mut idxs: Vec<usize> = (0..tev.len()).collect();
    idxs.sort_by(|&i, &j| tev[i].partial_cmp(&tev[j]).unwrap());

    let kk = std::cmp::min(k, tev.len());
    let mut eigvals: Vec<f64> = Vec::with_capacity(kk);
    let mut eigvecs: Vec<Vec<f64>> = vec![vec![0.0; kk]; n]; // rows x kk

    for (out_col, &ti) in idxs.iter().take(kk).enumerate() {
        eigvals.push(tev[ti]);
        // z = column ti of tevvecs (note: tevvecs[row][col])
        let mut z = vec![0.0; tevvecs.len()];
        for r in 0..tevvecs.len() {
            z[r] = tevvecs[r][ti];
        }
        // Map back to original space: y = V * z
        let mut y = vec![0.0; n];
        for (j, vj) in v_space.iter().enumerate() {
            let coeff = z[j];
            for i in 0..n {
                y[i] += coeff * vj[i];
            }
        }
        // normalize y
        let ny = norm(&y);
        if ny > 0.0 {
            for i in 0..n { y[i] /= ny; }
        }
        for i in 0..n {
            eigvecs[i][out_col] = y[i];
        }
    }

    Ok((eigvals, eigvecs))
}

#[pyfunction]
fn transposey(a: Vec<Vec<f64>>) -> PyResult<Vec<Vec<f64>>> {
    let n = a.len();
    let m = a[0].len();
    let mut t = vec![vec![0.0; n]; m];
    for i in 0..n {
        for j in 0..m {
            t[j][i] = a[i][j];
        }
    }
    Ok(t)
}


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
    Ok(())
}
