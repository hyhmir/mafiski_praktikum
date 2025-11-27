use pyo3::prelude::*;
use std::f64::consts::PI;
use pyo3::exceptions::PyValueError;


fn f(y: f64, k: f64, yz: f64) -> f64 {
    -k*(y - yz)
}

fn f_ext(t: f64, y: f64, k: f64, yz: f64, a: f64, d: f64) -> f64 {
    -k*(y - yz) + a*((2.0*PI/24.0) * (t - d)).sin()
}

fn analy(t: f64, y0: f64, k: f64, yz: f64) -> f64 {
    yz + (-k*t).exp()*(y0 - yz)
}


/// returns the analytical solution for the basic problem
#[pyfunction(signature = (t, y0, yz =-5.0, k=0.1))]
pub fn analyt(t: Vec<f64>, y0: f64, yz: f64, k: f64) -> PyResult<Vec<f64>> {
    Ok(t.iter().map(|&x| analy(x, y0, k, yz)).collect())
}


#[pyfunction(signature = (t, y0, yz =-5.0, k=0.1, extra=false, a=1.0, d=10.0))]
pub fn euler(t: Vec<f64>, y0: f64, yz: f64, k: f64, extra: bool, a: f64, d: f64) -> PyResult<Vec<f64>> {
    let mut out: Vec<f64> = vec![y0; t.len()];
    if extra {
        for i in 1..t.len() {
            out[i] = out[i-1] + (t[i] - t[i-1])*f_ext(t[i-1], out[i-1], k, yz, a, d)
        }
    }
    else {
        for i in 1..t.len() {
            out[i] = out[i-1] + (t[i] - t[i-1])*f(out[i-1], k, yz)
        }
    }
    Ok(out)
}


#[pyfunction(signature = (t, y0, yz =-5.0, k=0.1, extra=false, a=1.0, d=10.0))]
pub fn heun(t: Vec<f64>, y0: f64, yz: f64, k: f64, extra: bool, a: f64, d: f64) -> PyResult<Vec<f64>> {
    let mut out: Vec<f64> = vec![y0; t.len()];
    let mut h;
    let mut k1;
    let mut k2;
    if extra {
        for i in 1..t.len() {
            h = t[i] - t[i-1];
            k1 = h * f_ext(t[i-1], out[i-1], k, yz, a, d);
            k2 = h * f_ext(t[i], out[i-1] + k1, k, yz, a, d);
            out[i] = out[i-1] + (k1 + k2)/2.0
        }
    }
    else {
        for i in 1..t.len() {
            h = t[i] - t[i-1];
            k1 = h * f(out[i-1], k, yz);
            k2 = h * f(out[i-1] + k1, k, yz);
            out[i] = out[i-1] + (k1 + k2)/2.0
        }
    }
    Ok(out)
}


#[pyfunction(signature = (t, y0, yz =-5.0, k=0.1, extra=false, a=1.0, d=10.0))]
pub fn rk2a(t: Vec<f64>, y0: f64, yz: f64, k: f64, extra: bool, a: f64, d: f64) -> PyResult<Vec<f64>> {
    let mut out: Vec<f64> = vec![y0; t.len()];
    let mut h;
    let mut k1;
    if extra {
        for i in 1..t.len() {
            h = t[i] - t[i-1];
            k1 = h * f_ext(t[i-1], out[i-1], k, yz, a, d);
            out[i] = out[i-1] + h * f_ext(t[i-1] + h/2.0, out[i-1] + k1, k, yz, a, d)
        }
    }
    else {
        for i in 1..t.len() {
            h = t[i] - t[i-1];
            k1 = h * f(out[i-1], k, yz);
            out[i] = out[i-1] + h * f(out[i-1] + k1, k, yz)
        }
    }
    Ok(out)
}


#[pyfunction(signature = (t, y0, yz =-5.0, k=0.1, extra=false, a=1.0, d=10.0))]
pub fn rku4(t: Vec<f64>, y0: f64, yz: f64, k: f64, extra: bool, a: f64, d: f64) -> PyResult<Vec<f64>> {
    let mut out: Vec<f64> = vec![y0; t.len()];
    let mut h;
    let mut k1;
    let mut k2;
    let mut k3;
    let mut k4;
    if extra {
        for i in 1..t.len() {
            h = t[i] - t[i-1];
            k1 = h * f_ext(t[i-1], out[i-1], k, yz, a, d);
            k2 = h * f_ext(t[i-1] + 0.5*h, out[i-1] + 0.5*k1, k, yz, a, d);
            k3 = h * f_ext(t[i-1] + 0.5*h, out[i-1] + 0.5*k2, k, yz, a, d);
            k4 = h * f_ext(t[i], out[i-1] + k3, k, yz, a, d);
            out[i] = out[i-1] + ( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0
        }
    }
    else {
        for i in 1..t.len() {
            h = t[i] - t[i-1];
            k1 = h * f(out[i-1], k, yz);
            k2 = h * f(out[i-1] + 0.5*k1, k, yz);
            k3 = h * f(out[i-1] + 0.5*k2, k, yz);
            k4 = h * f(out[i-1] + k3, k, yz);
            out[i] = out[i-1] + ( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0
        }
    }
    Ok(out)
}


#[pyfunction(signature = (t, y0, yz =-5.0, k=0.1, extra=false, a=1.0, d=10.0))]
pub fn rk45(t: Vec<f64>, y0: f64, yz: f64, k: f64, extra: bool, a: f64, d: f64) -> PyResult<(Vec<f64>, Vec<f64>)> {
    let mut out: Vec<f64> = vec![y0; t.len()];
    let mut err: Vec<f64> = vec![0.0; t.len()];
    let c20 = 1.0/4.0;
    let c30 = 3.0/8.0;
    let c40 = 12.0/13.0;
    let c50 = 1.0;
    let c60 = 1.0/2.0;

    let c21 = 1.0/4.0;
    let c31 = 3.0/32.0;
    let c32 = 9.0/32.0;
    let c41 = 1932.0/2197.0;
    let c42 = -7200.0/2197.0;
    let c43 = 7296.0/2197.0;
    let c51 = 439.0/216.0;
    let c52 = -8.0;
    let c53 = 3680.0/513.0;
    let c54 = -845.0/4104.0;
    let c61 = -8.0/27.0;
    let c62 = 2.0;
    let c63 = -3544.0/2565.0;
    let c64 = 1859.0/4104.0;
    let c65 = -11.0/40.0;

    let a1 = 25.0/216.0;
    let a2 = 0.0;
    let a3 = 1408.0/2565.0;
    let a4 = 2197.0/4104.0;
    let a5 = -1.0/5.0;

    let b1 = 16.0/135.0;
    let b2 = 0.0;
    let b3 = 6656.0/12825.0;
    let b4 = 28561.0/56430.0;
    let b5 = -9.0/50.0;
    let b6 = 2.0/55.0;

    let mut h;
    let mut k1;
    let mut k2;
    let mut k3;
    let mut k4;
    let mut k5;
    let mut k6;
    let mut x5;
    if extra {
        for i in 1..t.len() {
            h = t[i] - t[i-1];
            k1 = h * f_ext(t[i-1], out[i-1], k, yz, a, d);
            k2 = h * f_ext(t[i-1] + c20*h, out[i-1] + c21*k1, k, yz, a, d);
            k3 = h * f_ext(t[i-1] + c30*h, out[i-1] + c31*k1+c32*k2, k, yz, a, d);
            k4 = h * f_ext(t[i-1] + c40*h, out[i-1] + c41*k1+c42*k2+c43*k3, k, yz, a, d);
            k5 = h * f_ext(t[i-1] + c50*h, out[i-1] + c51*k1+c52*k2+c53*k3+c54*k4, k, yz, a, d);
            k6 = h * f_ext(t[i-1] + c60*h, out[i-1] + c61*k1+c62*k2+c63*k3+c64*k4+c65*k5, k, yz, a, d);
            out[i] = out[i-1] + a1*k1 + a2*k2 + a3*k3 + a4*k4 + a5*k5;
            x5 = out[i-1] + b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6;
            err[i] = (x5 - out[i]).abs()
        }
    }
    else {
        for i in 1..t.len() {
            h = t[i] - t[i-1];
            k1 = h * f(out[i-1], k, yz);
            k2 = h * f(out[i-1] + c21*k1, k, yz);
            k3 = h * f(out[i-1] + c31*k1+c32*k2, k, yz);
            k4 = h * f(out[i-1] + c41*k1+c42*k2+c43*k3, k, yz);
            k5 = h * f(out[i-1] + c51*k1+c52*k2+c53*k3+c54*k4, k, yz);
            k6 = h * f(out[i-1] + c61*k1+c62*k2+c63*k3+c64*k4+c65*k5, k, yz);
            out[i] = out[i-1] + a1*k1 + a2*k2 + a3*k3 + a4*k4 + a5*k5;
            x5 = out[i-1] + b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6;
            err[i] = (x5 - out[i]).abs()
        }
    }
    Ok((out, err))
}


#[pyfunction(signature = (ts, te, y0, yz =-5.0, k=0.1, extra=false, a=1.0, d=10.0, tol=0.1, hmax=1.0, hmin=1e-14))]
pub fn rkf(ts: f64, te: f64, y0: f64, yz: f64, k: f64, extra: bool, a: f64, d: f64, tol: f64, hmax:f64, hmin:f64) -> PyResult<(Vec<f64>, Vec<f64>)> {
    let mut out: Vec<f64> = vec![y0];
    let mut t: Vec<f64> = vec![ts];
    let a2 = 1.0/4.0;
    let a3 = 3.0/8.0;
    let a4 = 12.0/13.0;
    let a5 = 1.0;
    let a6 = 1.0/2.0;

    let b21 = 1.0/4.0;
    let b31 = 3.0/32.0;
    let b32 = 9.0/32.0;
    let b41 = 1932.0/2197.0;
    let b42 = -7200.0/2197.0;
    let b43 = 7296.0/2197.0;
    let b51 = 439.0/216.0;
    let b52 = -8.0;
    let b53 = 3680.0/513.0;
    let b54 = -845.0/4104.0;
    let b61 = -8.0/27.0;
    let b62 = 2.0;
    let b63 = -3544.0/2565.0;
    let b64 = 1859.0/4104.0;
    let b65 = -11.0/40.0;

    let c1 = 25.0/216.0;
    // let c2 = 0.0;
    let c3 = 1408.0/2565.0;
    let c4 = 2197.0/4104.0;
    let c5 = -1.0/5.0;

    let r1 = 1.0/360.0;
    // let r2 = 0.0;
    let r3 = -128.0/4275.0;
    let r4 = -2197.0/75240.0;
    let r5 = 1.0/50.0;
    let r6 = 2.0/55.0;

    let mut k1;
    let mut k2;
    let mut k3;
    let mut k4;
    let mut k5;
    let mut k6;
    let mut r;
    let mut tic = ts;
    let mut y = y0;
    let mut h = hmax;
    if extra {
        while tic < te {
            if (tic + h) > te {
                h = te - tic;
            }
            k1 = h * f_ext(tic, y, k, yz, a, d);
            k2 = h * f_ext(tic + a2 * h, y + b21 * k1, k, yz, a, d);
            k3 = h * f_ext(tic + a3 * h, y + b31*k1 + b32*k2, k, yz, a, d);
            k4 = h * f_ext(tic + a4 * h, y + b41*k1 + b42*k2 + b43*k3, k, yz, a, d);
            k5 = h * f_ext(tic + a5 * h, y + b51*k1 + b52*k2 + b53*k3 + b54*k4, k, yz, a, d);
            k6 = h * f_ext(tic + a6 * h, y + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5, k, yz, a, d);

            r = (( r1 * k1 + r3 * k3 + r4 * k4 + r5 * k5 + r6 * k6 ) / h).abs();

            if r < tol {
                tic += h;
                y += c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5;
                out.push(y);
                t.push(tic);
            }

            h = h * ((0.84 * ( tol / r ).powf(0.25)).max(0.1)).min(4.0);

            if h > hmax {
                h = hmax
            }
            else if h < hmin {
                return Err(PyValueError::new_err("Error, stepsize should be smaller than the minimum allowed"));
            }
        }
    }
    else {
        while tic < te {
            if (tic + h) > te {
                h = te - tic;
            }
            k1 = h * f(y, k, yz);
            k2 = h * f(y + b21 * k1, k, yz);
            k3 = h * f(y + b31*k1 + b32*k2, k, yz);
            k4 = h * f(y + b41*k1 + b42*k2 + b43*k3, k, yz);
            k5 = h * f(y + b51*k1 + b52*k2 + b53*k3 + b54*k4, k, yz);
            k6 = h * f(y + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5, k, yz);

            r = (( r1 * k1 + r3 * k3 + r4 * k4 + r5 * k5 + r6 * k6 ) / h).abs();

            if r < tol {
                tic += h;
                y += c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5;
                out.push(y);
                t.push(tic);
            }

            h = h * ((0.84 * ( tol / r ).powf(0.25)).max(0.1)).min(4.0);

            if h > hmax {
                h = hmax
            }
            else if h < hmin {
                return Err(PyValueError::new_err("Error, stepsize should be smaller than the minimum allowed"));
            }
        }
    }
    Ok((out, t))
}


#[pyfunction(signature = (t, y0, yz =-5.0, k=0.1, extra=false, a=1.0, d=10.0))]
pub fn pc4(t: Vec<f64>, y0: f64, yz: f64, k: f64, extra: bool, a: f64, d: f64) -> PyResult<Vec<f64>> {
    let mut out: Vec<f64> = vec![y0; t.len()];
    let mut h;
    let mut k1;
    let mut k2;
    let mut k3;
    let mut k4;
    let mut f0;
    let mut f1 = 0.0;
    let mut f2 = 0.0;
    let mut f3 = 0.0;
    let mut w;
    let mut fw;
    if extra {
        for i in 0..(t.len() - 1).min(3) {
            h = t[i+1] - t[i];
            f0 = f_ext(t[i], out[i], k, yz, a, d);
            k1 = h * f0;
            k2 = h * f_ext(t[i] + 0.5*h, out[i] + 0.5*k1, k, yz, a, d);
            k3 = h * f_ext(t[i] + 0.5*h, out[i] + 0.5*k2, k, yz, a, d);
            k4 = h * f_ext(t[i+1], out[i] + k3, k, yz, a, d);
            out[i+1] = out[i] + ( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0;
            (f1, f2, f3) = (f0, f1, f2);
        }
        for i in 3..(t.len() - 1) {
            h = t[i+1] - t[i];
            f0 = f_ext(t[i], out[i], k, yz, a, d);
            w = out[i] + h * (55.0 * f0 - 59.0 * f1 + 37.0 * f2 - 9.0 * f3) / 24.0;
            fw = f_ext(t[i+1], w, k, yz, a, d);
            out[i+1] = out[i] + h * (9.0 * fw + 19.0 * f0 - 5.0 * f1 + f2) / 24.0;
            (f1, f2, f3) = (f0, f1, f2);
        }
    }
    else {
        for i in 0..(t.len() - 1).min(3) {
            h = t[i+1] - t[i];
            f0 = f(out[i], k, yz);
            k1 = h * f0;
            k2 = h * f(out[i] + 0.5*k1, k, yz);
            k3 = h * f(out[i] + 0.5*k2, k, yz);
            k4 = h * f(out[i] + k3, k, yz);
            out[i+1] = out[i] + ( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0;
            (f1, f2, f3) = (f0, f1, f2);
        }
        for i in 3..(t.len() - 1) {
            h = t[i+1] - t[i];
            f0 = f(out[i], k, yz);
            w = out[i] + h * (55.0 * f0 - 59.0 * f1 + 37.0 * f2 - 9.0 * f3) / 24.0;
            fw = f(w, k, yz);
            out[i+1] = out[i] + h * (9.0 * fw + 19.0 * f0 - 5.0 * f1 + f2) / 24.0;
            (f1, f2, f3) = (f0, f1, f2);
        }
    }
    Ok(out)
}