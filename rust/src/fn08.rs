use pyo3::prelude::*;
// use pyo3::types::PyAny;
// use std::error::Error;


// fn f(x: f64) -> f64 {
//     x
// }


// fn call(func, x: f64) -> f64 {
//     func(x)
// }

fn naloga(y: [f64; 2], e: f64) -> [f64; 2] {
    [y[1], -e * y[0]]
}

fn nalogica(y: [f64; 2], e: f64, t: f64) -> [f64; 2] {
    if t < 0.0 || t > 1.0 {
        [y[1], -(e - 100.0) * y[0]]
    }
    else {
        [y[1], -e * y[0]]
    }
}

#[pyfunction]
pub fn schrodinger(y: [f64; 2], e: f64) -> [f64; 2] {
    naloga(y, e)
}

fn rk4<F>(f: F, y0: [f64; 2], t: &Vec<f64>) -> Vec<[f64; 2]>
where
    F: Fn([f64; 2], f64) -> [f64; 2],
    {
    let mut out: Vec<[f64; 2]> = vec![y0; t.len()];
    let mut h;
    let mut k1;
    let mut k2;
    let mut k3;
    let mut k4;
    for i in 0..(t.len() - 1) {
        h = t[i+1] - t[i];
        k1 = f(out[i], t[i]);
        k1 = [h * k1[0], h * k1[1]];
        k2 = f([out[i][0] + 0.5 * k1[0], out[i][1] + 0.5 * k1[1]], t[i] + 0.5 * h);
        k2 = [h * k2[0], h * k2[1]];
        k3 = f([out[i][0] + 0.5 * k2[0], out[i][1] + 0.5 * k2[1]], t[i] + 0.5 * h);
        k3 = [h * k3[0], h * k3[1]];
        k4 = f([out[i][0] + k3[0], out[i][1] + k3[1]], t[i+1]);
        k4 = [h * k4[0], h * k4[1]];
        out[i+1] = [out[i][0] + ( k1[0] + 2.0 * ( k2[0] + k3[0] ) + k4[0] ) / 6.0, out[i][1] + ( k1[1] + 2.0 * ( k2[1] + k3[1] ) + k4[1] ) / 6.0];
        // // println!("{:?}",out)
    }
    out
}

// lahko bi bla ena funcija z dvema if stavkoma, ampak je bil development time tkole krajsi
fn shoot_valval<F>(f: &F, a: f64, b: f64, z1: f64, z2: f64, t: Vec<f64>, tol: f64) -> Vec<[f64; 2]>
where
    F: Fn([f64; 2], f64) -> [f64; 2],
    {
    let n = t.len();
    let max_iter: u16 = 30;
    let mut y: Vec<[f64; 2]> = rk4(f, [a, z1], &t);
    let mut w1 = y[n-1][0];
    let mut w2 = w1;
    let mut zz1 = z1;
    let mut zz2 = z2;

    // // println!("{:2}: z = {:10.3e}, error = {:10.3e}", 0, zz1, b - w1);

    for i in 0..max_iter {
        y = rk4(f, [a, zz2], &t);
        w2 = y[n-1][0];
        // println!("{:2}: z = {:10.3e}, error = {:10.3e}", i+1, zz1, b - w2);

        if (b - w2).abs() < tol {
            break;
        }

        (zz1, zz2) = ( zz2, zz2 + ( zz2 - zz1 ) / ( w2 - w1 ) * ( b - w2 ) );
        w1 = w2;
    }

    if (b - w2).abs() >= tol {
        // println!("Maximum iter num {} exceeded", max_iter);
        // println!("error estimate is {:10.3e}", (b - w2));
        //return Err("Napakaaaaaaaaaa!!!!".into());
    }
    

    y
}


fn shoot_valder<F>(f: &F, a: f64, b: f64, z1: f64, z2: f64, t: Vec<f64>, tol: f64) -> Vec<[f64; 2]>
where
    F: Fn([f64; 2], f64) -> [f64; 2],
    {
    let n = t.len();
    let max_iter: u16 = 30;
    let mut y: Vec<[f64; 2]> = rk4(f, [a, z1], &t);
    let mut w1 = y[n-1][1];
    let mut w2 = w1;
    let mut zz1 = z1;
    let mut zz2 = z2;

    // println!("{:2}: z = {:10.3e}, error = {:10.3e}", 0, zz1, b - w1);

    for i in 0..max_iter {
        y = rk4(f, [a, zz2], &t);
        w2 = y[n-1][1];
        // println!("{:2}: z = {:10.3e}, error = {:10.3e}", i+1, zz1, b - w2);

        if (b - w2).abs() < tol {
            break;
        }

        (zz1, zz2) = ( zz2, zz2 + ( zz2 - zz1 ) / ( w2 - w1 ) * ( b - w2 ) );
        w1 = w2;
    }

    if (b - w2).abs() >= tol {
        // println!("Maximum iter num {} exceeded", max_iter);
        // println!("error estimate is {:10.3e}", (b - w2));
        //return Err("Napakaaaaaaaaaa!!!!".into());
    }
    

    y
}


fn shoot_derval<F>(f: &F, a: f64, b: f64, z1: f64, z2: f64, t: Vec<f64>, tol: f64) -> Vec<[f64; 2]>
where
    F: Fn([f64; 2], f64) -> [f64; 2],
    {
    let n = t.len();
    let max_iter: u16 = 30;
    let mut y: Vec<[f64; 2]> = rk4(f, [z1, a], &t);
    let mut w1 = y[n-1][0];
    let mut w2 = w1;
    let mut zz1 = z1;
    let mut zz2 = z2;

    // println!("{:2}: z = {:10.3e}, error = {:10.3e}", 0, zz1, b - w1);

    for i in 0..max_iter {
        y = rk4(f, [zz2, a], &t);
        w2 = y[n-1][0];
        // println!("{:2}: z = {:10.3e}, error = {:10.3e}", i+1, zz1, b - w2);

        if (b - w2).abs() < tol {
            break;
        }

        (zz1, zz2) = ( zz2, zz2 + ( zz2 - zz1 ) / ( w2 - w1 ) * ( b - w2 ) );
        w1 = w2;
    }

    if (b - w2).abs() >= tol {
        // println!("Maximum iter num {} exceeded", max_iter);
        // println!("error estimate is {:10.3e}", (b - w2));
        //return Err("Napakaaaaaaaaaa!!!!".into());
    }
    

    y
}


fn shoot_derder<F>(f: &F, a: f64, b: f64, z1: f64, z2: f64, t: Vec<f64>, tol: f64) -> Vec<[f64; 2]>
where
    F: Fn([f64; 2], f64) -> [f64; 2],
    {
    let n = t.len();
    let max_iter: u16 = 30;
    let mut y: Vec<[f64; 2]> = rk4(f, [z1, a], &t);
    let mut w1 = y[n-1][1];
    let mut w2 = w1;
    let mut zz1 = z1;
    let mut zz2 = z2;

    // println!("{:2}: z = {:10.3e}, error = {:10.3e}", 0, zz1, b - w1);

    for i in 0..max_iter {
        y = rk4(f, [zz2, a], &t);
        w2 = y[n-1][1];
        // println!("{:2}: z = {:10.3e}, error = {:10.3e}", i+1, zz1, b - w2);

        if (b - w2).abs() < tol {
            break;
        }

        (zz1, zz2) = ( zz2, zz2 + ( zz2 - zz1 ) / ( w2 - w1 ) * ( b - w2 ) );
        w1 = w2;
    }

    if (b - w2).abs() >= tol {
        // println!("Maximum iter num {} exceeded", max_iter);
        // println!("error estimate is {:10.3e}", (b - w2));
        //return Err("Napakaaaaaaaaaa!!!!".into());
    }
    

    y
}


#[pyfunction(signature = (a,b,z1,z2,t,tol,e,der1,der2,extra=false))]
pub fn shooter(a: f64, b: f64, z1: f64, z2: f64, t: Vec<f64>, tol: f64, e:f64, der1: bool, der2: bool, extra: bool) -> PyResult<Vec<[f64; 2]>> {
    let f: Box<dyn Fn([f64; 2], f64) -> [f64; 2]> = if extra {
        Box::new(move |y: [f64; 2], t: f64| {
            nalogica(y, e, t)
        })
    } else {
        Box::new(move |y: [f64; 2], _t: f64| {
            naloga(y, e)
        })
    };
    if der1 && der2 {
        Ok(shoot_derder(&f, a, b, z1, z2, t, tol))
    }
    else if der1 && !der2 {
        Ok(shoot_derval(&f, a, b, z1, z2, t, tol))
    }
    else if !der1 && der2 {
        Ok(shoot_valder(&f, a, b, z1, z2, t, tol))
    }
    else {
        Ok(shoot_valval(&f, a, b, z1, z2, t, tol))
    }
}


#[pyfunction]
pub fn sch_rk4(e: f64, y0: [f64; 2], t: Vec<f64>) -> Vec<[f64; 2]> {
    let f = |y: [f64; 2], _t: f64| -> [f64; 2] {
        naloga(y, e)
    };
    let mut out: Vec<[f64; 2]> = vec![y0; t.len()];
    let mut h;
    let mut k1;
    let mut k2;
    let mut k3;
    let mut k4;
    for i in 0..(t.len() - 1) {
        h = t[i+1] - t[i];
        k1 = f(out[i], t[i]);
        k1 = [h * k1[0], h * k1[1]];
        k2 = f([out[i][0] + 0.5 * k1[0], out[i][1] + 0.5 * k1[1]], t[i] + 0.5 * h);
        k2 = [h * k2[0], h * k2[1]];
        k3 = f([out[i][0] + 0.5 * k2[0], out[i][1] + 0.5 * k2[1]], t[i] + 0.5 * h);
        k3 = [h * k3[0], h * k3[1]];
        k4 = f([out[i][0] + k3[0], out[i][1] + k3[1]], t[i+1]);
        k4 = [h * k4[0], h * k4[1]];
        out[i+1] = [out[i][0] + ( k1[0] + 2.0 * ( k2[0] + k3[0] ) + k4[0] ) / 6.0, out[i][1] + ( k1[1] + 2.0 * ( k2[1] + k3[1] ) + k4[1] ) / 6.0];
        // // println!("{:?}",out)
    }
    out
}


#[pyfunction]
pub fn fd(u: Vec<f64>, v:Vec<f64>, w:Vec<f64>, t: Vec<f64>, z: f64, k: f64) -> PyResult<Vec<f64>> {
    let n = t.len();
    let h = t[0] - t[1];
    let mut a: Vec<f64> = w[1..n]
        .iter()
        .map(|&wi| -(1.0 + wi * h / 2.0))
        .collect();

    // Set last element to 0.0
    // if let Some(last) = a.last_mut() {
    //     *last = 0.0;
    // };
    a[n-2] = 0.0;

    let mut c: Vec<f64> = w[0..(n-1)]
        .iter()
        .map(|&wi| -(1.0 + wi * h / 2.0))
        .collect();

    // Set last element to 0.0
    // if let Some(last) = c.last_mut() {
    //     *last = 0.0;
    // };
    c[0] = 0.0;

    let mut d: Vec<f64> = v
        .iter()
        .map(|&vi| 2.0 + vi * h *h )
        .collect();

    d[0] = 1.0;
    d[n-1] = 1.0;

    let mut b: Vec<f64> = u
        .iter()
        .map(|&ui| - ui * h *h )
        .collect();

    b[0] = z;
    b[n-1] = k;
    let mut xmult: f64;
    let mut x = vec![0.0; n];

    for i in 1..n {
        xmult = a[i-1] / d[i-1];
        d[i] = d[i] - xmult * c[i-1];
        b[i] = b[i] - xmult * b[i-1];
    }
    x[n-1] = b[n-1] / d[n-1];
    
    for i in (0..=n - 2).rev() {
        x[i] = (b[i] - c[i] * x[i + 1]) / d[i];
    }

    Ok(x)

}


#[pyfunction]
pub fn fdd(
    d_extra: Vec<f64>,
    e_extra: Vec<f64>,
) -> PyResult<(Vec<f64>, Vec<Vec<f64>>)> {
    let mut d = d_extra;
    let mut e = e_extra;
    e.push(0.0);
    let n = d.len();
    assert_eq!(e.len(), n);
    

    // Initialize eigenvector matrix as identity
    let mut z = vec![vec![0.0; n]; n];
    for i in 0..n {
        z[i][i] = 1.0;
    }

    // Extend e with a zero (algorithm convenience)
    // let mut e = {
    //     let mut tmp = e.to_vec();
    //     tmp.push(0.0);
    //     tmp
    // };

    for l in 0..n {
        let mut iter: i32 = 0;
        loop {
            // Find small subdiagonal element
            let mut m = l;
            while m + 1 < n && e[m].abs() > 1e-15 {
                m += 1;
            }

            if m == l {
                break;
            }

            if iter > 100 {
                panic!("QL algorithm failed to converge");
            }
            iter += 1;

            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let r = (g * g + 1.0).sqrt();
            g = d[m] - d[l] + e[l] / (g + r.copysign(g));

            let mut s = 1.0;
            let mut c = 1.0;
            let mut p = 0.0;

            for i in (l..m).rev() {
                let f = s * e[i];
                let b = c * e[i];

                if f.abs() >= g.abs() {
                    c = g / f;
                    let r = (c * c + 1.0).sqrt();
                    e[i + 1] = f * r;
                    s = 1.0 / r;
                    c *= s;
                } else {
                    s = f / g;
                    let r = (s * s + 1.0).sqrt();
                    e[i + 1] = g * r;
                    c = 1.0 / r;
                    s *= c;
                }

                let g2 = d[i + 1] - p;
                let r2 = (d[i] - g2) * s + 2.0 * c * b;
                p = s * r2;
                d[i + 1] = g2 + p;
                g = c * r2 - b;

                // Apply rotation to eigenvectors
                for k in 0..n {
                    let t = z[k][i + 1];
                    z[k][i + 1] = s * z[k][i] + c * t;
                    z[k][i] = c * z[k][i] - s * t;
                }
            }

            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }
    }

    // Sort eigenvalues and eigenvectors
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&i, &j| d[i].partial_cmp(&d[j]).unwrap());

    let eigenvalues: Vec<f64> = idx.iter().map(|&i| d[i]).collect();

    let mut eigenvectors = vec![vec![0.0; n]; n];
    for (new_i, &old_i) in idx.iter().enumerate() {
        for k in 0..n {
            eigenvectors[new_i][k] = z[k][old_i];
        }
    }

    Ok((eigenvalues, eigenvectors))
}


#[pyfunction]
pub fn ffd(
    mut d: Vec<f64>,
    mut e: Vec<f64>,
) -> PyResult<(Vec<f64>, Vec<Vec<f64>>)> {
    let n = d.len();
    assert_eq!(e.len(), n.saturating_sub(1));

    // Initialize eigenvector matrix as identity
    let mut q: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();

    let eps = f64::EPSILON;
    let max_it = 1000;

    for l in (0..n).rev() {
        let mut iter = 0;
        while iter < max_it {
            // check convergence of subdiagonal
            let mut m = l;
            while m > 0 && e[m - 1].abs() > eps * (d[m - 1].abs() + d[m].abs()) {
                m -= 1;
            }
            if m == l {
                break; // eigenvalue isolated
            }

            // compute Wilkinson shift
            let dd = (d[l - 1] - d[l]) / 2.0;
            let mu = d[l]
                - (e[l - 1].powi(2))
                    / (dd + dd.signum() * (dd.powi(2) + e[l - 1].powi(2)).sqrt());

            let mut x = d[m] - mu;
            let mut z = e[m];

            // perform QR sweep
            let mut c;
            let mut s;
            for k in m..=l - 1 {
                let r = (x.powi(2) + z.powi(2)).sqrt();
                if r == 0.0 {
                    c = 1.0;
                    s = 0.0;
                } else {
                    c = x / r;
                    s = z / r;
                }

                if k > m {
                    e[k - 1] = r;
                }

                let dk  = d[k];
                let ek  = e[k];
                let dk1 = d[k + 1];

                // apply Givens on [ d[k], e[k] ]
                let y  = c * dk + s * ek;
                let t  = -s * dk + c * ek;

                // apply Givens on [ t, d[k+1] ]
                let u  = c * t + s * dk1;
                let v  = -s * t + c * dk1;

                d[k]     = y;
                e[k]     = u;     // now using the intermediate
                d[k + 1] = v;

                for i in 0..n {
                    let tmp = q[i][k] * c + q[i][k + 1] * s;
                    q[i][k + 1] = -q[i][k] * s + q[i][k + 1] * c;
                    q[i][k] = tmp;
                }

                if k < l - 1 {
                    x = e[k];
                    z = -s * e[k + 1];
                    e[k + 1] *= c;
                }
            }
            e[l - 1] = x;
            d[l] -= mu;
            d[m] += mu;

            iter += 1;
        }
    }

    // Sort eigenvalues with eigenvectors
    let mut eig_pairs: Vec<(f64, Vec<f64>)> = d
        .into_iter()
        .zip((0..n).map(|j| (0..n).map(|i| q[i][j]).collect()))
        .collect();

    eig_pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let (values, vectors): (Vec<_>, Vec<_>) = eig_pairs.into_iter().unzip();
    Ok((values, vectors))
}