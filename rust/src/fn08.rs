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
        // println!("{:?}",out)
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

    println!("{:2}: z = {:10.3e}, error = {:10.3e}", 0, zz1, b - w1);

    for i in 0..max_iter {
        y = rk4(f, [a, zz2], &t);
        w2 = y[n-1][0];
        println!("{:2}: z = {:10.3e}, error = {:10.3e}", i+1, zz1, b - w2);

        if (b - w2).abs() < tol {
            break;
        }

        (zz1, zz2) = ( zz2, zz2 + ( zz2 - zz1 ) / ( w2 - w1 ) * ( b - w2 ) );
        w1 = w2;
    }

    if (b - w2).abs() >= tol {
        println!("Maximum iter num {} exceeded", max_iter);
        println!("error estimate is {:10.3e}", (b - w2));
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

    println!("{:2}: z = {:10.3e}, error = {:10.3e}", 0, zz1, b - w1);

    for i in 0..max_iter {
        y = rk4(f, [a, zz2], &t);
        w2 = y[n-1][1];
        println!("{:2}: z = {:10.3e}, error = {:10.3e}", i+1, zz1, b - w2);

        if (b - w2).abs() < tol {
            break;
        }

        (zz1, zz2) = ( zz2, zz2 + ( zz2 - zz1 ) / ( w2 - w1 ) * ( b - w2 ) );
        w1 = w2;
    }

    if (b - w2).abs() >= tol {
        println!("Maximum iter num {} exceeded", max_iter);
        println!("error estimate is {:10.3e}", (b - w2));
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

    println!("{:2}: z = {:10.3e}, error = {:10.3e}", 0, zz1, b - w1);

    for i in 0..max_iter {
        y = rk4(f, [zz2, a], &t);
        w2 = y[n-1][0];
        println!("{:2}: z = {:10.3e}, error = {:10.3e}", i+1, zz1, b - w2);

        if (b - w2).abs() < tol {
            break;
        }

        (zz1, zz2) = ( zz2, zz2 + ( zz2 - zz1 ) / ( w2 - w1 ) * ( b - w2 ) );
        w1 = w2;
    }

    if (b - w2).abs() >= tol {
        println!("Maximum iter num {} exceeded", max_iter);
        println!("error estimate is {:10.3e}", (b - w2));
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

    println!("{:2}: z = {:10.3e}, error = {:10.3e}", 0, zz1, b - w1);

    for i in 0..max_iter {
        y = rk4(f, [zz2, a], &t);
        w2 = y[n-1][1];
        println!("{:2}: z = {:10.3e}, error = {:10.3e}", i+1, zz1, b - w2);

        if (b - w2).abs() < tol {
            break;
        }

        (zz1, zz2) = ( zz2, zz2 + ( zz2 - zz1 ) / ( w2 - w1 ) * ( b - w2 ) );
        w1 = w2;
    }

    if (b - w2).abs() >= tol {
        println!("Maximum iter num {} exceeded", max_iter);
        println!("error estimate is {:10.3e}", (b - w2));
        //return Err("Napakaaaaaaaaaa!!!!".into());
    }
    

    y
}


#[pyfunction]
pub fn shooter(a: f64, b: f64, z1: f64, z2: f64, t: Vec<f64>, tol: f64, e:f64, der1: bool, der2: bool) -> PyResult<Vec<[f64; 2]>> {
    let f = |y: [f64; 2], _t: f64| -> [f64; 2] {
        naloga(y, e)
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
        // println!("{:?}",out)
    }
    out
}