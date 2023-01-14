use std::f64::consts::PI;

use nalgebra::{SMatrix, SVector};

type Matrix3x3 = SMatrix<f64, 3, 3>;
type Vector3f = SVector<f64, 3>;
pub fn main() {
    let emodul_ftm: f64 = 6.0e3;
    let kphi = 750.0;
    let laenge = 8.0;
    let kphi = kphi * laenge / emodul_ftm;
    let s: f64 = 0.7 * 1087.32;
    let eps: f64 = laenge * (s / emodul_ftm).sqrt();

    println!("{},{},{}", kphi, s, eps);
    let b = Matrix3x3::new(
        eps.cos() - 1.0,
        eps.sin(),
        1.0,
        eps.powi(2),
        eps * kphi,
        kphi,
        eps.powi(2) * eps.cos(),
        eps.powi(2) * eps.sin(),
        0.0,
    );

    let vm = 1.0e-2;
    let b2 = (PI * vm * eps.powi(2)) / (PI * PI - eps.powi(2)) * kphi;
    let inhom = Vector3f::new(0.0, b2, 0.0);
    println!("{}", b2);

    let full_piv = b.full_piv_lu();
    let sol = full_piv.solve(&inhom).unwrap();
    println!("{}", b);
    println!("{}", sol);

    let xi = 1.0;
    let (gunt, mom) = durchbieg(sol[0], sol[1], sol[2], -sol[0], eps, xi, emodul_ftm, laenge);
    println!("{},{}", gunt, mom);
}

fn durchbieg(
    c1: f64,
    c2: f64,
    c3: f64,
    c4: f64,
    eps: f64,
    xi: f64,
    emodul_ftm: f64,
    laenge: f64,
) -> (f64, f64) {
    let gunt = c1 * (eps * xi).cos() + c2 * (eps * xi).sin() + c3 * xi + c4;
    let mom = -emodul_ftm / laenge
        * (-c1 * eps.powi(2) * (eps * xi).cos() - c2 * eps.powi(2) * (eps * xi).sin());
    return (gunt, mom);
}
