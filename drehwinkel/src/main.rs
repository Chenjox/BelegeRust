use std::f64::consts::PI;

use plotters::prelude::*;
use nalgebra::{SMatrix, SVector};

type Matrix3x3 = SMatrix<f64, 3, 3>;
type Vector3f = SVector<f64, 3>;

fn plot_moment_to_file(path: &str, c1: f64,
    c2: f64,
    c3: f64,
    c4: f64,
    eps: f64,
    emodul_ftm: f64,
    laenge: f64){
    let root = BitMapBackend::new(path, (640, 480)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut chart = ChartBuilder::on(&root)
        .caption("Momentenverlauf", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0f64..8f64, -200f64..100f64)
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    chart
        .draw_series(LineSeries::new(
            (0..=80)
            .map(|x| x as f64 / 80.0)
            .map(|x| (x*8.0, -durchbieg(c1, c2, c3, c4, eps, x, emodul_ftm, laenge).1))
            ,
            &RED,
        )).unwrap()
        .label("Momentenverlauf")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw().unwrap();

    root.present();
}

fn plot_biege_to_file(path: &str, c1: f64,
    c2: f64,
    c3: f64,
    c4: f64,
    eps: f64,
    emodul_ftm: f64,
    laenge: f64){
    let root = BitMapBackend::new(path, (640, 480)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut chart = ChartBuilder::on(&root)
        .caption("Biegelinie", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0f64..8f64, -0.1f64..0.1f64)
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    chart
        .draw_series(LineSeries::new(
            (0..=80)
            .map(|x| x as f64 / 80.0)
            .map(|x| (x*8.0, -durchbieg(c1, c2, c3, c4, eps, x, emodul_ftm, laenge).0))
            ,
            &RED,
        )).unwrap()
        .label("Biegelinie")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw().unwrap();

    root.present();
}


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
    plot_biege_to_file("vm001Bieg.png", sol[0], sol[1], sol[2], -sol[0], eps, emodul_ftm, laenge);
    
    plot_moment_to_file("vm001mom.png", sol[0], sol[1], sol[2], -sol[0], eps, emodul_ftm, laenge);

    let vm = 3.0e-2;
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
    plot_biege_to_file("vm003Bieg.png", sol[0], sol[1], sol[2], -sol[0], eps, emodul_ftm, laenge);
    
    plot_moment_to_file("vm003mom.png", sol[0], sol[1], sol[2], -sol[0], eps, emodul_ftm, laenge);
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
