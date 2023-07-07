use std::{error, fs::File, io::Write, ops::Mul, path::Path};

use linfa_linalg::{
    cholesky::InverseC,
    qr::{QRInto, QR},
};
use ndarray::{array, Array1, Array2};

macro_rules! float_interval {
    ($e:expr,$a:expr,$b:expr) => {
        $e >= $a && $e <= $b
    };
}

struct GebrochenRational {
    p_vector: Vec<f64>,
    q_vector: Vec<f64>,
    n: usize,
}

impl GebrochenRational {
    fn new(p_vector: Vec<f64>, q_vector: Vec<f64>) -> Self {
        if p_vector.len() - 1 != q_vector.len() {
            panic!(
                "Lenghts of given Vectors do not match! {} != {}",
                p_vector.len(),
                q_vector.len()
            );
        }
        let mut q_vector = q_vector;
        let n = p_vector.len();
        Self {
            p_vector,
            q_vector,
            n,
        }
    }

    fn get_linear_isomorphism(&self) -> (Array2<f64>, Array2<f64>) {
        let mut p = self.p_vector.clone();
        let mut q = self.q_vector.clone();
        q.push(0.0);
        let n = self.n;
        let mut A: Array2<f64> = Array2::zeros([n - 1, n - 1]);
        let mut B: Array2<f64> = Array2::zeros([n - 1, n - 1]);

        println!("{:?},\n{:?}", p, q);
        let mut r_2 = p;
        let mut r_1 = q;
        // FIXME: Funktioniert noch nicht.
        let mut j = 0;
        loop {
            // Abbruchbedingung
            if j == n - 1 {
                break;
            }
            let s_1 = r_2[n - j - 1] / r_1[n - j - 2];
            // Hier wird im letzten Schitt gepanict!
            let s_0 = if j == n - 2 {
                (r_2[n - j - 2]) / r_1[n - j - 2]
            } else {
                (r_2[n - j - 2] - s_1 * r_1[n - j - 3]) / r_1[n - j - 2]
            };
            // r_j
            let mut r_0 = vec![0.0; n - j - 1];
            if j < n - 3 {
                // Wir zählen von unten nach oben!
                let mut m = n - j - 3;
                loop {
                    if m == 0 {
                        break;
                    }
                    r_0[m] = r_2[m] - s_0 * r_1[m] - s_1 * r_1[m - 1];
                    m -= 1;
                }
            }
            r_0[0] = r_2[0] - s_0 * r_1[0];
            println!("{:?}", r_0);

            A[[j, j]] = s_1 * (-1.0_f64).powi(j as i32);
            B[[j, j]] = s_0 * (-1.0_f64).powi(j as i32);
            if j < n - 2 {
                B[[j + 1, j]] = (-1.0_f64).powi(j as i32);
                B[[j, j + 1]] = (-1.0_f64).powi(j as i32);
            }
            // Anpassen der Vectoren
            // r_0 liegt eins in der vergangenheit
            r_2 = r_1;
            r_1 = r_0;

            // Inkrementiere um 1
            j += 1;
        }
        return (A, B);
    }
}

fn belast_function(time: f64) -> f64 {
    match time {
        dt if float_interval![dt, 0.0e-2, 1.0e-2] => (dt / (1.0e-2)) * 100.0 * 1000.0,
        dt if float_interval![dt, 1.0e-2, 2.0e-2] => 100.0 * 1000.0,
        dt if float_interval![dt, 2.0e-2, 3.0e-2] => -((dt - 3.0e-2) / (1.0e-2)) * 100.0 * 1000.0,
        _ => 0.0,
    }
}

fn write_string_to_file(nam: &str, l: Vec<[f64; 2]>) {
    let r = l
        .iter()
        .map(|f| format!("{} {}", f[0], f[1]))
        .collect::<Vec<String>>()
        .join("\n");
    let path = Path::new(&nam);
    let mut file = match File::create(&path) {
        Err(why) => {
            //error!("Couldn't create {}: {}", nam, why);
            return;
        }
        Ok(file) => file,
    };
    file.write_all(r.as_bytes()).unwrap();
    //match file.write_all(r.as_bytes()) {
    //    Err(why) => error!("couldn't write to {}: {}", nam, why),
    //    Ok(_) => info!("successfully wrote to {}", nam),
    //}
}

fn main() {
    let g = GebrochenRational::new(
        vec![
            10337706.5,
            48744.172,
            455.205728,
            1.13950676,
            0.00483245447,
            3.11054119e-6,
        ],
        vec![
            1.0,
            0.00490648418,
            2.71147593e-5,
            5.42912842e-8,
            6.94958499e-12,
        ],
    );
    let (a, b) = g.get_linear_isomorphism();

    let timestep = 1e-5;
    let mut time = 0.0;
    let max_time = 1.5;
    let v_0 = 0.000;
    let mut z_i = Array1::zeros(a.shape()[0]);
    z_i[0] = v_0;

    let mut hist = Vec::with_capacity((max_time / timestep) as usize);
    hist.push([0.0, v_0]);

    let timestep_mat_1 = &a + (&b * timestep / 2.0);
    let inv = timestep_mat_1.clone().qr_into().unwrap().inverse().unwrap();

    let timestep_mat_2 = &b * timestep / 2.0 - &a;
    //println!("{}",inv);
    //println!("{}",inv.dot(&timestep_mat));
    println!("{},\n{}", a, b);
    loop {
        let mut z_ip1 = -timestep_mat_2.dot(&z_i);
        let integral =
            ((belast_function(time + timestep) + belast_function(time)) / 2.0) * timestep;

        z_ip1[0] = z_ip1[0] + integral;
        z_i = inv.dot(&z_ip1);
        time += timestep;
        //
        hist.push([time, z_i[0]]);

        if time >= max_time {
            break;
        }

        //println!("{},{}",current[0],belast_function(time) );
    }
    write_string_to_file("result.csv", hist);
    //
}
