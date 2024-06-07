use std::{f64::consts::PI, fs::File};

use compensated_summation::KahanBabuskaNeumaier;
use compute::{distributions::{Continuous, Gumbel, Normal}, integrate::trapz};
use faer::{mat, Mat};

const EULER_MASCHERONI: f64 = 0.577215664901532860606512090082402431042159335939923598805767234884867726777664670936947063291746749;
// Lognormal Distribution
struct LogNormal {
  mean_nn: f64,
  sigma_nn: f64,
  lower_value: f64,
}

impl LogNormal {
  fn new(mean: f64, sigma: f64, lower_bound: f64) -> Self {
    let sigma_nn = (1.0 + (sigma / (mean - lower_bound)).powi(2)).ln().sqrt();
    let mean_nn = (mean - lower_bound).ln() - sigma_nn.powi(2) / 2.0;

    Self {
      mean_nn,
      sigma_nn,
      lower_value: lower_bound,
    }
  }

  fn cdf(&self, x: f64) -> f64 {
    if x <= self.lower_value {
      return 0.0;
    }
    let n = Normal::new(self.mean_nn, self.sigma_nn);

    let trans_x = (x - self.lower_value).ln();

    return n.cdf(trans_x);
  }
}

fn double_exponential_transformation(c: f64, x: f64) -> [f64; 2] {
  let trans_x = (c * (x).sinh()).sinh();
  let trans_x_deriv = (c * (x).sinh()).cosh() * c * x.cosh();

  return [trans_x, trans_x_deriv];
}
fn inverse_double_exponential_transformation(c: f64, x: f64) -> f64 {
  (c * (x.asinh())).asinh()
}

trait IntegrableFunction {
  fn function_value(&self, x: f64) -> f64;
}

struct Task1Function {
  load_mat: Mat<f64>,
}

const simpson_weights: [f64; 3] = [1.0 / 6.0, 4.0 / 6.0, 1.0 / 6.0];
const simpson_places: [f64; 3] = [0.0, 0.5, 1.0];

fn integrate_function<T: IntegrableFunction>(
  func: &T,
  start_point: f64,
  step_size: f64,
  epsilon: f64,
) -> Option<f64> {
  // trapezoidal rule

  //transform the start point
  let start = inverse_double_exponential_transformation(1.0, start_point);
  let inkre = inverse_double_exponential_transformation(1.0, start_point+step_size);
  let step_size = (inkre - start).abs();
  //let start = double_exponential_transformation(1.0, start_point);
  // check it's NaN ness
  //if !start[0].is_finite() && !start[1].is_finite() {
  //  println!("{},{}",start_point,start[0]);
  //  return None;
  //}

  // get the first interval
  let mut lower_val = start;
  let mut upper_val = start + step_size;

  let mut sum = KahanBabuskaNeumaier::new();
  loop {
    let mut small_sum = 0.0;
    for i in 0..3 {
      let x_val = double_exponential_transformation(1.0, lower_val + step_size * simpson_places[i]);
      let func_val = step_size * func.function_value(x_val[0]) * x_val[1] * simpson_weights[i];
      if !func_val.is_finite() {
        break;
      }
      small_sum += func_val
    }
    if small_sum.abs() < epsilon {
      break;
    }
    //println!("{},{},{}",small_sum,lower_val,upper_val);
    sum += small_sum;
    lower_val = upper_val;
    upper_val += step_size;
  }
  // now everything behind
  let mut lower_val = start - step_size;
  let mut upper_val = start;

  loop {
    let mut small_sum = 0.0;
    for i in 0..3 {
      let x_val = double_exponential_transformation(1.0, upper_val - step_size * simpson_places[i]);
      let func_val = step_size * func.function_value(x_val[0]) * x_val[1] * simpson_weights[i];
      if !func_val.is_finite() {
        break;
      }
      small_sum += func_val
    }
    if small_sum.abs() < epsilon {
      break;
    }
    //println!("{},{},{}",small_sum,lower_val,upper_val);
    sum += small_sum;
    upper_val = lower_val;
    lower_val -= step_size;
  }

  return Some(sum.total());
}

struct TestFun {}

impl IntegrableFunction for TestFun {
  fn function_value(&self, x: f64) -> f64 {
    1.0 / (2.0 * std::f64::consts::PI).sqrt() * (-0.5 * x.powi(2)).exp()
  }
}

struct Task1 {
  load_vec: Vec<f64>,
  area_vec: Vec<f64>,
}

impl IntegrableFunction for Task1 {
  fn function_value(&self, x: f64) -> f64 {
    let fy = LogNormal::new(30.2e4, 1.44e4, 19.9e4);

    let mean = 410.0;
    let std_dev = 70.0;

    let beta = std_dev*6.0_f64.sqrt()/std::f64::consts::PI;
    let mu = mean - beta*EULER_MASCHERONI;

    let load = Gumbel::new(mu, beta);

    let mut prod = 1.0;
    for i in 0..self.load_vec.len() {
      prod *= 1.0 - fy.cdf((self.load_vec[i] / self.area_vec[i]).abs() * x)
    }

    return prod * load.pdf(x);
  }
}

impl Task1 {
  fn erg(&self) -> Option<f64> {
    return integrate_function(self, 410.0, 0.1, 1e-15).map(|f| 1.0 - f);
  }
}

struct Task2 {
  load_vec: Vec<f64>,
  area_vec: Vec<f64>,
  i: usize
}

impl IntegrableFunction for Task2 {
  fn function_value(&self, x: f64) -> f64 {
    let fy = LogNormal::new(30.2e4, 1.44e4, 19.9e4);

    let mean = 410.0;
    let std_dev = 70.0;

    let beta = std_dev*6.0_f64.sqrt()/std::f64::consts::PI;
    let mu = mean - beta*EULER_MASCHERONI;

    let load = Gumbel::new(mu, beta);

    return (1.0 - fy.cdf((self.load_vec[self.i] / self.area_vec[self.i]).abs() * x) )* load.pdf(x);
  }
}

impl Task2 {
  fn erg() -> f64 {
    let mut res: f64 = 0.;
    let mut k = 0;
    for i in 0..10 {
      let taks2 = Task2 { load_vec: vec![
        0.0,
        1.0,
        0.0,
        5.0 / 4.0,
        -3.0 / 4.0,
        -2.0_f64.sqrt(),
        -1.0,
        3.0 / 4.0,
        0.0,
        -7.0 / 4.0
      ],
      area_vec: vec![
        3.77e-3, 3.77e-3, 3.77e-3, 3.77e-3, 3.77e-3, 4.7e-3, 3.77e-3, 3.77e-3, 3.77e-3, 5.74e-3
      ],
      i};
      if let Some( r) = integrate_function(&taks2, 410.0, 0.01, 1e-15).map(|f| 1.0 - f) {
        if r > res {
          res = r;
          k = i;
        }
      };
    }
    println!("{}",k);
    return res;
  }
}

struct Task3 {
  load_vec: Vec<f64>,
  area: Vec<f64>,
  length: Vec<f64>,
  ftm: Vec<f64>,
  youngs_modulus: Vec<f64>
}

struct Task3Buckling {
  load_factor: f64,
  youngs_modulus: f64,
  ftm: f64,
  length: f64,
  buckling_length: f64
}

impl IntegrableFunction for Task3Buckling {
  fn function_value(&self, x: f64) -> f64 {
    
    let mean = 410.0;
    let std_dev = 70.0;

    let beta = std_dev*6.0_f64.sqrt()/std::f64::consts::PI;
    let mu = mean - beta*EULER_MASCHERONI;

    let load = Gumbel::new(mu, beta);

    let critical_load = self.youngs_modulus * self.ftm* std::f64::consts::PI.powi(2) / (self.length * self.buckling_length).powi(2);
    let point = critical_load/self.load_factor;

    return load.cdf
  }
}

fn main() {
  let task1 = Task1 {
    load_vec: vec![
      0.0,
      1.0,
      0.0,
      5.0 / 4.0,
      -3.0 / 4.0,
      -2.0_f64.sqrt(),
      -1.0,
      3.0 / 4.0,
      0.0,
      -7.0 / 4.0
    ],
    area_vec: vec![
      3.77e-3, 3.77e-3, 3.77e-3, 3.77e-3, 3.77e-3, 4.7e-3, 3.77e-3, 3.77e-3, 3.77e-3, 5.74e-3
    ],
  };


  //task1.test1();
  if let Some(erg) = task1.erg() {
    println!("{}", erg);
  }

  println!("{}", Task2::erg());

  

  //println!("{}",load.pdf(60.0));
  //println!("{}",fy.cdf(50.0e4));
}
