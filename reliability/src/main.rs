use std::{f64::consts::PI, fs::File};

use compute::distributions::{Continuous, Gumbel, Normal};
use faer::{mat, Mat};

// Lognormal Distribution
struct LogNormal {
  mean_nn: f64,
  sigma_nn: f64,
  lower_value: f64
}

impl LogNormal {
  fn new(mean:f64, sigma:f64, lower_bound: f64) -> Self {
    let sigma_nn = (1.0 + (sigma/(mean - lower_bound)).powi(2)).ln().sqrt();
    let mean_nn = (mean - lower_bound).ln() - sigma_nn.powi(2)/2.0;
  
    Self {
      mean_nn,
      sigma_nn,
      lower_value: lower_bound
    }
  }

  fn cdf(&self,x: f64) -> f64 {
    if x <= self.lower_value {
      return 0.0;
    }
    let n = Normal::new(self.mean_nn,self.sigma_nn);

    let trans_x = (x - self.lower_value).ln();

    return n.cdf(trans_x)
  }
}

fn double_exponential_transformation(c: f64,x: f64) -> [f64;2] {
  let trans_x = (c*(x).sinh()).sinh();
  let trans_x_deriv = (c*(x).sinh()).cosh() * c * x.cosh();

  return [trans_x,trans_x_deriv]
}
fn inverse_double_exponential_transformation(c:f64, x: f64) -> f64 {
  todo!()
}

trait IntegrableFunction {
  fn function_value(&self,x: f64) -> f64;
}

struct Task1Function {
  load_mat: Mat<f64>,
}

const simpson_weights: [f64;3] = [1.0/6.0,4.0/6.0,1.0/6.0];
const simpson_places: [f64;3] = [0.0, 0.5, 1.0];

fn integrate_function<T: IntegrableFunction>(func: &T, start_point: f64, step_size: f64, epsilon: f64) -> Option<f64> {
  // trapezoidal rule
  
  //transform the start point
  let start = double_exponential_transformation(1.0, start_point);
  // check it's NaN ness
  if !start[0].is_finite() && !start[1].is_finite() {
    return None;
  }

  // get the first interval
  let mut lower_val = start_point;
  let mut upper_val = start_point + step_size;
  
  let mut large_sum = 0.0;
  loop {
    let mut small_sum = 0.0;
    for i in 0..3 {
      let x_val = double_exponential_transformation(1.0, lower_val + step_size * simpson_places[i]);
      let func_val = step_size*func.function_value(x_val[0])*x_val[1]*simpson_weights[i];
      if !func_val.is_finite() {
        break;
      }
      small_sum += func_val
    }
    if small_sum.abs() < epsilon {
      break
    }
    //println!("{},{},{}",small_sum,lower_val,upper_val);
    large_sum += small_sum;
    lower_val = upper_val;
    upper_val += step_size;
  }
  // now everything behind
  let mut lower_val = start_point - step_size;
  let mut upper_val = start_point;
  
  loop {
    let mut small_sum = 0.0;
    for i in 0..3 {
      let x_val = double_exponential_transformation(1.0, upper_val - step_size * simpson_places[i]);
      let func_val = step_size*func.function_value(x_val[0])*x_val[1]*simpson_weights[i];
      if !func_val.is_finite() {
        break;
      }
      small_sum += func_val
    }
    if small_sum.abs() < epsilon {
      break
    }
    //println!("{},{},{}",small_sum,lower_val,upper_val);
    large_sum += small_sum;
    upper_val = lower_val;
    lower_val -= step_size;
  }
  
  return Some(large_sum)
}

struct TestFun {

}

impl IntegrableFunction for TestFun {
  fn function_value(&self,x: f64) -> f64 {
      1.0/(2.0 * std::f64::consts::PI).sqrt() *(-0.5*x.powi(2)).exp()
  }
}

struct Task1 {
  load_vec: Vec<f64>,
  area_vec: Vec<f64>
}

impl IntegrableFunction for Task1 {
  fn function_value(&self,x: f64) -> f64 {
    let fy = LogNormal::new(30.20e4,1.44e4,19.9e4);
    let load = Gumbel::new(410.0,70.0);
    
    let mut prod = 1.0;
    for i in 0..self.load_vec.len() {
      prod *= 1.0 - fy.cdf(self.load_vec[i]/self.area_vec[i]*x)
    }

    return prod*load.pdf(x);
  }
}

fn main(){

  let task1 = Task1 {
    load_vec: vec![],
    area_vec: vec![]
  };

  let test = TestFun{};

  if let Some(result) = integrate_function(&test, 0.0, 0.01, 1e-15) {
    println!("{}",result);
  };

  //println!("{}",load.pdf(60.0));
  //println!("{}",fy.cdf(50.0e4));
}