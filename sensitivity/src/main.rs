fn rosenbrock(input_vec: &Vec<f64>) -> f64 {
  let mut running_sum = 0.0;
  for i in 0..input_vec.len() - 1 {
    running_sum +=
      100.0 * (input_vec[i + 1] - input_vec[i].powi(2)).powi(2) + (1.0 - input_vec[i]).powi(2);
  }
  return running_sum;
}

fn rosenbrock_gradient(input_vec: &Vec<f64>) -> Vec<f64> {
  let mut result = vec![0.0; input_vec.len()];
  for i in 0..input_vec.len() {
    if i == 0 {
      result[i] = 400.0 * input_vec[0].powi(3) - 400.0 * input_vec[0] * input_vec[1]
        + 2.0 * input_vec[0]
        - 2.0;
    } else if i < input_vec.len() - 1 {
      result[i] = 400.0 * input_vec[i].powi(3) - 400.0 * input_vec[i] * input_vec[i + 1]
        + 202.0 * input_vec[i]
        - 200.0 * input_vec[i - 1].powi(2)
        - 2.0;
    } else {
      result[i] = 200.0 * input_vec[i] - 200.0 * input_vec[i - 1].powi(2);
    }
  }
  return result;
}

fn diff_of_squares(min: f64, max: f64, power: i32) -> f64 {
  return max.powi(power) - min.powi(power);
}

fn gradient_sensitivity_s1(min: f64, max: f64, order: i32) -> f64 {
  let c0 = 160000.0 / 7.0 * diff_of_squares(min, max, 7)
    + 1600.0 / 5.0 * diff_of_squares(min, max, 5)
    - 1600.0 / 4.0 * diff_of_squares(min, max, 4)
    + 4.0 / 3.0 * diff_of_squares(min, max, 3)
    - 4.0 * diff_of_squares(min, max, 2)
    + 4.0 * diff_of_squares(min, max, 1);
  let c1 = -320000.0 / 5.0 * diff_of_squares(min, max, 5)
    - 1600.0 / 3.0 * diff_of_squares(min, max, 3)
    + 1600.0 * diff_of_squares(min, max, 2);
  let c2 = 1.0 / 3.0 * diff_of_squares(min, max, 3);

  let a = c0 * diff_of_squares(min, max, 1)
    + c1 * 0.5 * diff_of_squares(min, max, 2)
    + c2 / 3.0 * diff_of_squares(min, max, 3);

  return a * diff_of_squares(min, max, 1).powi(order - 2);
}

fn gradient_sensitivity_sj(min: f64, max: f64, order: i32) -> f64 {
  let a1 = 160_000.0 / 7.0 * diff_of_squares(min, max, 7) * diff_of_squares(min, max, 1).powi(2);
  let a2 = 1.0 / 5.0
    * (161_000.0 * diff_of_squares(min, max, 1) - 320_000.0 / 2.0 * diff_of_squares(min, max, 2))
    * diff_of_squares(min, max, 5)
    * diff_of_squares(min, max, 1);
  let a3 = 1.0 / 4.0
    * (-160_000.0 / 3.0 * diff_of_squares(min, max, 3) - 1600.0 * diff_of_squares(min, max, 1))
    * diff_of_squares(min, max, 4)
    * diff_of_squares(min, max, 1);
  let a4 = 1.0 / 3.0
    * (40804.0 * diff_of_squares(min, max, 1) - 161_600.0 / 2.0 * diff_of_squares(min, max, 2)
      + 160_000.0 / 3.0 * diff_of_squares(min, max, 3))
    * diff_of_squares(min, max, 3)
    * diff_of_squares(min, max, 1);
  let a5 = 1.0 / 2.0
    * (-80800.0 / 3.0 * diff_of_squares(min, max, 3) * diff_of_squares(min, max, 1)
      - 808.0 * diff_of_squares(min, max, 1).powi(2)
      + 160_000.0 / 6.0 * diff_of_squares(min, max, 2) * diff_of_squares(min, max, 3)
      + 800.0 * diff_of_squares(min, max, 2) * diff_of_squares(min, max, 1))
    * diff_of_squares(min, max, 2);
  let a6 = (40000.0 / 5.0 * diff_of_squares(min, max, 5)
    + 800.0 / 3.0 * diff_of_squares(min, max, 3)
    + 4.0 * diff_of_squares(min, max, 1))
    * diff_of_squares(min, max, 1).powi(2);
  return (a1 + a2 + a3 + a4 + a5 + a6) * diff_of_squares(min, max, 1).powi(order - 3);
}

fn gradient_sensitivity_sn(min: f64, max: f64, order: i32) -> f64 {
  return (40000.0 / 3.0 * diff_of_squares(min, max, 3) * diff_of_squares(min, max, 1)
    - 40000.0 / 3.0 * diff_of_squares(min, max, 2) * diff_of_squares(min, max, 3)
    + 40000.0 / 5.0 * diff_of_squares(min, max, 1) * diff_of_squares(min, max, 5))
    * diff_of_squares(min, max, order - 2);
}

fn gaussian_integration() {
  let min = -2.046;
  let max = 2.046;
  let scale = (max - min) * 0.5;
  let shift = (max + min) * 0.5;

  let n = 5;
  let m = 3;
  let mvec = [-(3.0f64 / 5.0).sqrt(), 0.0, (3.0f64 / 5.0).sqrt()];
  let wvec = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];
}

fn main() {
  let min = -2.046;
  let max = 2.046;
  let order = 5;
  let s1 = gradient_sensitivity_s1(min, max, order);
  let sj = gradient_sensitivity_sj(min, max, order);
  let sn = gradient_sensitivity_sn(min, max, order);
  println!("{},{},{}", s1, sj, sn);
  let a = s1.max(sj).max(sn);
  println!("{},{},{}", s1 / a, sj / a, sn / a);
}
