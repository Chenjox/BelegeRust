use std::ops::Mul;

use compute::distributions::{Continuous, Gumbel, Normal};
use faer::{Mat, Scale};
use realiability::{LogNormal, EULER_MASCHERONI};

mod realiability;



fn system_matrix_1(l: f64, h: f64, wplast: f64) -> Mat<f64> {

  let mut system_matrix = Mat::<f64>::zeros(3, 3);

  system_matrix[(0,0)] = 5./3. * l;
  system_matrix[(1,1)] = -h;
  system_matrix[(2,2)] = -(1. + 7./3.)*wplast;

  return system_matrix;
}

fn system_matrix_12(l: f64, h: f64, wplast: f64) -> Mat<f64> {

  let mut system_matrix = Mat::<f64>::zeros(3, 3);

  system_matrix[(0,0)] = -5./3. * l;
  system_matrix[(1,1)] = h;
  system_matrix[(2,2)] = -(1. + 7./3.)*wplast;

  return system_matrix;
}

fn system_matrix_2(_l: f64, h: f64, wplast: f64) -> Mat<f64> {

  let mut system_matrix = Mat::<f64>::zeros(3, 3);

  system_matrix[(0,0)] = 0.0;
  system_matrix[(1,1)] = h;
  system_matrix[(2,2)] = -(5./3.)*wplast;

  return system_matrix;
}

fn system_matrix_22(_l: f64, h: f64, wplast: f64) -> Mat<f64> {

  let mut system_matrix = Mat::<f64>::zeros(3, 3);

  system_matrix[(0,0)] = -0.0;
  system_matrix[(1,1)] = -h;
  system_matrix[(2,2)] = -(5./3.)*wplast;

  return system_matrix;
}

fn get_system_matrix(l: f64, h: f64, wplast: f64, index: usize) -> Mat<f64> {
  return match index {
      0 => system_matrix_1(l, h, wplast),
      1 => system_matrix_12(l, h, wplast),
      2 => system_matrix_2(l, h, wplast),
      3 => system_matrix_22(l, h, wplast),
      _ => panic!()
  };
}

fn dependency_1() -> (Mat<f64>,Mat<f64>) {
  let mut dependency_matrix = Mat::<f64>::zeros(3, 2);
  let mut affine_vector = Mat::<f64>::zeros(3,1);

  dependency_matrix[(0,0)] = 6.;
  dependency_matrix[(1,0)] = 2.;
  dependency_matrix[(2,1)] = 1.;

  affine_vector[(0,0)] = 40.;

  return (dependency_matrix,affine_vector);
}

fn dependency_2() -> (Mat<f64>,Mat<f64>) {
  let mut dependency_matrix = Mat::<f64>::zeros(3, 3);
  let mut affine_vector = Mat::<f64>::zeros(3,1);

  dependency_matrix[(0,0)] = 6.;
  dependency_matrix[(1,1)] = 2.;
  dependency_matrix[(2,2)] = 1.;

  affine_vector[(0,0)] = 60.;

  return (dependency_matrix,affine_vector);
}

fn dependency_3() -> (Mat<f64>,Mat<f64>) {
  let mut dependency_matrix = Mat::<f64>::zeros(3, 2);
  let mut affine_vector = Mat::<f64>::zeros(3,1);

  dependency_matrix[(0,0)] = 13.;
  dependency_matrix[(1,0)] = 6.;
  dependency_matrix[(2,1)] = 1.;

  affine_vector[(0,0)] = 40.;

  return (dependency_matrix,affine_vector);
}

fn get_dependency_matrix(index: usize) -> (Mat<f64>,Mat<f64>) {
  return match index {
      0 => dependency_1(),
      1 => dependency_2(),
      2 => dependency_3(),
      _ => panic!()
  };
}

fn z_trafo_1() -> (Mat<f64>,Mat<f64>) {
  let mut scaling_matrix = Mat::<f64>::zeros(2, 2);
  let mut affine_vector = Mat::<f64>::zeros(2,1);

  scaling_matrix[(0,0)] = 5.;
  scaling_matrix[(1,1)] = 2.64e4;

  affine_vector[(0,0)] = 30.;
  affine_vector[(1,0)] = 28.8e4;

  return (scaling_matrix,affine_vector);
}

fn z_trafo_2() -> (Mat<f64>,Mat<f64>) {
  let mut scaling_matrix = Mat::<f64>::zeros(3, 3);
  let mut affine_vector = Mat::<f64>::zeros(3,1);

  scaling_matrix[(0,0)] = 5.;
  scaling_matrix[(1,1)] = 5.;
  scaling_matrix[(2,2)] = 2.64e4;

  affine_vector[(0,0)] = 30.;
  affine_vector[(1,0)] = 30.;
  affine_vector[(2,0)] = 28.8e4;

  return (scaling_matrix,affine_vector);
}

fn get_z_trafo_matrix(index: usize) -> (Mat<f64>,Mat<f64>) {
  return match index {
      0 => z_trafo_1(),
      1 => z_trafo_2(),
      _ => panic!()
  };
}

fn task12() {
  let l = 3.5;
  let h = 4.0;
  
  let wplast = 1.628e-3;

  //every constellation
  for i in 0..2 {
    let (dependency, affine_dependency) = get_dependency_matrix(i);
    let (sigma, affine) = get_z_trafo_matrix(i);

    // every failure plane
    for j in 0..4 {
      let system = get_system_matrix(l, h, wplast, j);
      //println!("{:?}",dependency);

      let linear_part = &system*&dependency*&sigma;
      let affine_part = &system* &affine_dependency + &system*&dependency*&affine;
      
      let one_vector = Mat::<f64>::full(1, affine_part.nrows(), 1.0);
      

      let v = &one_vector*linear_part;
      let u = &one_vector*affine_part;

      let index = u.norm_l2()/v.norm_l2();

      println!("Task {}: b_{} = {}",i+1,j+1,index);
    }
  }
}


fn Gumbel_inverse_nataf_transformation(beta: f64, mu: f64, norm: &Normal, y: f64) -> f64 {
  return mu - (-norm.cdf(y).ln()).ln()*beta;
}

fn task3() {
  let l = 3.5;
  let h = 4.0;
  
  let wplast = 1.628e-3;

  let (dependency, affine_dependency) = get_dependency_matrix(2);

  let mean = 25.0;
  let std_dev = 5.0;

  let beta = std_dev*6.0_f64.sqrt()/std::f64::consts::PI;
  let mu = mean - beta*EULER_MASCHERONI;

  let load = Gumbel::new(mu, beta);

  let festigkeit = LogNormal::new(28.8e4, 2.64e4, 19.9e4);

  let normed_normal = Normal::default();

  // every failure plane
  for j in 0..4 {
    let system = get_system_matrix(l, h, wplast, j);
    
    let linear_part = &system*&dependency;
    let affine_part = &system* &affine_dependency;
    //println!("{:?}",linear_part);
    
    let one_vector = Mat::<f64>::full(1, affine_part.nrows(), 1.0);

    let mut guess = Mat::<f64>::zeros(2, 1);
    let mut delta : f64 = 0.0;
    loop {
        let last_guess = guess.clone();

        let mut base_space = Mat::<f64>::zeros(2, 1);
        base_space[(0,0)] = Gumbel_inverse_nataf_transformation(beta, mu, &normed_normal, guess[(0,0)]);
        base_space[(1,0)] = festigkeit.inverse_nataf_transformation(guess[(1,0)]);

        if base_space.norm_l2().is_infinite() {
          println!("Base Space was not finite!");
          break;
        }

        let failure_function = &one_vector*&linear_part*&base_space + &one_vector*&affine_part;

        //println!("{:?}",failure_function);

        let mut prob_gradient = Mat::<f64>::zeros(2, 2);
        prob_gradient[(0,0)] = normed_normal.pdf(guess[(0,0)])/load.pdf(base_space[(0,0)]);
        prob_gradient[(1,1)] = festigkeit.inverse_nataf_transformation_dx(guess[(1,0)]);

        //println!("{:?}",prob_gradient);

        prob_gradient = &one_vector*&linear_part*prob_gradient;
        
        //println!("{:?}",prob_gradient);

        let gradient_normalizer = 1.0/prob_gradient.norm_l2();
        // now its alpha
        prob_gradient *= Scale(-gradient_normalizer);

        //println!("{:?},{:?},{:?},{:?}",failure_function,prob_gradient,guess,base_space);
        // delta
        delta = (failure_function * Scale(gradient_normalizer) + (&prob_gradient*guess).transpose()).norm_l2();

        guess = prob_gradient.transpose() * Scale(delta);
        //println!("{:?}",guess);
        if (last_guess - guess.clone()).norm_l2() < 1e-10 {
          break;
        }
    }
    

    println!("Task {}: b_{} = {}",3,j+1,delta);
  }
}


fn main() {
  task12();
  task3()



  //println!("{}",load.pdf(60.0));
  //println!("{}",fy.cdf(50.0e4));
}
