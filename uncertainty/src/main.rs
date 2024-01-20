use std::{f64::consts::PI, fs::File};

use fuzzy::{FuzzyAnalysis, FuzzyTriangularNumber};
use rand::{Rng, thread_rng};
use stochastic::{NormalDistributedVariable, StochasticAnalysis, StochasticVariable, LogNormalDistributedVariable};

use crate::fuzzy::{EmpiricalFuzzyNumber, FuzzyNumber};
use std::io::Write;

mod fuzzy;
mod optimizer;
mod stochastic;

struct ExcentricBeam {
  levels: Vec<f64>,
}

impl FuzzyAnalysis for ExcentricBeam {
  fn get_alpha_levels(&self) -> Vec<f64> {
    return self.levels.clone();
  }

  fn get_fuzzy_numbers(&self) -> Vec<Box<dyn FuzzyNumber>> {
    return vec![
      Box::new(FuzzyTriangularNumber::new(0.15, 0.2, 0.3)), // Excentrizität
      Box::new(FuzzyTriangularNumber::new(200.0, 220.0, 230.0)), // Kraft F
      Box::new(FuzzyTriangularNumber::new(0.39, 0.4, 0.41)), // Stützenradius
    ];
  }

  fn deterministic_solution_function(&self, input_parameters: &Vec<f64>) -> f64 {
    let excetricity = input_parameters[0];
    let force = input_parameters[1];
    let radius = input_parameters[2];

    let bending_moment = -force * excetricity;

    let area = PI * radius * radius;
    let ftm = PI * radius * radius * radius * radius;

    let result = force / area - radius * bending_moment / ftm; // In kN/m^2
    return result;
  }
}

fn get_inters(low: f64, high: f64, num_samples: u32) -> Vec<f64> {
  let step = 1. / num_samples as f64;
  let mut cum_step = 0.0;
  let mut result = Vec::new();
  for i in 0..=num_samples {
    let interpolant = low * (1. - cum_step) + high * cum_step;
    cum_step += step;
    result.push(interpolant.min(high));
  }
  result
}

struct Grenzzustandsfunktion {
  param: f64
}

impl<R: Rng> StochasticAnalysis<R> for Grenzzustandsfunktion {
  fn get_distributions(&self) -> Vec<Box<dyn StochasticVariable<R>>> {
      return vec![Box::new(NormalDistributedVariable::new(230.0, 50.0)),Box::new(LogNormalDistributedVariable::new(27.70e4,2.01e4,19.95e4))];
  }
  fn output_function(&self, input_vec: &Vec<f64>) -> f64 {
      return input_vec[0] * self.param - input_vec[1];
  }
}

fn main() {
  let mut rng = thread_rng();

  let stoch = Grenzzustandsfunktion {
    param: 1052.85319
  };

  let erg = stoch.stochastics_analysis(0.66, 1e-3, 0.95, &mut rng);

  println!("{:?}",erg.get_samples().len());

  for i in 1..20 {
    let qua = 1.0/20.0 * i as f64;
    println!("{}, {}", qua, erg.quantile(qua));
  }

  println!("a) = {}",erg.get_probability_of(0.0));

  let stoch = Grenzzustandsfunktion {
    param: 654.309
  };

  let erg = stoch.stochastics_analysis(0.95, 1e-3, 0.95, &mut rng);

  println!("{:?}",erg.get_samples().len());

  println!("b) = {}",erg.get_probability_of(0.0));

}

fn main2() {
  let ex = ExcentricBeam {
    levels: vec![0.0, 0.5, 1.0], //get_inters(0.0, 1.0, 3),
  };

  let eFN = ex.fuzzy_analysis();

  /*
    let eFNSamples = eFN.get_discrete_alpha_levels();

    println!("{:?}",eFNSamples);

    let alph1 = eFN.get_lower_discrete_alpha_level(0.6);
    let alph2 = eFN.get_upper_discrete_alpha_level(0.6);
    println!("{:?},{:?}",alph1, alph2);
  */
  let mut p = 0.0;
  for n in 0..=10 {
    let a = eFN.alpha_level_interval(p);
    println!("{:3.2},{:3.2},{:?}", p, a[1] - a[0], a);
    p += 0.1;
  }

  //println!("{:?}",get_inters(0., 1., 50));

  let sampling = eFN.membership_function_samplings(get_inters(0., 1., 100));

  let mut file = File::create("result.csv").unwrap();
  for record in sampling {
    write!(&mut file, "{},{}\n", record[0], record[1]).unwrap();
  }
}
