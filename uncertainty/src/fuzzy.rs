
use crate::optimizer::{self, simplex_optimization};

/// A public trait generalising the notion of a T1 Fuzzy Set.
/// Will be used in the context of the Alpha-Level-Optimization.
/// This considers only convex fuzzy numbers, non convex fuzzy sets are not considered
pub trait FuzzyNumber {
  /// Returns the interval for a igven alpha level.
  /// [level] is considered to be in the range [0,1].
  /// Should return a array containing [lower, upper].
  /// alpha level 0 should return the support of the Fuzzy Number.
  fn alpha_level_interval(&self, level: f64) -> [f64; 2];

  /// Should return the lower bound of the given alpha level
  fn min(&self, level: f64) -> f64;
  /// Should return the upper bound of the given alpha level
  fn max(&self, level: f64) -> f64;

  /// Returns the support of the Fuzzy Number.
  /// Should be same as alpha level 0.
  fn support(&self) -> [f64; 2] {
    return self.alpha_level_interval(0.0);
  }

  /// Returns the Area of the Fuzzy Number
  fn area(&self) -> f64;

  /// Returns the Variance of the Fuzzy Number
  fn variance(&self) -> f64;

  /// Returns the SHANNON Entropy of the Fuzzy Number
  fn entropy(&self) -> f64;

  /// Returns the Centroid of the Fuzzy Number
  fn centroid(&self) -> f64;

  /// Returns the Alpha Level Cut of the Fuzzy Number
  fn level_set_cut(&self) -> f64;

  /// Returns the Membership-Function of the Fuzzy Number, sampled for certain alpha_levels
  /// Useful for plotting a Fuzzy Number
  fn membership_function_samplings(&self, alphas: Vec<f64>) -> Vec<[f64; 2]> {

    let mut result = Vec::new();
    for alpha in alphas {
      let alpha_level_interval = self.alpha_level_interval(alpha);
      result.push([alpha_level_interval[0],alpha]);
      result.push([alpha_level_interval[1],alpha]);
    }
    // Sortieren
    result.sort_by(|a, b| a[0].partial_cmp(&b[0]).unwrap());
    return result;
  }
}

pub struct FuzzyTriangularNumber {
  lower: f64,
  middle: f64,
  upper: f64,
}

impl FuzzyTriangularNumber {
    pub fn new(lower: f64, middle: f64, upper: f64) -> Self {
      if lower <= middle && middle <= upper {
        return FuzzyTriangularNumber { lower, middle, upper}
      } else {
          panic!("The Values are not sorted in ascending order! {},{},{}", lower, middle, upper);
      }

    }
}

impl FuzzyNumber for FuzzyTriangularNumber {
  fn alpha_level_interval(&self, level: f64) -> [f64; 2] {
    return [
      (1. - level) * (self.lower) + level * self.middle,
      level * self.middle + (1.0 - level) * self.upper,
    ];
  }

  fn min(&self, level: f64) -> f64 {
    //level * (self.middle - self.lower) + self.lower
    (1. - level) * (self.lower) + level * self.middle
  }

  fn max(&self, level: f64) -> f64 {
    level * self.middle + (1.0 - level) * self.upper
  }

  fn area(&self) -> f64 {
    0.5 * (self.upper - self.lower)
  }

  fn centroid(&self) -> f64 {
    1.0 / 6.0
      * ((2.0 * self.middle + self.lower) * (self.middle - self.lower)
        + (2.0 * self.middle + self.upper) * (self.upper - self.middle))
      * 1.0
      / self.area()
  }

  fn variance(&self) -> f64 {
    todo!("Implement Variance")
  }

  fn entropy(&self) -> f64 {
    todo!("Schannon's Entropy for the Triangular Fuzzy Number")
  }

  fn level_set_cut(&self) -> f64 {
    (self.lower + 2.0 * self.middle + self.upper) / 4.0
  }
}

pub struct FuzzyTrapezoidalNumber {
  lower: f64,
  middle_lower: f64,
  middle_upper: f64,
  upper: f64,
}

impl FuzzyNumber for FuzzyTrapezoidalNumber {
  fn alpha_level_interval(&self, level: f64) -> [f64; 2] {
    return [
      level * (self.middle_lower - self.lower) + self.lower,
      (1.0 - level) * (self.upper - self.middle_upper) + self.middle_upper,
    ];
  }

  fn min(&self, level: f64) -> f64 {
    level * (self.middle_lower - self.lower) + self.lower
  }

  fn max(&self, level: f64) -> f64 {
    (1.0 - level) * (self.upper - self.middle_upper) + self.middle_upper
  }

  fn area(&self) -> f64 {
    0.5 * (self.upper - self.lower + self.middle_upper - self.middle_lower)
  }

  fn centroid(&self) -> f64 {
    todo!("Implement Centroid")
  }

  fn variance(&self) -> f64 {
    todo!("Implement Variance")
  }

  fn entropy(&self) -> f64 {
    todo!("Schannon's Entropy for the Triangular Fuzzy Number")
  }

  fn level_set_cut(&self) -> f64 {
    todo!("Implement level set cut")
  }
}

type AlphaLevel = (f64,[f64; 2]);

pub struct EmpiricalFuzzyNumber {
  samples: Vec<AlphaLevel>,
}

impl EmpiricalFuzzyNumber {

  pub fn new(samples: Vec<AlphaLevel>) -> Self {
    let mut samples = samples;
    samples.sort_by(|a,b| a.0.partial_cmp(&b.0).unwrap());
    return EmpiricalFuzzyNumber {
      samples
    };
  }

  pub fn get_discrete_alpha_levels(&self) -> Vec<f64> {
    self.samples.iter().map(|f| f.0).collect()
  }

  /// Gives the lower discrete alpha level, that is closest to the given one.
  pub fn get_lower_discrete_alpha_level(&self, alpha: f64) -> AlphaLevel {
    let levels = &self.samples;

    let mut result_level = (f64::NAN,[f64::NAN,f64::NAN]);
    for level in levels.iter().rev() {
      if level.0 <= alpha {
        result_level = *level;
        break;
      }
    }
    return result_level;
  }

  /// Gives the upper discrete alpha level, that is closest to the given one.
  pub fn get_upper_discrete_alpha_level(&self, alpha: f64) -> AlphaLevel {
    let levels = &self.samples;

    let mut result_level= (f64::NAN,[f64::NAN,f64::NAN]);
    for level in levels {
      if level.0 >= alpha {
        result_level = *level;
        break;
      }
    }
    return result_level;
  }

}

impl FuzzyNumber for EmpiricalFuzzyNumber {
  fn alpha_level_interval(&self, level: f64) -> [f64; 2] {
    // finden des Alpha-Level
    let alph1 = self.get_lower_discrete_alpha_level(level);
    let alph2 = self.get_upper_discrete_alpha_level(level);

    if alph1.0 == alph2.0 {
      return alph1.1;
    }

    let percent = ( level - alph1.0 )/(alph2.0 - alph1.0);

    let downer = alph2.1[0] * percent + (1. - percent) * alph1.1[0];
    let upper = alph2.1[1] * percent + (1. - percent) * alph1.1[1];

    return [downer, upper];
  }

  fn min(&self, level: f64) -> f64 {
    self.alpha_level_interval(level)[0]
  }

  fn max(&self, level: f64) -> f64 {
    self.alpha_level_interval(level)[1]
  }

  fn area(&self) -> f64 {
    let samples = &self.samples;

    let mut result_area = 0.0;

    for i in 0..samples.len()-1 {
      let alpha1 = samples[i];
      let alpha2 = samples[i+1];

      let diff_alpha = alpha2.0    - alpha1.0;
      let diff_lower = alpha1.1[1] - alpha1.1[0];
      let diff_upper = alpha2.1[1] - alpha2.1[0];

      //println!("{},{},{}",diff_alpha,diff_lower,diff_upper);

      let area = (diff_lower + diff_upper)*diff_alpha * 0.5;
      result_area += area;
    }

    return result_area;
  }

  fn centroid(&self) -> f64 {
    todo!("Implement Centroid")
  }

  fn variance(&self) -> f64 {
    todo!("Implement Variance")
  }

  fn entropy(&self) -> f64 {
    todo!("Schannon's Entropy for the Triangular Fuzzy Number")
  }

  fn level_set_cut(&self) -> f64 {
    todo!("Implement level set cut")
  }
}


pub trait FuzzyAnalysis {
    
    fn deterministic_solution_function(&self,input_parameters: &Vec<f64>) -> f64;

    fn get_fuzzy_numbers(&self) -> Vec<Box<dyn FuzzyNumber>>;

    fn get_alpha_levels(&self) -> Vec<f64>;

    fn fuzzy_analysis(&self) -> EmpiricalFuzzyNumber {
      let fuzz = self.get_fuzzy_numbers();
      let n_fuzz = fuzz.len();
      let mut alphas = self.get_alpha_levels();
      alphas.sort_by(|a, b| a.partial_cmp(b).unwrap());
      alphas.reverse();
      let alphas = alphas;
      println!("{:?}",alphas);

      let mut alpha_level_samples = Vec::<AlphaLevel>::new();

      let mut has_run_once = false;
      let mut initial_min = vec![0.0; n_fuzz];
      let mut initial_max = vec![0.0; n_fuzz];
      // Jedes Alpha Level
      for alpha_level in alphas {
        // Bekommt die Schranken
        let mut mins = vec![0.; n_fuzz];
        let mut maxs = vec![0.; n_fuzz];
        for fz_num in 0..n_fuzz {
          mins[fz_num] =fuzz[fz_num].min(alpha_level);
          maxs[fz_num] =fuzz[fz_num].max(alpha_level);
        }
        // Startpunkt in der Mitte, auf dem Alpha Level 1
        if !has_run_once {
          let mut initial = vec![0.0; n_fuzz];
          for i in 0..n_fuzz {
              initial[i] = 0.5 *mins[i] + 0.5 * maxs[i];
          }
          initial_min = initial.clone();
          initial_max = initial;
          has_run_once = true;
        }

        let f = |vl: &Vec<f64>| (self.deterministic_solution_function(vl));
        let f2 = |vl: &Vec<f64>| (-self.deterministic_solution_function(vl));

        // Minimum
        let min_point = simplex_optimization(&mins, &maxs, &initial_min, f);
        // Maximum
        let max_point = simplex_optimization(&mins, &maxs, &initial_max, f2);

        let min = self.deterministic_solution_function(&min_point);
        let max = self.deterministic_solution_function(&max_point);

        initial_min = min_point;
        initial_max = max_point;

        let e_alpha_level = (alpha_level, [min, max]);
        alpha_level_samples.push(e_alpha_level);
      }

      return EmpiricalFuzzyNumber::new(alpha_level_samples);
    }
}