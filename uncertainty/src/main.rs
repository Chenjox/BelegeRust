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
}

struct FuzzyTriangularNumber {
  lower: f64,
  middle: f64,
  upper: f64,
}

impl FuzzyNumber for FuzzyTriangularNumber {
  fn alpha_level_interval(&self, level: f64) -> [f64; 2] {
    return [
      level * (self.middle - self.lower) + self.lower,
      (1.0 - level) * (self.upper - self.middle) + self.middle,
    ];
  }

  fn min(&self, level: f64) -> f64 {
    level * (self.middle - self.lower) + self.lower
  }

  fn max(&self, level: f64) -> f64 {
    (1.0 - level) * (self.upper - self.middle) + self.middle
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

struct FuzzyTrapezoidalNumber {
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

struct EmpiricalFuzzyNumber {
  samples: Vec<AlphaLevel>,
}

impl EmpiricalFuzzyNumber {

  fn new(samples: Vec<AlphaLevel>) -> Self {
    let mut samples = samples;
    samples.sort_by(|a,b| a.0.partial_cmp(&b.0).unwrap());
    return EmpiricalFuzzyNumber {
      samples
    };
  }

  fn get_discrete_alpha_levels(&self) -> Vec<f64> {
    self.samples.iter().map(|f| f.0).collect()
  }

  /// Gives the lower discrete alpha level, that is closest to the given one.
  fn get_lower_discrete_alpha_level(&self, alpha: f64) -> AlphaLevel {
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
  fn get_upper_discrete_alpha_level(&self, alpha: f64) -> AlphaLevel {
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

fn main() {

  let vec = vec![
    (0.0,[1.0,4.0]),
    (1.0,[2.0,2.5]),
    (0.5,[1.5,3.5]),
  ];

  let eFN = EmpiricalFuzzyNumber::new(vec);

  let eFNSamples = eFN.get_discrete_alpha_levels();

  println!("{:?}",eFNSamples);

  let alph1 = eFN.get_lower_discrete_alpha_level(0.6);
  let alph2 = eFN.get_upper_discrete_alpha_level(0.6);
  println!("{:?},{:?}",alph1, alph2);

  let mut p = 0.0;
  for n in 0..=10 {
    let a = eFN.alpha_level_interval(p);
    println!("{:3.2},{:3.2},{:?}",p,a[1] - a[0],a);
    p += 0.1;
  }

  println!("{}",eFN.area());
  

  println!("Hello World");
}
