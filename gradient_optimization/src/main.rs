use num_dual::*;

trait GradientOptimization2D {
  /// The target function in question.
  /// Must be differentialable at all points inside the search domain
  fn target_function<D: DualNum<f64>>(&self, x1: D, x2: D) -> D;

  /// The contraint function(s)
  /// Must return a boolean which indicates a violation of the constraints.
  fn violates_constraint_functions<D: DualNum<f64>>(&self, x1: D, x2: D) -> bool;

  /// the start point for the optimization.
  /// must not violate the contraint function
  fn start_point(&self) -> (f64,f64);

  /// A ordered list of search points, starting with greatest.
  /// all values must be strictly positive.
  fn steps_sizes(&self) -> &Vec<f64>;

  /// The accepted error between iterations
  fn eps(&self) -> f64;

  fn gradient_minimize(&self) -> (f64, f64) {
    let start_point = self.start_point();
    let mut current_point = StaticVec::new_vec([start_point.0,start_point.1]);
    let current_dual_point = current_point.map(DualVec64::<2>::from).derive();
    let mut current = self.target_function(current_dual_point[0], current_dual_point[1]);

    let mut stepsize_index = 0;
    let max_stepsize_index = self.steps_sizes().len();
    let max_iter = 1000;
    let mut iter = 0;
    loop { // Ein Iterationsschritt
        iter = iter+1;
        if iter > max_iter {
          break;
        }
        let new_point = current_point - current.eps * self.steps_sizes()[stepsize_index] /current.eps.norm();
        let new_dual_point = new_point.map(DualVec64::<2>::from).derive();
        let new = self.target_function(new_dual_point[0], new_dual_point[1]);
        println!("{},{}",new_point,stepsize_index);
        if (new.re - current.re).abs() < self.eps() {
          break;
        }
        // Nebenbedingungen
        if self.violates_constraint_functions(new_dual_point[0], new_dual_point[1]) {
          stepsize_index = if stepsize_index +1 == max_stepsize_index {stepsize_index} else  {stepsize_index + 1}; // Verkleinerung der Schrittweite, sofern möglich.
          continue;// Neuer Iterationsschritt
        }
        if new.re < current.re { // Zielfunktion wurde verbessert
          current_point = new_point;
          current = new; 
          stepsize_index = if stepsize_index == 0 {stepsize_index} else  {stepsize_index - 1}; // Vergrößerung der Schrittweite, sofern möglich.
          continue;// Neuer Iterationsschritt
        } else {
          stepsize_index = if stepsize_index +1 == max_stepsize_index {stepsize_index} else  {stepsize_index + 1}; // Verkleinerung der Schrittweite, sofern möglich.
          continue;// Neuer Iterationsschritt
        }
        
    }
    
    //println!("{}",current_point);
    return (current_point[0],current_point[1]);
  }
}

struct BelegFunction {
  pub start_point: (f64, f64),
  pub search_steps: Vec<f64>,
  pub eps: f64
}

impl GradientOptimization2D for BelegFunction {
  fn target_function<D: DualNum<f64>>(&self, x1: D, x2: D) -> D {
      (x1 * x2).mul(100.0) + x1 / x2 + (x1 * x2).powi(-1)
  }
  fn violates_constraint_functions<D: DualNum<f64>>(&self, x1: D, x2: D) -> bool {
      !(x1.re() >= 0.25 && x1.re() <= 1.0 && x2.re() >= 0.25 && x2.re() <= 1.0)
  }
  fn start_point(&self) -> (f64,f64) {
      self.start_point
  }
  fn steps_sizes(&self) -> &Vec<f64> {
      &self.search_steps
  }
  fn eps(&self) -> f64 {
      self.eps
  }
}

fn main() {

  let opt = BelegFunction {
    start_point: (0.5,0.5),
    search_steps: vec![0.01,0.005,0.0025,0.001],
    eps: 0.01
  };
  let res = opt.gradient_minimize();
  println!("{},{}",res.0,res.1);
}