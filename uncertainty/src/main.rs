
/// A public trait generalising the notion of a T1 Fuzzy Set.
/// Will be used in the context of the Alpha-Level-Optimization.
/// This considers only convex fuzzy numbers, non convex fuzzy sets are not considered
pub trait FuzzyNumber {
    /// Returns the interval for a igven alpha level.
    /// [level] is considered to be in the range [0,1].
    /// Should return a array containing [lower, upper].
    /// alpha level 0 should return the support of the Fuzzy Number.
    fn alpha_level_interval(&self,level: f64) -> [f64; 2];

    /// Should return the lower bound of the given alpha level
    fn min(&self,level: f64) -> f64;
    /// Should return the upper bound of the given alpha level
    fn max(&self,level: f64) -> f64;

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
    upper: f64
}

struct FuzzyTrapezoidalNumber {
    lower: f64,
    middle_lower: f64,
    middle_upper: f64,
    upper: f64
}

struct EmpiricalFuzzyNumber {
    samples: Vec<[f64; 2]>
}

fn main() {
    println!("Hello World");
}