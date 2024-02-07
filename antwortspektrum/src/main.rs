use nalgebra::DVector;


type FVector = DVector<f64>;

#[inline(always)]
fn lagrange_n1(xi: f64) -> f64 {
  0.5 + xi * 0.5
}

#[inline(always)]
fn lagrange_n2(xi: f64) -> f64 {
  0.5 - xi * 0.5
}

fn get_lagrange_vector(xi: f64, eta: f64) -> FVector {
  let mut result = FVector::zeros(4);
  result[0] = lagrange_n1(xi)*lagrange_n1(eta);
  result[1] = lagrange_n1(xi)*lagrange_n2(eta);
  result[2] = lagrange_n2(xi)*lagrange_n2(eta);
  result[3] = lagrange_n2(xi)*lagrange_n1(eta);
  return result;
}

fn main() {

  
  println!("Hello World");
}
