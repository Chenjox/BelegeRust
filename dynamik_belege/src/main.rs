use nalgebra::{SMatrix, SVector};

// Ein paar Typen für Operationen
type Matrix3x3 = SMatrix<f64, 3, 3>;
type Vector3 = SVector<f64, 3>;

type Matrix4x4 = SMatrix<f64, 4, 4>;
type Vector4 = SVector<f64, 4>;

fn main() {
  let M = Matrix3x3::new(1.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 1.0);
  let K = Matrix3x3::new(4.0, -2.0, 0.0, -2.0, 6.0, -4.0, 0.0, -4.0, 8.0);
  let V = Matrix3x3::new(
    0.6318, 1.0000, 0.0849, 1.0000, -0.1748, -0.2039, 0.5582, -0.1916, 1.0000,
  );

  let Mg = V.transpose() * M * V;
  let Kg = V.transpose() * K * V;
  println!("{},{}", Mg, Kg);

  let mut eigenfreq = vec![0.0; 3];
  for i in 0..3 {
    eigenfreq[i] = (Kg[(i, i)] / Mg[(i, i)]).sqrt()
  }
  println!("{:?}", eigenfreq);

  // Aufgabe 2

  let n = Vector3::new(1.0, 0.0, 1.0);
  let omega = 9.0;
  let v_0 = Vector3::new(0.0, 0.0, 0.0);
  let vdot_0 = Vector3::new(1.0, 0.0, 0.0);

  let mut diagonalmat = Matrix3x3::zeros();
  for i in 0..3 {
    diagonalmat[(i, i)] = 1.0 / (Kg[(i, i)] - omega * omega * Mg[(i, i)]);
  }

  let v_parti = V * diagonalmat * V.transpose() * n;

  println!("{}", v_parti);

  let mut m_inv = Matrix3x3::zeros();
  let mut m_inv_omega = Matrix3x3::zeros();
  for i in 0..3 {
    m_inv[(i, i)] = 1.0 / (Mg[(i, i)]);
    m_inv_omega[(i, i)] = 1.0 / (Mg[(i, i)] * eigenfreq[i]);
  }

  let a = m_inv * V.transpose() * M * (v_0 - v_parti);
  let b = m_inv_omega * V.transpose() * M * (vdot_0);

  println!("a = {}, b = {}", a, b);
}
