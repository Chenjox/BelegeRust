use nalgebra::{SMatrix, SVector};

// Ein paar Typen für Operationen
type Matrix2x2 = SMatrix<f64, 2, 2>;
type Vector2 = SVector<f64, 2>;

// Berechnen des ersten Eigenwerts (also der kleinste Eigenwert)
fn get_first_eigenvalue(matA: Matrix2x2, matB: Matrix2x2, start: Vector2) -> f64 {
    let mut v_i = start;
    let MAX_ITER = 4; // nach maximal n Schritten Abbrechen
    let mut iter = 0;
    loop {
      iter += 1;
      let z = (v_i.transpose() * matA * v_i)[0];
      let h = matB * v_i;
      let n = (v_i.transpose() * h)[0];

      let rayleigh = z / n;

      let h_norm = 1.0 / n.sqrt() * h;
      println!("{},{},{},{},{}",v_i[0], v_i[1], h_norm[0], h_norm[1], rayleigh);
      v_i = matA.try_inverse().unwrap() * h_norm;

      if iter > MAX_ITER {
        break;
      }
    }

    return 0.0;
}

fn main() {
    let A = Matrix2x2::new( 1.0, -2.0, -2.0 , 5.0);
    let B = Matrix2x2::new( 1.0, 1.0, 1.0, 2.0);
    let start = Vector2::new( 1.0, 1.0 );
    get_first_eigenvalue(A, B, start);
}
