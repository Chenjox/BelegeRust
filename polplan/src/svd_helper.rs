use std::cmp::Ordering;

use nalgebra::DMatrix;

type GenMatrix = DMatrix<f64>;


fn get_eigenmatrix(mat: &GenMatrix) -> Option<GenMatrix> {
  if !mat.is_square() {
    panic!("Matrix is not square! {:?}", mat.shape())
  }
  let eigen = (mat.clone()).symmetric_eigen();

  //println!("{:1.2}", eigen.eigenvalues);
  //println!("{:1.2}", eigen.eigenvectors);

  let eigenvectors = &eigen.eigenvectors;
  let eigenvalues = &eigen.eigenvalues;

  let mut count = 0;
  let mut eigenval_vec: Vec<(f64, usize)> = eigenvalues
    .iter()
    .map(|f| {
      count += 1;
      (*f, count - 1)
    })
    .collect();

  let mut has_errored = false;
  eigenval_vec.sort_by(|a, b| {
    a.0.partial_cmp(&b.0).unwrap_or_else(|| {
      has_errored = true;
      Ordering::Equal
    })
  });
  eigenval_vec.reverse();
  if has_errored {
    println!("{:?}", eigenval_vec);
    return None;
  }

  let eigenval_vec = eigenval_vec;
  //println!("{:?}", eigenval_vec);

  let mut sorted_eigenvectors = GenMatrix::zeros(mat.nrows(), mat.nrows());

  for i in 0..mat.nrows() {
    sorted_eigenvectors.set_column(i, &eigenvectors.column(eigenval_vec[i].1))
  }
  //println!("{:1.2}", sorted_eigenvectors);

  return Some(sorted_eigenvectors);
}

pub fn get_nullspace(mat: &GenMatrix) -> Option<(GenMatrix, GenMatrix)> {
  let btb = mat.transpose() * mat;
  let bbt = mat * mat.transpose();

  let X = get_eigenmatrix(&btb);
  let Y = get_eigenmatrix(&bbt);
  //println!("{:1.2}", sorted_eigenvectors);

  let result = Y.zip(X);

  return result;
}
