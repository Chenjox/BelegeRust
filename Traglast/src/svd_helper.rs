use std::cmp::Ordering;

use faer_svd::compute_svd;
use nalgebra::DMatrix;

use dyn_stack::{DynStack, GlobalMemBuffer, ReborrowMut};
use faer_core::{Mat, Parallelism};

type GenMatrix = DMatrix<f64>;

pub fn get_svd_decomp(mat: &GenMatrix) -> (GenMatrix, GenMatrix, GenMatrix) {
  let n = mat.nrows();
  let m = mat.ncols();

  let mut c = Mat::<f64>::zeros(n, m);

  for i in 0..n {
    for j in 0..m {
      c.write(i, j, mat[(i, j)]);
    }
  }

  let mut s = Mat::<f64>::zeros(n, m);
  let mut u = Mat::<f64>::zeros(n, n);
  let mut v = Mat::<f64>::zeros(m, m);

  let mut mem = GlobalMemBuffer::new(
    faer_svd::compute_svd_req::<f64>(
      n.max(n),
      m.max(m),
      faer_svd::ComputeVectors::No,
      faer_svd::ComputeVectors::Full,
      Parallelism::None,
      faer_svd::SvdParams::default(),
    )
    .unwrap(),
  );
  let mut stack = DynStack::new(&mut mem);

  compute_svd(
    c.as_ref(),
    s.as_mut().diagonal(),
    None,
    Some(v.as_mut()),
    1e-10,
    1e-10,
    Parallelism::None,
    stack.rb_mut(),
    faer_svd::SvdParams::default(),
  );

  let mut s_res = GenMatrix::zeros(n, m);
  for i in 0..n {
    for j in 0..m {
      s_res[(i, j)] = s.read(i, j);
    }
  }
  let mut u_res = GenMatrix::zeros(n, n);
  for i in 0..n {
    for j in 0..n {
      u_res[(i, j)] = u.read(i, j);
    }
  }
  let mut v_res = GenMatrix::zeros(m, m);
  for i in 0..m {
    for j in 0..m {
      v_res[(i, j)] = v.read(i, j);
    }
  }
  return (u_res, s_res, v_res);
}

fn get_eigenmatrix(mat: &GenMatrix) -> Option<GenMatrix> {
  if !mat.is_square() {
    panic!("Matrix is not square! {:?}",mat.shape())
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

pub fn get_nullspace(mat: &GenMatrix) -> Option<(GenMatrix,GenMatrix)> {
  let btb = mat.transpose() * mat;
  let bbt = mat * mat.transpose();

  let X = get_eigenmatrix(&btb);
  let Y = get_eigenmatrix(&bbt);
  //println!("{:1.2}", sorted_eigenvectors);

  let result =Y.zip(X);

  return result;
}
