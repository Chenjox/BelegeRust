use faer_svd::compute_svd;
use nalgebra::{Dyn, OMatrix};

use dyn_stack::{DynStack, GlobalMemBuffer, ReborrowMut};
use faer_core::{Mat, Parallelism};

type GenMatrix = OMatrix<f64, Dyn, Dyn>;

pub fn get_svd_decomp(mat: GenMatrix) -> (GenMatrix, GenMatrix, GenMatrix) {
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
      n,
      m,
      faer_svd::ComputeVectors::Full,
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
    Some(u.as_mut()),
    Some(v.as_mut()),
    1e-16,
    1e-16,
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
