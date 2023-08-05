use faer_svd::compute_svd;
use nalgebra::{Dyn, OMatrix};

use dyn_stack::{DynStack, GlobalMemBuffer, ReborrowMut};
use faer_core::{Mat, Parallelism};

type GenMatrix = OMatrix<f64, Dyn, Dyn>;

fn get_svd_decomp(mat: GenMatrix) -> (u, s, v) {
  let n = mat.nrows();
  let m = mat.ncols();

  let mut c = Mat::<f64>::zeros(n, m);

  for i in 0..n {
    for j in 0..m {
      c.write(i, j, mat[(i, j)]);
    }
  }

  let mut s = Mat::<f64>::zeros(n, m);
  let mut u = Mat::<f64>::zeros(m, m);
  let mut v = Mat::<f64>::zeros(n, n);

  let mut mem = GlobalMemBuffer::new(
    faer_svd::compute_svd_req::<f64>(
      4096,
      n,
      faer_svd::ComputeVectors::Thin,
      faer_svd::ComputeVectors::Thin,
      Parallelism::None,
      faer_svd::SvdParams::default(),
    )
    .unwrap(),
  );
  let mut stack = DynStack::new(&mut mem);

  compute_svd(
    c.as_ref(),
    s.as_mut(),
    Some(u.as_mut()),
    Some(v.as_mut()),
    1e-10,
    1e-16,
    Parallelism::None,
    stack.rb_mut(),
    faer_svd::SvdParams::default(),
  )
}
