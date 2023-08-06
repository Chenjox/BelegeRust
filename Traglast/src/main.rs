mod visualisation;

mod polplan;
mod svd_helper;

use nalgebra::{DMatrix, Dyn, OMatrix, SVector, U2, U3};

pub const ZERO_THRESHHOLD: f64 = 1e-10;
type S2Vec = SVector<f64, 2>;
type S3Vec = SVector<f64, 3>;
type Point2DMatrix = OMatrix<f64, U2, Dyn>;
type Beam2DMatrix = OMatrix<usize, U2, Dyn>;
type DofMatrix = OMatrix<i32, U3, Dyn>;

type RigidMatrix = OMatrix<f64, Dyn, Dyn>;

#[derive(Debug)]
pub enum CONSTRAINS {
  XRestrain,
  YRestrain,
  PhiRestrain,
}

impl CONSTRAINS {
  pub fn number(con: Self) -> usize {
    match con {
      Self::XRestrain => 1,
      Self::YRestrain => 2,
      Self::PhiRestrain => 4,
    }
  }
}

#[derive(Debug)]
pub struct Point2D {
  num: usize,
  coordinates: S2Vec,
  constraints: Vec<CONSTRAINS>,
}

impl Point2D {
  pub fn new(number: usize, x: f64, y: f64, constraints: Vec<CONSTRAINS>) -> Self {
    return Point2D {
      num: number,
      coordinates: S2Vec::new(x, y),
      constraints,
    };
  }
  pub fn new_unconstrained(number: usize, x: f64, y: f64) -> Self {
    return Point2D {
      num: number,
      coordinates: S2Vec::new(x, y),
      constraints: vec![],
    };
  }
  pub fn is_x_constraint(&self) -> bool {
    for i in &self.constraints {
      match i {
        CONSTRAINS::XRestrain => return true,
        _ => {}
      }
    }
    false
  }
  pub fn is_y_constraint(&self) -> bool {
    for i in &self.constraints {
      match i {
        CONSTRAINS::YRestrain => return true,
        _ => {}
      }
    }
    false
  }
  pub fn get_coordinates(&self) -> S2Vec {
    return self.coordinates;
  }
  pub fn get_constraints(&self) -> &Vec<CONSTRAINS> {
    return &self.constraints;
  }
  pub fn get_homogenious_coordinates(&self) -> S3Vec {
    return S3Vec::new(self.coordinates.x, self.coordinates.y, 1.0);
  }
}

#[derive(Debug)]
pub struct Beam2D {
  from: usize,
  to: usize,
  nplast: f64,
  mplast: f64,
}

impl Beam2D {
  pub fn new(from: usize, to: usize, nplast: f64, mplast: f64) -> Self {
    if from == to {
      panic!(
        "Illegal Beam created! From and To Point are the same! from: {} == to: {}",
        from, to
      );
    }
    if nplast < 0.0 {
      panic!("Illegal Beam created! Nplast is negative: {}", nplast);
    }
    if mplast < 0.0 {
      panic!("Illegal Beam created! Mplast is negative: {}", mplast);
    }
    return Beam2D {
      from,
      to,
      nplast,
      mplast,
    };
  }

  pub fn get_from_num(&self) -> usize {
    return self.from;
  }

  pub fn get_to_num(&self) -> usize {
    return self.to;
  }

  pub fn get_nplast(&self) -> f64 {
    return self.nplast;
  }

  pub fn get_mplast(&self) -> f64 {
    return self.mplast;
  }
}

#[derive(Debug)]
pub struct Fachwerk2D {
  points: Vec<Point2D>,
  beams: Vec<Beam2D>,
}

impl Fachwerk2D {
  pub fn new(points: Vec<Point2D>, beams: Vec<Beam2D>) -> Self {
    return Fachwerk2D { points, beams };
  }

  pub fn matrix_form(&self) -> (Point2DMatrix, Beam2DMatrix, DofMatrix, i32) {
    let map: Vec<usize> = self.points.iter().map(|f| f.num).collect();
    let mut p_matrix = Point2DMatrix::zeros(self.points.len());
    let mut b_matrix = Beam2DMatrix::zeros(self.beams.len());
    let mut dof_matrix = DofMatrix::from_element(self.points.len(), -1);
    //for i in 0..map.len() {
    //    if res == map[i] {
    //        break;
    //    }
    //}
    let mut dof_num = 0;
    for p in 0..self.points.len() {
      p_matrix[(0, p)] = self.points[p].coordinates.x;
      p_matrix[(1, p)] = self.points[p].coordinates.y;

      if !self.points[p].is_x_constraint() {
        dof_matrix[(0, p)] = dof_num;
        dof_num += 1;
      }
      if !self.points[p].is_y_constraint() {
        dof_matrix[(1, p)] = dof_num;
        dof_num += 1;
      }
    }
    for b in 0..self.beams.len() {
      let from = self.beams[b].from;
      let to = self.beams[b].to;
      let mut new_from = 0;
      let mut new_to = 0;
      for i in 0..map.len() {
        if from == map[i] {
          new_from = i;
        }
        if to == map[i] {
          new_to = i;
        }
      }
      b_matrix[(0, b)] = new_from;
      b_matrix[(1, b)] = new_to;
    }
    //println!("{:?},{},{}",map, p_matrix, b_matrix);
    return (p_matrix, b_matrix, dof_matrix, dof_num);
  }
}

fn remove_beams(beams: &mut Vec<Beam2D>, indeces: &Vec<usize>) -> Vec<Beam2D> {
  let mut res = Vec::new();
  let beams_total = beams.len();
  let mut remove_total = 0;
  for i in 0..beams_total {
      if indeces.contains(&i) {
          let i = i - remove_total;
          let b = beams.remove(i);
          remove_total += 1;
          res.push(b);
      }
  }
  return res;
}

fn get_tragwerk() -> Fachwerk2D {
  let mut v = Vec::new();

  // Knoten (noch nicht in Koordinaten die ich brauch)
  for i in 0..7 {
    v.push(Point2D::new_unconstrained(i + 1, 3.0 * i as f64, 7.0));
  }
  for i in 0..7 {
    if (10..=12).contains(&(i + 8)) {
      v.push(Point2D::new_unconstrained(i + 8, 3.0 * i as f64, 4.5));
    } else {
      v.push(Point2D::new_unconstrained(i + 8, 3.0 * i as f64, 3.0));
    }
  }
  for i in 0..2 {
    v.push(Point2D::new(
      i + 15,
      3.0 * i as f64,
      0.0,
      vec![CONSTRAINS::XRestrain, CONSTRAINS::YRestrain],
    ));
    v.push(Point2D::new(
      18 - i,
      3.0 * (1 - i) as f64 + 15.0,
      0.0,
      vec![CONSTRAINS::XRestrain, CONSTRAINS::YRestrain],
    ));
  }
  // Kanten
  let mut kant = Vec::new();
  for i in 0..7 {
    kant.push(Beam2D::new(i + 1, i + 8, 1.0, 1.0));
  }
  for i in 0..6 {
    kant.push(Beam2D::new(i + 1, i + 2, 1.0, 1.0));
    kant.push(Beam2D::new(i + 8, i + 9, 1.0, 1.0));
  }
  for i in 0..3 {
    kant.push(Beam2D::new(i + 1, i + 9, 1.0, 1.0));
    kant.push(Beam2D::new(7 - i, 13 - i, 1.0, 1.0));
  }
  for i in 0..2 {
    kant.push(Beam2D::new(i + 8, i + 15, 1.0, 1.0));
    kant.push(Beam2D::new(14 - i, 18 - i, 1.0, 1.0));
  }
  kant.push(Beam2D::new(8, 16, 1.0, 1.0));
  kant.push(Beam2D::new(14, 17, 1.0, 1.0));

  return Fachwerk2D::new(v, kant);
}

fn get_testtragwerk() -> Fachwerk2D {
  let points = vec![
    Point2D::new(
      0,
      0.,
      0.,
      vec![CONSTRAINS::XRestrain, CONSTRAINS::YRestrain],
    ),
    Point2D::new_unconstrained(1, 1., 0.),
    Point2D::new_unconstrained(2, 2., 0.),
    Point2D::new(
      3,
      3.,
      0.,
      vec![CONSTRAINS::XRestrain, CONSTRAINS::YRestrain],
    ),
  ];
  let beams = vec![
    Beam2D::new(0, 1, 10., 1.0),
    Beam2D::new(1, 2, 10., 1.0),
    Beam2D::new(2, 3, 10., 1.0),
  ];

  return Fachwerk2D::new(points, beams);
}

fn get_testtragwerk2() -> Fachwerk2D {
  let points = vec![
    Point2D::new(
      1,
      0.,
      0.,
      vec![CONSTRAINS::XRestrain, CONSTRAINS::YRestrain],
    ),
    Point2D::new_unconstrained(2, 0., 1.),
    Point2D::new_unconstrained(3, 1., 1.),
    Point2D::new(4, 1., 0., vec![]),
  ];
  let beams = vec![
    Beam2D::new(1, 2, 10., 1.0),
    Beam2D::new(2, 3, 10., 1.0),
    Beam2D::new(3, 4, 10., 1.0),
    Beam2D::new(1, 4, 10., 1.0),
    Beam2D::new(1, 3, 10., 1.0),
    Beam2D::new(2, 4, 10., 1.0),
  ];

  return Fachwerk2D::new(points, beams);
}

fn main() {
  let trag = get_testtragwerk2();
  let (points, beams, dofmatrix, dof_num) = trag.matrix_form();

  //let erdscheiben_connections = binomial(erdscheibe.len(), 2);
  let mut mat = RigidMatrix::zeros(beams.ncols(), points.ncols() * 2);

  let num_beams = beams.ncols();
  let num_points = points.ncols();
  for i in 0..num_beams {
    let from_point = beams[(0, i)];
    let to_point = beams[(1, i)];

    let x_diff = points[(0, from_point)] - points[(0, to_point)];
    let y_diff = points[(1, from_point)] - points[(1, to_point)];
    let length = (x_diff * x_diff + y_diff * y_diff).sqrt();
    //println!("{}",length);

    mat[(i, 2 * from_point)] = x_diff / length;
    mat[(i, 2 * from_point + 1)] = y_diff / length;

    mat[(i, 2 * to_point)] = -x_diff / length;
    mat[(i, 2 * to_point + 1)] = -y_diff / length;
  }
  println!("{}", dofmatrix);

  println!("{:3.5}", mat);

  let mut bool_mat = RigidMatrix::zeros(points.ncols() * 2, dof_num as usize);
  let mut counter = 0;
  for i in 0..num_points {
    if dofmatrix[(0, i)] != -1 {
      bool_mat[(2 * i, counter)] = 1.;
      counter += 1;
    }
    if dofmatrix[(1, i)] != -1 {
      bool_mat[(2 * i + 1, counter)] = 1.;
      counter += 1;
    }
  }
  println!("{}", bool_mat);
  // Matrix aller Freiheitsgrade
  let mat = mat * &bool_mat;
  let rank = mat.rank(ZERO_THRESHHOLD);

  println!("{:?}", mat.shape());

  visualisation::visualise("test.png", 300, 300, &points, &beams);

  //let (_Y, _N, X) = svd_helper::get_svd_decomp(&mat);

  let X = svd_helper::get_nullspace(&mat);
  println!(
    "{:3.5}",
    X.view((0, rank), (dof_num as usize, dof_num as usize - rank))
  );

  let nullspace = X.view((0, rank), (dof_num as usize, dof_num as usize - rank));
  let other = nullspace.transpose() * bool_mat.transpose();

  println!("{:3.5}", other);

  let mut points_defo = points.clone();

  for i in 0..num_points {
    points_defo[(0, i)] = points[(0, i)] + 1. * &other.row(0)[(2 * i)];
    points_defo[(1, i)] = points[(1, i)] + 1. * &other.row(0)[(2 * i + 1)];
  }
  println!("{:3.5}", points_defo);

  visualisation::visualise("test1.png", 300, 300, &points_defo, &beams);

  // rigid body physics!

  //let p = add_displacements(&v, &displace);
  //visualise(
  //    &"2.png",
  //    720,
  //    720,
  //    &p,
  //    &kant,
  //    &rigid,
  //    &erd,
  //    &DCoordMat::from_element(0, 0, Pol::new_not_existing()),
  //    &vec![],
  //);
}
