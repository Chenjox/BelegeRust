mod visualisation;

mod polplan;
mod svd_helper;

use nalgebra::{Dyn, OMatrix, SVector, U2};

pub const ZERO_THRESHHOLD: f64 = 1e-10;
type S2Vec = SVector<f64, 2>;
type S3Vec = SVector<f64, 3>;
type Point2DMatrix = OMatrix<f64, U2, Dyn>;
type Beam2DMatrix = OMatrix<usize, U2, Dyn>;

type RigidMatrix = OMatrix<f64, Dyn, Dyn>;

#[derive(Debug)]
pub struct Point2D {
  num: usize,
  coordinates: S2Vec,
}

impl Point2D {
  pub fn new(number: usize, x: f64, y: f64) -> Self {
    return Point2D {
      num: number,
      coordinates: S2Vec::new(x, y),
    };
  }
  pub fn get_coordinates(&self) -> S2Vec {
    return self.coordinates;
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
}

impl Beam2D {
  pub fn new(from: usize, to: usize, nplast: f64) -> Self {
    if from == to {
      panic!(
        "Illegal Beam created! From and To Point are the same! from: {} == to: {}",
        from, to
      );
    }
    if nplast < 0.0 {
      panic!("Illegal Beam created! Nplast is negative: {}", nplast);
    }
    return Beam2D { from, to, nplast };
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
}

#[derive(Debug)]
pub struct Fachwerk2D {
  points: Vec<Point2D>,
  beams: Vec<Beam2D>,
  erdscheibe: Vec<usize>,
}

impl Fachwerk2D {
  pub fn new(points: Vec<Point2D>, beams: Vec<Beam2D>, erdscheibe: Vec<usize>) -> Self {
    return Fachwerk2D {
      points,
      beams,
      erdscheibe,
    };
  }

  pub fn matrix_form(&self) -> (Point2DMatrix, Beam2DMatrix, Vec<usize>) {
    let map: Vec<usize> = self.points.iter().map(|f| f.num).collect();
    let mut p_matrix = Point2DMatrix::zeros(self.points.len());
    let mut b_matrix = Beam2DMatrix::zeros(self.beams.len());
    //for i in 0..map.len() {
    //    if res == map[i] {
    //        break;
    //    }
    //}
    for p in 0..self.points.len() {
      p_matrix[(0, p)] = self.points[p].coordinates.x;
      p_matrix[(1, p)] = self.points[p].coordinates.y;
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
    let mut erd_vec = vec![0; self.erdscheibe.len()];
    for i in 0..self.erdscheibe.len() {
      let old_index = self.erdscheibe[i];
      let mut new_index = 0;
      for j in 0..map.len() {
        if map[j] == old_index {
          new_index = j;
        }
      }
      erd_vec[i] = new_index;
    }
    //println!("{:?},{},{}",map, p_matrix, b_matrix);
    return (p_matrix, b_matrix, erd_vec);
  }
}


fn get_tragwerk() -> Fachwerk2D {
  let mut v = Vec::new();

  // Knoten (noch nicht in Koordinaten die ich brauch)
  for i in 0..7 {
    v.push(Point2D::new(i + 1, 3.0 * i as f64, 7.0));
  }
  for i in 0..7 {
    if (10..=12).contains(&(i + 8)) {
      v.push(Point2D::new(i + 8, 3.0 * i as f64, 4.5));
    } else {
      v.push(Point2D::new(i + 8, 3.0 * i as f64, 3.0));
    }
  }
  for i in 0..2 {
    v.push(Point2D::new(i + 15, 3.0 * i as f64, 0.0));
    v.push(Point2D::new(18 - i, 3.0 * (1 - i) as f64 + 15.0, 0.0));
  }
  // Kanten
  let mut kant = Vec::new();
  for i in 0..7 {
    kant.push(Beam2D::new(i + 1, i + 8, 1.0));
  }
  for i in 0..6 {
    kant.push(Beam2D::new(i + 1, i + 2, 1.0));
    kant.push(Beam2D::new(i + 8, i + 9, 1.0));
  }
  for i in 0..3 {
    kant.push(Beam2D::new(i + 1, i + 9, 1.0));
    kant.push(Beam2D::new(7 - i, 13 - i, 1.0));
  }
  for i in 0..2 {
    kant.push(Beam2D::new(i + 8, i + 15, 1.0));
    kant.push(Beam2D::new(14 - i, 18 - i, 1.0));
  }
  kant.push(Beam2D::new(8, 16, 1.0));
  kant.push(Beam2D::new(14, 17, 1.0));

  let erd = vec![18, 17, 16, 15];

  return Fachwerk2D::new(v, kant, erd);
}


fn get_testtragwerk() -> Fachwerk2D {
  let points = vec![
    Point2D::new(0, 10. * 0., 10. * 0.),
    Point2D::new(1, 10. * 2., 10. * 0.),
    Point2D::new(2, 10. * 3., 10. * 0.)
  ];
  let beams = vec![
    Beam2D::new(0, 1, 10.),
    Beam2D::new(1, 2, 10.)
  ];
  let erd = vec![0,2];

  return Fachwerk2D::new(points, beams, erd);
}

fn main() {
  let trag = get_testtragwerk();
  let (points, beams, erdscheibe) = trag.matrix_form();

  //let erdscheiben_connections = binomial(erdscheibe.len(), 2);
  let mut mat = RigidMatrix::zeros(beams.ncols(), points.ncols() * 2);

  let num_beams = beams.ncols();
  let num_points = points.ncols();
  for i in 0..num_beams {
    let from_point = beams[(0, i)];
    let to_point = beams[(1, i)];

    let x_diff = points[(0, from_point)] - points[(0, to_point)];
    let y_diff = points[(1, from_point)] - points[(1, to_point)];
    let length = (x_diff*x_diff + y_diff * y_diff).sqrt();
    println!("{}",length);
    if !erdscheibe.contains(&from_point) {
      mat[(i, 2 * from_point)] = x_diff/length;
      mat[(i, 2 * from_point + 1)] = y_diff/length;
    }
    if !erdscheibe.contains(&to_point) { 
      mat[(i, 2 * to_point)] = -x_diff/length;
      mat[(i, 2 * to_point + 1)] = -y_diff/length;
    }
  }

  visualisation::visualise("test.png", 300, 300, &points, &beams);

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
