mod visualisation;

mod polplan;
mod svd_helper;

use faer_core::ComplexField;
use itertools::Itertools;
use nalgebra::{DMatrix, Dyn, OMatrix, SVector, U1, U2, U3};

pub const ZERO_THRESHHOLD: f64 = 1e-10;
type S2Vec = SVector<f64, 2>;
type S3Vec = SVector<f64, 3>;
type Point2DMatrix = OMatrix<f64, U2, Dyn>;
type Beam2DMatrix = OMatrix<usize, U2, Dyn>;
type Load2DMatrix = OMatrix<f64, U2, Dyn>;
type Load2DFlattenMatrix = OMatrix<f64, U1, Dyn>;
type DofMatrix = OMatrix<i32, U3, Dyn>;

type RigidMatrix = OMatrix<f64, Dyn, Dyn>;

#[derive(Debug, Clone)]
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

include!("main_tragwerke.rs");

#[derive(Debug, Clone)]
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

#[derive(Debug, Clone)]
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

  pub fn matrix_form(&self) -> (Point2DMatrix, Beam2DMatrix, DofMatrix, i32, Indexmap) {
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
    let indexmap = Indexmap::new_from_orig_to_comp(map);
    //println!("{:?},{},{}",map, p_matrix, b_matrix);
    return (p_matrix, b_matrix, dof_matrix, dof_num, indexmap);
  }
}

fn remove_beams(beams: &Vec<Beam2D>, indeces: &Vec<usize>) -> (Vec<Beam2D>, Vec<Beam2D>) {
  let mut res = Vec::new();
  let mut beams = beams.clone();
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
  return (beams, res);
}

pub struct Indexmap {
  orig_to_computed: Vec<usize>,
}

impl Indexmap {
  fn new_from_orig_to_comp(map: Vec<usize>) -> Self {
    return Self {
      orig_to_computed: map,
    };
  }

  fn map_orig_to_comp(&self, num: usize) -> usize {
    let mut comp_index = 0;
    for i in 0..self.orig_to_computed.len() {
      if num == self.orig_to_computed[i] {
        comp_index = i;
      }
    }
    return comp_index;
  }
}

pub struct Load2D {
  node: usize,
  loads: S2Vec,
}

impl Load2D {
  fn new_from_components(node: usize, x_load: f64, y_load: f64) -> Self {
    return Self {
      node,
      loads: S2Vec::new(x_load, y_load),
    };
  }
}

/// Returns a Matrix of every Loads x and y component, similar to the point matrix.
/// also returns a Vector indicating the node of the load
fn get_loading_matrix_form(loads: Vec<Load2D>, map: &Indexmap) -> (Vec<usize>, Load2DMatrix) {
  let nloads = loads.len();
  let mut result = Load2DMatrix::zeros(nloads);
  let mut map_vec = vec![0; nloads];

  for (index, load) in loads.iter().enumerate() {
    let comp_index = map.map_orig_to_comp(load.node);
    result[(0, index)] += load.loads[0];
    result[(1, index)] += load.loads[1];

    map_vec[index] = comp_index;
  }

  return (map_vec, result);
}

/// Flattens a 2xN matrix into a 2Nx1 Matrix more suitable for multiplication.
fn get_loading_flatten_form(
  map_vec: &Vec<usize>,
  load_matrix: &Load2DMatrix,
  num_points: usize,
) -> Load2DFlattenMatrix {
  let mut mat = Load2DFlattenMatrix::zeros(num_points * 2);
  for (index, point_index) in map_vec.iter().enumerate() {
    mat[(2 * point_index)] += load_matrix[(0, index)];
    mat[(2 * point_index + 1)] += load_matrix[(1, index)];
  }
  return mat;
}

fn get_inner_loading_matrix(
  removed_beams: &Vec<Beam2D>,
  indexmap: &Indexmap,
  points: &Point2DMatrix,
  defo_points: &Point2DMatrix,
) -> Load2DFlattenMatrix {
  // zuerst alle Beams in jeweils 2 Lasten umrechnen:
  let mut res = Vec::<Load2D>::new();
  for beam in removed_beams {
    let from_index = indexmap.map_orig_to_comp(beam.from);
    let to_index = indexmap.map_orig_to_comp(beam.to);
    // Coordinaten der Punkte
    let from_coords = points.column(from_index);
    let to_coords = points.column(to_index);

    let diff = to_coords - from_coords;
    let length = (diff.component_mul(&diff)).sum().sqrt();
    let angle = diff.y.atan2(diff.x);

    let cos = angle.cos();
    let sin = angle.sin();

    let from_deformation = defo_points.column(from_index);
    let to_deformation = defo_points.column(to_index);
    let diff_deformation = to_deformation - from_deformation;
    let length_change = (diff_deformation.component_mul(&diff_deformation))
      .sum()
      .sqrt();
    let signum = if (length_change - length).abs() > ZERO_THRESHHOLD && length >= length_change {
      -1.
    } else if (length_change - length).abs() > ZERO_THRESHHOLD && length <= length_change {
      1.
    } else {
      0.
    };

    let load_from = Load2D::new_from_components(
      beam.from,
      signum * cos * beam.nplast,
      signum * sin * beam.nplast,
    );
    let load_to = Load2D::new_from_components(
      beam.to,
      -signum * cos * beam.nplast,
      -signum * sin * beam.nplast,
    );
    res.push(load_from);
    res.push(load_to);
  }

  let (inde, mat) = get_loading_matrix_form(res, indexmap);
  let flatten_mat = get_loading_flatten_form(&inde, &mat, points.ncols());

  return flatten_mat;
}

fn get_rigidity_matrix(points: &Point2DMatrix, beams: &Beam2DMatrix) -> RigidMatrix {
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
  return mat;
}

fn get_boolean_matrix(dofs: &DofMatrix, dof_num: usize, num_points: usize) -> RigidMatrix {
  let mut bool_mat = RigidMatrix::zeros(num_points * 2, dof_num);
  let mut counter = 0;
  for i in 0..num_points {
    if dofs[(0, i)] != -1 {
      bool_mat[(2 * i, counter)] = 1.;
      counter += 1;
    }
    if dofs[(1, i)] != -1 {
      bool_mat[(2 * i + 1, counter)] = 1.;
      counter += 1;
    }
  }
  return bool_mat;
}

struct EigenvectorCache {
  found_eigenvectors: Vec<OMatrix<f64, Dyn, Dyn>>,
}

impl EigenvectorCache {
  pub fn new() -> Self {
    return Self {
      found_eigenvectors: Vec::new(),
    };
  }

  pub fn is_contained(&self, other: &OMatrix<f64, Dyn, Dyn>) -> (usize, bool) {
    for (i, vec) in self.found_eigenvectors.iter().enumerate() {
      if (vec.normalize() - other.normalize()).norm() <= ZERO_THRESHHOLD*1e5 {
        return (i, true);
      }
    }
    return (usize::MAX, false);
  }

  pub fn add_eigenvector(&mut self, other: &OMatrix<f64, Dyn, Dyn>) -> usize {
    self.found_eigenvectors.push(other.normalize());
    self.found_eigenvectors.len() - 1
  }
}

fn main() {
  let trag = get_tragwerk();
  let beams = trag.beams;
  let points = trag.points;

  let beam_count = beams.len();

  let mut count = [0; 10];
  let mut eigencache = EigenvectorCache::new();
  let mut min_traglast = f64::MAX;
  for j in 1..=3 {
    for i in (0..beam_count).combinations(j) {
      let (beams, removed_beams) = remove_beams(&beams, &i);
      let points = points.clone();
      let trag = Fachwerk2D {
        beams,
        points: points,
      };

      let (points, beams, dofmatrix, dof_num, indexmap) = trag.matrix_form();

      let num_points = points.ncols();
      //let erdscheiben_connections = binomial(erdscheibe.len(), 2);

      let mat = get_rigidity_matrix(&points, &beams);

      let bool_mat = get_boolean_matrix(&dofmatrix, dof_num as usize, num_points);

      let mat = mat * &bool_mat;
      let rank = mat.rank(ZERO_THRESHHOLD);

      if let Some(YX) = svd_helper::get_nullspace(&mat) {
        let x = YX.1; // Mechanism
        let nullspace = x.view((0, rank), (dof_num as usize, dof_num as usize - rank));
        let deformations = nullspace.transpose() * bool_mat.transpose();

        count[nullspace.ncols()] += 1;

        if nullspace.ncols() == 1 {
          let loads = get_loading();
          let (map_vec, matrix_loads) = get_loading_matrix_form(loads, &indexmap);
          let outer_work_matrix = get_loading_flatten_form(&map_vec, &matrix_loads, num_points);

          let norm_deformations = 1.
            / (deformations
              .iter()
              .map(|f| f.abs())
              .max_by(|a, b| a.partial_cmp(&b).unwrap()))
            .unwrap()
            * &deformations;

          

          let outer_work = outer_work_matrix.dot(&norm_deformations);

          let norm_deformations = outer_work.signum() * norm_deformations;
          let outer_work = outer_work.signum() * outer_work;

          let m = eigencache.is_contained(&norm_deformations);
          let mut k = 0;
          if m.1 {
            k = m.0;
            continue;
          } else {
            k = eigencache.add_eigenvector(&norm_deformations);
          }
          let k = k;

          let point_deformations = {
            let mut mut_norm_deformations = Point2DMatrix::zeros(num_points);
            for i in 0..num_points {
              mut_norm_deformations[(0, i)] = points[(0, i)] + norm_deformations.row(0)[(2 * i)];
              mut_norm_deformations[(1, i)] =
                points[(1, i)] + norm_deformations.row(0)[(2 * i + 1)];
            }
            mut_norm_deformations
          };

          let inner_loads =
            get_inner_loading_matrix(&removed_beams, &indexmap, &points, &point_deformations);

          //println!("{:2.2}",inner_loads.transpose());

          //println!("{:2.2}",norm_deformations.transpose());

          let inner_work = inner_loads.dot(&norm_deformations);

          let traglast = if outer_work.abs() > ZERO_THRESHHOLD {
            -inner_work / outer_work
          } else {
            f64::INFINITY
          };
          //println!("{:2.2}", traglast);

          if outer_work > ZERO_THRESHHOLD && traglast.is_finite() {
            min_traglast = traglast.min(min_traglast);
            let path = format!(
              "Z-Resultfiles\\test{}-{}-{:2.3}.png",
              k,
              count[nullspace.ncols()],
              traglast
            );

            visualisation::visualise(&path, 400, 300, &points, &beams, &inner_loads);

            let points_defo = point_deformations.clone();

            let path = format!("Z-Resultfiles\\test{}v-{}.png", k, count[nullspace.ncols()],);
            visualisation::visualise(&path, 400, 300, &points_defo, &beams, &outer_work_matrix);
          }
        }
      }
    }
  }

  println!("{:?},{:2.4}", count, min_traglast);

  //

  //let (_Y, _N, X) = svd_helper::get_svd_decomp(&mat);

  //println!("{:3.5}", other);

  //let mut points_defo = points.clone();

  //for i in 0..num_points {
  //  points_defo[(0, i)] = points[(0, i)] - 0.3 * &other.row(0)[(2 * i)];
  //  points_defo[(1, i)] = points[(1, i)] - 0.3 * &other.row(0)[(2 * i + 1)];
  //}
  //println!("{:3.5}", points_defo);

  //visualisation::visualise("test1.png", 300, 300, &points_defo, &beams);

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
