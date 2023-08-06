use std::fmt::{self, Display, Formatter};

use nalgebra::{Dyn, OMatrix, SVector, U2};

use crate::ZERO_THRESHHOLD;

type IncidenceMatrix = OMatrix<usize, Dyn, Dyn>;
type PolplanMatrix = OMatrix<Pol, Dyn, Dyn>;

type Point2DMatrix = OMatrix<f64, U2, Dyn>;
type Beam2DMatrix = OMatrix<usize, U2, Dyn>;

type S2Vec = SVector<f64, 2>;
type S3Vec = SVector<f64, 3>;

#[derive(Debug, Clone)]
pub struct RigidBody {
  points: Vec<usize>,
}

impl RigidBody {
  pub fn new(points: Vec<usize>) -> Self {
    return RigidBody { points };
  }
  fn is_joined_with(&self, other: &RigidBody) -> bool {
    let mut count = 0;
    for i in 0..self.points.len() {
      for j in 0..other.points.len() {
        if self.points[i] == other.points[j] {
          //println!("{},{},{},{}",i,j,self.points[i],other.points[j]);
          count += 1;
        }
      }
    }
    return count >= 2;
  }

  fn join_with(&self, other: RigidBody) -> RigidBody {
    let mut n_p = self.points.clone();
    for i in other.points {
      if !n_p.contains(&i) {
        n_p.push(i);
      }
    }
    return RigidBody { points: n_p };
  }

  fn contains_point(&self, pointnum: usize) -> bool {
    return self.points.contains(&pointnum);
  }
  fn get_connecting_point(&self, other: &Self) -> Option<usize> {
    for i in 0..self.points.len() {
      for j in 0..other.points.len() {
        if self.points[i] == other.points[j] {
          return Some(self.points[i]);
        }
      }
    }
    return None;
  }
}

/// Der zweite Index wird zum ersten hinzugefügt
/// Wenn die Erdscheibe als letztes gelistet ist, und ein körper zur erdscheibe hinzugefügt werden soll dann
/// ```
/// erdscheiben_index, index_der_scheibe, liste
/// ```
fn join_rigid_bodies(index1: usize, index2: usize, bodies: &mut Vec<RigidBody>) {
  let bod = bodies.remove(index2);
  if index1 < index2 {
    bodies[index1] = bodies[index1].join_with(bod);
  } else {
    bodies[index1 - 1] = bodies[index1 - 1].join_with(bod);
  }
}

#[derive(Debug, PartialEq, Clone)]
struct Pol {
  coordinates: S3Vec,
}

impl Pol {
  pub fn new(x: f64, y: f64, z: f64) -> Self {
    return Pol {
      coordinates: S3Vec::new(x, y, z),
    };
  }

  pub fn infer(p1: &Pol, p2: &Pol, q1: &Pol, q2: &Pol) -> Self {
    //Zuerst die Ebenen
    let p_plane = p1.coordinates.cross(&p2.coordinates);
    let q_plane = q1.coordinates.cross(&q2.coordinates);
    // Schnittgerade der Ebenen:
    let new_pol = p_plane.cross(&q_plane);
    // Ist der Schnitt vllt der Nullschnitt?
    if new_pol.norm() < ZERO_THRESHHOLD {
      return Pol::new(0.0, 0.0, 0.0);
    } else {
      return Pol {
        coordinates: new_pol.normalize(),
      };
    }
  }

  pub fn exists(&self) -> bool {
    return !self.coordinates.x.is_nan()
      && !self.coordinates.y.is_nan()
      && !self.coordinates.z.is_nan();
  }

  pub fn is_same(&self, other: &Pol) -> bool {
    return self.coordinates.cross(&other.coordinates).norm() < ZERO_THRESHHOLD;
  }
}

impl Display for Pol {
  fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
    let coords = &self.coordinates;
    return if coords.z.abs() > ZERO_THRESHHOLD {
      write!(
        f,
        "({},{},{})",
        coords.x / coords.z,
        coords.y / coords.z,
        1.0
      )
    } else {
      write!(f, "({},{},{})", coords.x, coords.y, 0.0)
    };
  }
}

fn beams_to_rigid_bodies(beams: &Beam2DMatrix) -> Vec<RigidBody> {
  let mut res = Vec::new();
  for beam in 0..beams.shape().1 {
    let points = vec![beams[(0, beam)], beams[(1, beam)]];
    res.push(RigidBody::new(points));
  }
  return res;
}

fn polplan(rigid_bodies: &Vec<RigidBody>, erdscheibe: &Vec<usize>, points: &Point2DMatrix) {
  let mut rigid_bodies = rigid_bodies.clone();
  rigid_bodies.push(RigidBody::new(erdscheibe.clone())); // letzter Rigidbody ist erdscheibe
                                                         //
  let erdscheiben_index = rigid_bodies.len() - 1;
  let scheiben = rigid_bodies.len() - 1;
  // Erdscheibe ist letzter Index!
  let mut polplan =
    PolplanMatrix::from_element(scheiben, scheiben, Pol::new(f64::NAN, f64::NAN, f64::NAN));
  // Alle Scheiben adjazent zur Erdscheibe bekommen einen Hauptpol
  let erdscheibe = &rigid_bodies[erdscheiben_index];
  for i in 0..scheiben {
    // offenkundige Hauptpole
    let scheibe = &rigid_bodies[i];
    if erdscheibe.is_joined_with(scheibe) {
      // sollte das der Fall sein, sind beide Scheiben Erdscheiben
      todo!("Two Rigidbodies share two points, are therefore a single rigid body!");
    }
    if let Some(point_index) = erdscheibe.get_connecting_point(scheibe) {
      let coords = points.column(point_index);
      polplan[(i, i)] = Pol::new(coords.x, coords.y, 1.0);
      //println!("({},{}) at {}",i,i,coords);
    }
  }
  // alle offenkundigen Nebenpole
  for i in 0..scheiben {
    let i_scheibe = &rigid_bodies[i];
    for j in i + 1..scheiben {
      let j_scheibe = &rigid_bodies[j];
      if i_scheibe.is_joined_with(j_scheibe) {
        todo!("Two Rigidbodies share two points, are therefore a single rigid body!");
      }
      if let Some(point_index) = i_scheibe.get_connecting_point(j_scheibe) {
        let coords = points.column(point_index);
        let new_pol = Pol::new(coords.x, coords.y, 1.0);
        if polplan[(i, j)].exists() && !polplan[(i, j)].is_same(&new_pol) {
          todo!("Widerspruch im Nebenpol, beide müssen zusammengeführt werden!")
        } else {
          polplan[(i, j)] = new_pol.clone();
          polplan[(j, i)] = new_pol;
        }
      }
    }
  }
  println!("{}", polplan);
  // Das offenkundige war Widerspruchsfrei, also kommt jetzt die Inzidenzmatrix
  // Für n Hauptpole gibt es ( n 2 ) hauptpolregeln und ( n 3 ) nebenpolregeln..
  // https://stackoverflow.com/a/57444840
}
