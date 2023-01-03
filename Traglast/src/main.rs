mod visualisation;

use std::fmt::{Display,Formatter};
use std::fmt;
use nalgebra::{Dynamic, OMatrix, SMatrix, SVector, U2};

const ZERO_THRESHHOLD: f64 = 1e-10;

type S2Vec = SVector<f64, 2>;
type S3Vec = SVector<f64, 3>;

type Point2DMatrix = OMatrix<f64, U2, Dynamic>;
type Beam2DMatrix = OMatrix<usize, U2, Dynamic>;
type IncidenceMatrix = OMatrix<usize, Dynamic, Dynamic>;
type PolplanMatrix = OMatrix<Pol, Dynamic, Dynamic>;

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
    erdscheibe: Vec<usize>
}

impl Fachwerk2D {
    pub fn new(points: Vec<Point2D>, beams: Vec<Beam2D>, erdscheibe: Vec<usize>) -> Self {
        return Fachwerk2D { points, beams, erdscheibe };
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
            p_matrix[(0,p)] = self.points[p].coordinates.x;
            p_matrix[(1,p)] = self.points[p].coordinates.y;
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
            b_matrix[(0,b)] = new_from;
            b_matrix[(1,b)] = new_to;
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
        return (p_matrix,b_matrix,erd_vec);
    }
}

#[derive(Debug, Clone)]
pub struct RigidBody {
    points: Vec<usize>,
}

impl RigidBody {
    pub fn new(points: Vec<usize>) -> Self {
        return RigidBody {
            points
        }
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
    coordinates: S3Vec
}

impl Pol {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        return Pol {
            coordinates: S3Vec::new(x,y,z)
        }
    }

    pub fn infer(p1: &Pol, p2: &Pol, q1: &Pol, q2: &Pol) -> Self {
        //Zuerst die Ebenen
        let p_plane = p1.coordinates.cross(&p2.coordinates);
        let q_plane = q1.coordinates.cross(&q2.coordinates);
        // Schnittgerade der Ebenen:
        let new_pol = p_plane.cross(&q_plane);
        // Ist der Schnitt vllt der Nullschnitt?
        if new_pol.norm() < ZERO_THRESHHOLD {
            return Pol::new(0.0,0.0,0.0);
        }else {
            return Pol {
                coordinates: new_pol.normalize()
            }
        }
    }

    pub fn exists(&self) -> bool {
        return !self.coordinates.x.is_nan() && !self.coordinates.y.is_nan() && !self.coordinates.z.is_nan();
    }

    pub fn is_same(&self, other: &Pol) -> bool {
        return self.coordinates.cross(&other.coordinates).norm() < ZERO_THRESHHOLD;
    }
}

impl Display for Pol {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let coords = &self.coordinates;
        return if coords.z.abs() > ZERO_THRESHHOLD {
            write!(f, "({},{},{})", coords.x/coords.z, coords.y/coords.z, 1.0)
        } else {
            write!(f, "({},{},{})", coords.x, coords.y, 0.0)
        }
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

    let erd = vec![18,17,16,15];

    return Fachwerk2D::new(v, kant, erd);
}

fn beams_to_rigid_bodies(beams: &Beam2DMatrix) -> Vec<RigidBody> {
    let mut res = Vec::new();
    for beam in 0..beams.shape().1 {
        let points = vec![beams[(0,beam)],beams[(1,beam)]];
        res.push(RigidBody::new(points));
    }
    return res;
}

fn polplan(rigid_bodies: &Vec<RigidBody>, erdscheibe: &Vec<usize>, points: &Point2DMatrix) {
    let mut rigid_bodies = rigid_bodies.clone();
    rigid_bodies.push(RigidBody::new(erdscheibe.clone())); // letzter Rigidbody ist erdscheibe
    //
    let erdscheiben_index = rigid_bodies.len()-1;
    let scheiben = rigid_bodies.len()-1;
    // Erdscheibe ist letzter Index!
    let mut polplan = PolplanMatrix::from_element(scheiben, scheiben,Pol::new(f64::NAN,f64::NAN,f64::NAN));
    // Alle Scheiben adjazent zur Erdscheibe bekommen einen Hauptpol
    let erdscheibe = &rigid_bodies[erdscheiben_index];
    for i in 0..scheiben { // offenkundige Hauptpole
        let scheibe = &rigid_bodies[i];
        if erdscheibe.is_joined_with(scheibe) { // sollte das der Fall sein, sind beide Scheiben Erdscheiben
            todo!("Two Rigidbodies share two points, are therefore a single rigid body!");
        }
        if let Some(point_index) = erdscheibe.get_connecting_point(scheibe) {
            let coords = points.column(point_index);
            polplan[(i,i)] = Pol::new(coords.x, coords.y, 1.0);
            //println!("({},{}) at {}",i,i,coords);
        }
    }
    // alle offenkundigen Nebenpole
    for i in 0..scheiben {
        let i_scheibe = &rigid_bodies[i];
        for j in i+1..scheiben {
            let j_scheibe = &rigid_bodies[j];
            if i_scheibe.is_joined_with(j_scheibe) {
                todo!("Two Rigidbodies share two points, are therefore a single rigid body!");
            }
            if let Some(point_index) = i_scheibe.get_connecting_point(j_scheibe) {
                let coords = points.column(point_index);
                let new_pol = Pol::new(coords.x, coords.y, 1.0);
                if polplan[(i,j)].exists() && !polplan[(i,j)].is_same(&new_pol) {
                    todo!("Widerspruch im Nebenpol, beide müssen zusammengeführt werden!")
                } else {
                    polplan[(i,j)] = new_pol.clone();
                    polplan[(j,i)] = new_pol;
                }
            }
        }
    }
    println!("{}",polplan);
    // Das offenkundige war Widerspruchsfrei, also kommt jetzt die Inzidenzmatrix
    // Für n Hauptpole gibt es ( n 2 ) hauptpolregeln und ( n 3 ) nebenpolregeln..
    // https://stackoverflow.com/a/57444840
}

fn main() {

    let trag = get_tragwerk();
    let (points,beams,erdscheibe) = trag.matrix_form();
    let rigid_bodies = beams_to_rigid_bodies(&beams);
    //println!("{:?}",rigid_bodies);
    polplan(&rigid_bodies, &erdscheibe, &points);
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
