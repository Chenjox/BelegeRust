mod visualisation;

use nalgebra::{Dynamic, OMatrix, SMatrix, SVector, U2};

const ZERO_THRESHHOLD: f64 = 1e-10;

type S2Vec = SVector<f64, 2>;
type S3Vec = SVector<f64, 3>;

type Point2DMatrix = OMatrix<f64, U2, Dynamic>;
type Beam2DMatrix = OMatrix<usize, U2, Dynamic>;

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
}

impl Fachwerk2D {
    pub fn new(points: Vec<Point2D>, beams: Vec<Beam2D>) -> Self {
        return Fachwerk2D { points, beams };
    }

    pub fn matrix_form(&self) -> (Point2DMatrix, Beam2DMatrix) {
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
        //println!("{:?},{},{}",map, p_matrix, b_matrix);
        return (p_matrix,b_matrix);
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

    return Fachwerk2D::new(v, kant);
}

fn beams_to_rigid_bodies(beams: &Beam2DMatrix) -> Vec<RigidBody> {
    let mut res = Vec::new();
    for beam in 0..beams.shape().1 {
        let points = vec![beams[(0,beam)],beams[(1,beam)]];
        res.push(RigidBody::new(points));
    }
    return res;
}

fn main() {
    println!("Hello World");

    let trag = get_tragwerk();
    let (points,beams) = trag.matrix_form();
    let rigidBodies = beams_to_rigid_bodies(&beams);
    println!("{:?}",rigidBodies);
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
