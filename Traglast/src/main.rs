mod visualisation;

use nalgebra::{Dynamic, OMatrix, OVector, SVector};
use std::collections::hash_map::DefaultHasher;
use std::fmt;
use std::hash::{Hash, Hasher};
use visualisation::visualise;

type DAdjUsize = OMatrix<usize, Dynamic, Dynamic>;
type DMatrixf64 = OMatrix<f64, Dynamic, Dynamic>;
type DVectorf64 = OVector<f64, Dynamic>;
type S2 = SVector<f64, 2>;
type S3 = SVector<f64, 3>;

fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}

#[derive(Debug, Clone)]
pub struct Point {
    pub num: usize,
    pub name: String,
    pub x: f64,
    pub y: f64,
}

impl Point {
    pub fn to_vector(&self) -> S2 {
        return S2::new(self.x, self.y);
    }
}

#[derive(Clone, PartialEq, Debug)]
pub struct Pol {
    // Angegeben in homogenen Koordiaten
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Clone, Hash, PartialEq)]
struct Polschlussregel([usize; 2], [usize; 2], [usize; 2]);

impl Polschlussregel {
    fn new(first: (usize, usize), second: (usize, usize), result: (usize, usize)) -> Self {
        let mut first_pair = if first.0 > first.1 {
            [first.1, first.0]
        } else {
            [first.0, first.1]
        };
        let mut second_pair = if second.0 > second.1 {
            [second.1, second.0]
        } else {
            [second.0, second.1]
        };
        let result = if result.0 > result.1 {
            [result.1, result.0]
        } else {
            [result.0, result.1]
        };
        //
        if first_pair[0] > second_pair[0] {
            let temp = first_pair;
            first_pair = second_pair;
            second_pair = temp;
        } else if first_pair[1] > second_pair[1] {
            let temp = first_pair;
            first_pair = second_pair;
            second_pair = temp;
        }

        return Polschlussregel(first_pair, second_pair, result);
    }
    fn first(&self) -> (usize, usize) {
        return (self.0[0], self.0[1]);
    }
    fn first_is_hauptpol(&self) -> bool {
        return self.0[0] == self.0[1];
    }
    fn second(&self) -> (usize, usize) {
        return (self.1[0], self.1[1]);
    }
    fn second_is_hauptpol(&self) -> bool {
        return self.1[0] == self.1[1];
    }
    fn result(&self) -> (usize, usize) {
        return (self.2[0], self.2[1]);
    }
    fn result_is_hauptpol(&self) -> bool {
        return self.2[0] == self.2[1];
    }
    fn rule_type(&self) -> u8 {
        return if self.first_is_hauptpol() && !self.second_is_hauptpol()
            || !self.first_is_hauptpol() && self.second_is_hauptpol()
        {
            2
        } else if self.first_is_hauptpol() && self.second_is_hauptpol() {
            1
        } else {
            0
        };
    }
}

impl fmt::Display for Pol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({0:.4}, {1:.4}, {2:.4})", self.x, self.y, self.z)
    }
}

impl Pol {
    fn new_not_existing() -> Self {
        return Pol {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
    }
    fn new(infty: bool, x: f64, y: f64) -> Self {
        if infty {
            Pol { x, y, z: 0.0 }
        } else {
            Pol { x, y, z: 1.0 }
        }
    }
    fn exists(&self) -> bool {
        return !((self.x == 0.0) && (self.y == 0.0) && (self.z == 0.0));
    }
    fn is_same(&self, other: &Pol) -> bool {
        let sv = S3::new(self.x, self.y, self.z);
        let ov = S3::new(other.x, other.y, other.z);
        let sv_max = sv.max();
        let ov_max = ov.max();
        let max = sv_max.max(ov_max);

        let cv = sv.cross(&ov);
        //println!("{},{}",cv.norm()/max < 1e-16, cv);
        return cv.norm() / max < 1e-10;
    }
}
impl Pol {
    fn infer(p1: &Pol, p2: &Pol, q1: &Pol, q2: &Pol) -> Option<Self> {
        let p1v = S3::new(p1.x, p1.y, p1.z);
        let p2v = S3::new(p2.x, p2.y, p2.z);
        let q1v = S3::new(q1.x, q1.y, q1.z);
        let q2v = S3::new(q2.x, q2.y, q2.z);
        // deriving the planes of possible points
        let eb1 = p1v.cross(&p2v);
        let eb2 = q1v.cross(&q2v);
        // schnitt der ebenen
        let result = eb1.cross(&eb2);
        if result[0] == 0.0 && result[1] == 0.0 && result[2] == 0.0 {
            // es gibt keinen Schnittpunkt (bzw. ist der unendliche Fernpunkt)
            return None;
        }
        let result = result.normalize();
        return Some(Pol {
            x: result[0],
            y: result[1],
            z: result[2],
        });
    }
    fn is_at_infinity(&self) -> bool {
        return self.z.abs() < 1e-12;
    }
    fn get_real_coordinates(&self) -> (f64, f64) {
        if self.z != 0.0 {
            return (self.x / self.z, self.y / self.z);
        } else if self.x != 0.0 {
            return (self.y / self.x, f64::INFINITY);
        } else {
            return (f64::INFINITY, f64::INFINITY);
        }
    }
}

type DCoordMat = OMatrix<Pol, Dynamic, Dynamic>;

#[derive(Debug)]
pub struct Beam {
    pub from: usize,
    pub to: usize,
    pub nplast: f64,
}

#[derive(Debug)]
struct Triangle(usize, usize, usize);

#[derive(Debug, Clone)]
pub struct RigidBody {
    pub points: Vec<usize>,
}

pub struct Mapping {
    pub map: Vec<usize>,
}

impl Mapping {
    // 0 -> 1
    pub fn map_to(&self, index: usize) -> usize {
        return self.map[index];
    }
    //
    pub fn map_from(&self, index: usize) -> usize {
        for i in 0..self.map.len() {
            if index == self.map[i] {
                return i;
            }
        }
        return usize::MAX;
    }
}

impl RigidBody {
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

impl Triangle {
    fn is_the_same(&self, other: &Triangle) -> bool {
        let s_elem = [self.0, self.1, self.2];
        let o_elem = [other.0, other.1, other.2];
        //Ist jedes Element der einen Menge in der anderen Enthalten?
        //
        let mut eq = true;
        for i in 0..3 {
            let mut is_member = false;
            for j in 0..3 {
                is_member = s_elem[i] == o_elem[j] || is_member;
            }
            eq = eq && is_member;
        }
        return eq;
    }

    fn contains_beam(&self, other: &Beam) -> bool {
        let s_elem = [self.0, self.1, self.2];
        let o_elem = [other.from, other.to];
        let mut is_contained = [false, false];
        // sind 2 elemente in der anderen enthalten?
        for i in 0..2 {
            let mut is_member = false;
            for j in 0..3 {
                is_member = s_elem[j] == o_elem[i] || is_member;
            }
            is_contained[i] = is_member;
        }
        return is_contained[0] && is_contained[1];
    }

    fn to_rigid_body(&self) -> RigidBody {
        return RigidBody {
            points: vec![self.0, self.1, self.2],
        };
    }
}

#[derive(Debug)]
pub struct Load {
    x: f64,
    y: f64,
    point: usize,
}

impl Load {
    fn new(x: f64, y: f64, point: usize) -> Self {
        Load { x, y, point }
    }
    pub fn to_vector(&self) -> S2 {
        return S2::new(self.x, self.y);
    }
}

fn get_loading() -> Vec<Load> {
    let mut res = Vec::new();

    res.push(Load::new(3.0, 0.0, 1));
    res.push(Load::new(0.0, -2.0, 4));
    res.push(Load::new(0.0, -1.0, 5));

    return res;
}

fn get_tragwerk() -> (Vec<Point>, Vec<Beam>) {
    let mut v = Vec::new();

    // Knoten (noch nicht in Koordinaten die ich brauch)
    for i in 0..7 {
        v.push(Point {
            num: i + 1,
            name: format!("A{}", i + 1),
            x: (3.0 * i as f64),
            y: 7.0,
        });
    }
    for i in 0..7 {
        if (10..=12).contains(&(i + 8)) {
            v.push(Point {
                num: i + 8,
                name: format!("A{}", i + 1),
                x: (3.0 * i as f64),
                y: 4.5,
            });
        } else {
            v.push(Point {
                num: i + 8,
                name: format!("A{}", i + 1),
                x: (3.0 * i as f64),
                y: 3.0,
            });
        }
    }
    for i in 0..2 {
        v.push(Point {
            num: i + 15,
            name: format!("B{}", i + 1),
            x: (3.0 * i as f64),
            y: 0.0,
        });
        v.push(Point {
            num: 18 - i,
            name: format!("B{}", 4 - i),
            x: (3.0 * (1 - i) as f64 + 15.0),
            y: 0.0,
        });
    }
    // Kanten
    let mut kant = Vec::new();
    for i in 0..7 {
        kant.push(Beam {
            from: i + 1,
            to: i + 8,
            nplast: 1.0,
        });
    }
    for i in 0..6 {
        kant.push(Beam {
            from: i + 1,
            to: i + 2,
            nplast: 1.0,
        });
        kant.push(Beam {
            from: i + 8,
            to: i + 9,
            nplast: 1.0,
        });
    }
    for i in 0..3 {
        kant.push(Beam {
            from: i + 1,
            to: i + 9,
            nplast: 1.0,
        });
        kant.push(Beam {
            from: 7 - i,
            to: 13 - i,
            nplast: 1.0,
        });
    }
    for i in 0..2 {
        kant.push(Beam {
            from: i + 8,
            to: i + 15,
            nplast: 1.0,
        });
        kant.push(Beam {
            from: 14 - i,
            to: 18 - i,
            nplast: 1.0,
        });
    }
    kant.push(Beam {
        from: 8,
        to: 16,
        nplast: 1.0,
    });
    kant.push(Beam {
        from: 14,
        to: 17,
        nplast: 1.0,
    });

    //kant.swap_remove(3);
    kant.swap_remove(4);
    kant.swap_remove(7);
    kant.swap_remove(27);

    return (v, kant);
}

fn get_rigid_bodies(v: &Vec<Point>, kant: &Vec<Beam>, erdscheibe: &RigidBody) -> Vec<RigidBody> {
    let vsort = sortPoints(v);
    let mapping = Mapping { map: vsort.clone() };
    let mut adj = get_adjacency_matrix(v, &vsort, kant);
    {
        //Einarbeiten der Erdscheibe
        for i in 0..erdscheibe.points.len() {
            for j in i + 1..erdscheibe.points.len() {
                let to = mapping.map_from(erdscheibe.points[i]);
                let from = mapping.map_from(erdscheibe.points[j]);
                adj[(to, from)] = 1;
                adj[(from, to)] = 1;
            }
        }
    }
    let adj = adj;

    let adj_square = &adj * &adj;

    //println!("{},{}", adj, adjSquare);

    let mut triangle_list: Vec<Triangle> = Vec::new();

    for i in 0..v.len() {
        for j in i + 1..v.len() {
            if adj[(i, j)] > 0 && adj_square[(i, j)] > 0 {
                // das finden der dritten Node des Dreiecks
                for k in i + 1..v.len() {
                    if adj[(i, k)] > 0 && adj[(k, j)] > 0 {
                        let tri = Triangle(vsort[i], vsort[j], vsort[k]);
                        // aber wir könnten es bereits gefunden haben!
                        let mut already_found = false;
                        for t in 0..triangle_list.len() {
                            already_found = already_found || triangle_list[t].is_the_same(&tri);
                        }
                        if !already_found {
                            //println!("Triangle found: {:?}", tri);
                            triangle_list.push(tri);
                        }
                        //else {
                        //    println!("Triangle already found: {:?}", tri);
                        //}
                    }
                }
            }
        }
    }
    // Jetzt alle in den Dreiecken enthaltenen Kanten filtern,
    let mut kant_not = Vec::new();
    for i in 0..kant.len() {
        let mut is_in_triangle = false;
        for t in &triangle_list {
            is_in_triangle = is_in_triangle || t.contains_beam(&kant[i]);
        }
        if !is_in_triangle {
            kant_not.push(i);
        }
    }
    //println!("{:?}", kant_not);
    // Abschließend alle Dreiecke zu Starrkörpern zusammenfügen.
    let mut rigid: Vec<RigidBody> = triangle_list
        .iter()
        .map(|tri| tri.to_rigid_body())
        .collect();
    // Alle Rigid Bodys die zusammengehören, zusammenbasteln

    // alle Rigidbodys zusammenführen
    //println!("{:?}\n", rigid);

    for _i in 0..rigid.len() {
        let mut anz_rem = 0;
        for cur_i in 0..rigid.len() - 1 {
            let cur_i = if cur_i > anz_rem {
                cur_i - anz_rem
            } else {
                continue;
            };
            for j in 1..rigid.len() - cur_i {
                if j > anz_rem && rigid[cur_i].is_joined_with(&rigid[cur_i + j - anz_rem]) {
                    polplan_join_bodies(cur_i, cur_i + j - anz_rem, &mut rigid);
                    anz_rem += 1;
                }
            }
        }
    }
    // Und einzelne Kanten zu Rigidbodys hinzufügen
    for i in kant_not {
        let from = kant[i].from;
        let to = kant[i].to;
        rigid.push(RigidBody {
            points: vec![from, to],
        })
    }
    //println!("{:?}\n", rigid);
    return rigid;
}

fn get_rigid_body_connectivity(
    points: &Vec<Point>,
    bodies: &Vec<RigidBody>,
) -> Vec<(usize, Vec<usize>)> {
    let mut res = Vec::new();
    for i in 0..points.len() {
        let mut bod_res = Vec::new();
        for j in 0..bodies.len() {
            if bodies[j].contains_point(points[i].num) {
                bod_res.push(j);
            }
        }
        if bod_res.len() > 1 {
            res.push((i, bod_res));
        }
    }
    return res;
}

// Eine optimierte indizierung
pub fn sortPoints(pvec: &Vec<Point>) -> Vec<usize> {
    let mut sortVec = vec![0; pvec.len()];

    // hier könnte ich McKee Cuthill machen
    // Einsortieren
    for i in 0..pvec.len() {
        sortVec[i] = pvec[i].num; // nummer von Punkt in pvec ist in sortVec enthalten
    }
    // Check, dass keine Duplikate enthalten sind.
    {
        let mut check_vec = Vec::with_capacity(pvec.len());
        for i in 0..pvec.len() {
            let test = sortVec[i];
            for j in 0..check_vec.len() {
                if test == check_vec[j] {
                    panic!("Double Indices found: {} is already given!", test);
                }
            }
            check_vec.push(test);
        }
    }
    return sortVec;
}

// erstellen der Adjazenzmatrix
fn get_adjacency_matrix(pvec: &Vec<Point>, sort_vec: &Vec<usize>, bvec: &Vec<Beam>) -> DAdjUsize {
    let mut adj_mat = DAdjUsize::zeros(pvec.len(), pvec.len());

    for i in 0..sort_vec.len() {
        for j in 0..bvec.len() {
            let b = &bvec[j];
            if sort_vec[i] == b.from {
                // index des von Punktes finden
                let f = b.to;
                for k in 0..sort_vec.len() {
                    if sort_vec[k] == f {
                        adj_mat[(i, k)] = 1;
                        adj_mat[(k, i)] = 1;
                    }
                }
            }
        }
    }
    //println!("{}",adj_mat);
    return adj_mat;
}
pub fn filter_erdscheibe_ans_ende(
    bodies: &Vec<RigidBody>,
    erdscheibe: &RigidBody,
) -> Vec<RigidBody> {
    let mut rigid = bodies.clone();
    let end = rigid.len() - 1;
    // finden der Erdscheibe,
    let mut erd = 0;
    'erdloop: for r_body_index in 0..rigid.len() {
        let r_body = &rigid[r_body_index];
        for i in &r_body.points {
            for j in &erdscheibe.points {
                if i == j {
                    erd = r_body_index;
                    break 'erdloop;
                }
            }
        }
        //rigid.push(erdscheibe.clone());
        // tauschen der Erdscheibe ans ende
    }
    rigid.swap(erd, end);
    //rigid.remove(end);
    return rigid;
}
/*
    Alle offenkundigen Hauptpole werden bestimmt.
    Die Erdscheibe ist bereits am Ende der Liste
*/
fn get_trivial_hauptpole(
    matrix: &mut DCoordMat,
    rigid: &Vec<RigidBody>,
    points: &Vec<Point>,
) -> (bool, usize) {
    let connect_points = get_rigid_body_connectivity(points, &rigid);
    // Mittels der Erdscheibe werden die offenkundigen Hauptpole bestimmt!
    let end = rigid.len() - 1;
    let erd = end;
    for i in 0..connect_points.len() {
        let con_bodies = &connect_points[i].1;
        // First pass, ist die Erdscheibe enthalten?
        let mut is_erde = false;
        for j in 0..con_bodies.len() {
            if con_bodies[j] == erd {
                is_erde = true;
            }
        }
        if is_erde {
            for j in 0..con_bodies.len() {
                if con_bodies[j] != erd {
                    let po = &points[connect_points[i].0];
                    let npol = Pol::new(false, po.x, po.y);
                    if matrix[(con_bodies[j], con_bodies[j])].exists()
                        && !matrix[(con_bodies[j], con_bodies[j])].is_same(&npol)
                    {
                        println!("Widerspruch im Hauptpol {}", con_bodies[j]);
                        return (false, con_bodies[j]);
                    } else {
                        matrix[(con_bodies[j], con_bodies[j])] = npol;
                    }
                }
            }
        }
    }
    return (true, usize::MAX);
}

fn get_trivial_nebenpole(
    matrix: &mut DCoordMat,
    rigid: &Vec<RigidBody>,
    points: &Vec<Point>,
) -> (bool, usize, usize) {
    let end = rigid.len() - 1;
    let connect_points = {
        let mut rigid = rigid.clone();
        rigid.remove(end);
        get_rigid_body_connectivity(points, &rigid)
    };
    //println!("{:?}", rigid);
    //println!("{:?}", connect_points);
    for i in 0..connect_points.len() {
        let con_bodies = &connect_points[i].1; // alle indezes der verbundenen Rigidbodies
        let po = &points[connect_points[i].0]; // der index des punktes in point
        for j in 0..con_bodies.len() {
            // mindestens 2 Elemente müssen enthalten sein
            for k in j + 1..con_bodies.len() {
                let npol = Pol::new(false, po.x, po.y);
                if matrix[(con_bodies[j], con_bodies[k])].exists()
                    && !matrix[(con_bodies[j], con_bodies[k])].is_same(&npol)
                {
                    println!(
                        "Widerspruch im Nebenpol {},{}",
                        con_bodies[j], con_bodies[k]
                    );
                    return (false, con_bodies[j], con_bodies[k]);
                } else {
                    matrix[(con_bodies[j], con_bodies[k])] = npol.clone();
                    matrix[(con_bodies[k], con_bodies[j])] = npol;
                }
            }
        }
    }
    return (true, usize::MAX, usize::MAX);
}
fn polplan_schlussregeln(
    mat: &DCoordMat,
    anzahl_pole: usize,
    regeln: &mut Vec<Polschlussregel>,
    already_done_rules: &Vec<u64>,
) {
    for i in 0..anzahl_pole {
        // Iterieren durch alle möglichen
        for j in i..anzahl_pole {
            // Der i,j Pol mit
            if mat[(i, j)].exists() {
                for k in 0..anzahl_pole {
                    for l in 0..anzahl_pole {
                        // mit dem k,l ten pol
                        if i == j && i == k && i != l && mat[(i, l)].exists() {
                            // (i,j) + (k,l) Remap (i,i) + (i,l) -> (l,l) HP+NP -> (HP)
                            regeln.push(Polschlussregel::new((i, i), (i, l), (l, l)));
                        }
                        if i == j && k == l && i != k && mat[(k, k)].exists() {
                            // (i,j) + (k,l) Remap (i,i) + (k,k) -> (i,k) HP+HP -> (NP)
                            regeln.push(Polschlussregel::new((i, i), (k, k), (i, k)));
                        }
                        if i != j && k != l && i == k && j != l && mat[(i, l)].exists() {
                            // (i,j) + (k,l) Remap (i,j) + (i,l) -> (j,l)
                            regeln.push(Polschlussregel::new((i, j), (i, l), (j, l)));
                        }
                    }
                }
            }
        }
    }
    // Entfernen von Duplikaten
    regeln.sort_by_key(|f| f.rule_type());
    //regeln.dedup();
    let mut cur = 0;
    loop {
        if cur >= regeln.len() {
            break;
        }
        let mut rm = 0;
        for i in 0..regeln.len() {
            let i = i - rm;
            if cur != i && calculate_hash(&regeln[i]) == calculate_hash(&regeln[cur]) {
                regeln.remove(i);
                rm += 1;
            }
        }
        cur += 1;
    }

    // Entfernen von befolgten regelen
    let anz_reg = regeln.len();
    let mut anz_remov = 0;
    for re in 0..anz_reg {
        let re = re - anz_remov;
        if already_done_rules.contains(&calculate_hash(&regeln[re])) {
            regeln.remove(re);
            anz_remov += 1;
        }
    }
}

/// Der zweite Index wird zum ersten hinzugefügt
/// Wenn die Erdscheibe als letztes gelistet ist, und ein körper zur erdscheibe hinzugefügt werden soll dann
/// ```
/// erdscheiben_index, index_der_scheibe, liste
/// ```
fn polplan_join_bodies(index1: usize, index2: usize, bodies: &mut Vec<RigidBody>) {
    let bod = bodies.remove(index2);
    if index1 < index2 {
        bodies[index1] = bodies[index1].join_with(bod);
    } else {
        bodies[index1 - 1] = bodies[index1 - 1].join_with(bod);
    }
}

// Aufstellen des Polplans
fn polplan(
    points: &Vec<Point>,
    bodies: &Vec<RigidBody>,
    erdscheibe: &RigidBody,
) -> (Vec<RigidBody>, DCoordMat) {
    let _mapping = Mapping {
        map: sortPoints(&points),
    };
    let mut rigid = filter_erdscheibe_ans_ende(bodies, erdscheibe);

    let mut mat;

    'se_big_loop: loop {
        let end = rigid.len() - 1;
        mat = DCoordMat::from_element(end, end, Pol::new_not_existing());

        //println!("{:?}", rigid);

        // Hier könnten schon widersprüche entstehen
        let (widerspruchsfrei, index) = get_trivial_hauptpole(&mut mat, &rigid, points);
        if !widerspruchsfrei {
            polplan_join_bodies(end, index, &mut rigid);
            continue 'se_big_loop;
        };
        let (widerspruchsfrei, ite, jte) = get_trivial_nebenpole(&mut mat, &rigid, points);
        if !widerspruchsfrei {
            polplan_join_bodies(ite, jte, &mut rigid);
            continue 'se_big_loop;
        };

        //println!("{:?},{}", rigid);
        let anzahl_pole = rigid.len() - 1;
        // Vector der bekannten Schlussregeln
        let mut regeln = Vec::new();
        let mut already_done_rules = Vec::new();

        polplan_schlussregeln(&mat, anzahl_pole, &mut regeln, &already_done_rules);
        //
        loop {
            // Der Polplan selbst
            //println!("{:?}",regeln);
            if regeln.len() <= 0 {
                // es gibt nur eine
                break;
            }
            let mut ig = usize::MAX;
            let mut jg = usize::MAX;
            for i in 0..regeln.len() {
                for j in i + 1..regeln.len() {
                    if regeln[i].2[0] == regeln[j].2[0] && regeln[i].2[1] == regeln[j].2[1] {
                        ig = i;
                        jg = j;
                    }
                }
            }
            if ig == usize::MAX || jg == usize::MAX {
                break;
            }

            //println!("{},{},{:?},{:?}", ig,jg,regeln[ig],regeln[jg]);
            let regel1 = &regeln[ig];
            let regel2 = &regeln[jg];

            // Diese Regeln werden befolgt, daher sind sie nicht doppelt zu befolgen.

            let (newpol_i, newpol_j) = regel1.result();

            let new_pol = Pol::infer(
                &mat[regel1.first()],
                &mat[regel1.second()],
                &mat[regel2.first()],
                &mat[regel2.second()],
            );
            let new_pol = match new_pol {
                Some(p) => p,
                None => {
                    already_done_rules.push(calculate_hash(&regel1));
                    already_done_rules.push(calculate_hash(&regel2));
                    polplan_schlussregeln(&mat, anzahl_pole, &mut regeln, &already_done_rules);
                    continue;
                }
            };
            //println!("{}", mat);
            //println!(
            //    "{:?},{:?} -> [{},{}] = {}",
            //    regel1, regel2, newpol_i, newpol_j, new_pol
            //);
            // Fallunterscheidung!
            if mat[(newpol_i, newpol_j)].exists() {
                // Der pol existiert schon!
                if newpol_i == newpol_j && !mat[(newpol_i, newpol_j)].is_same(&new_pol) {
                    // wir haben einen hauptpol
                    //println!(
                    //    "Widerspruch im Polplan, Hauptpol {} doppelt gefunden!\n{} != {}",
                    //    newpol_i,
                    //    new_pol,
                    //    mat[(newpol_i, newpol_j)]
                    //);
                    polplan_join_bodies(end, newpol_i, &mut rigid);
                    continue 'se_big_loop;
                    // betreffende Scheibe ist erdscheibe
                }
                if newpol_i != newpol_j && !mat[(newpol_i, newpol_j)].is_same(&new_pol) {
                    // wir haben einen nebenpol
                    //println!(
                    //    "Widerspruch im Nebenpol, Nebenpol {},{} doppelt gefunden!\n{} != {}",
                    //    newpol_i,
                    //    newpol_j,
                    //    new_pol,
                    //    mat[(newpol_i, newpol_j)]
                    //);
                    polplan_join_bodies(newpol_i, newpol_j, &mut rigid);
                    continue 'se_big_loop;
                    // Die Scheiben gehören zusammen.
                }
            }
            mat[(newpol_i, newpol_j)] = new_pol.clone();
            mat[(newpol_j, newpol_i)] = new_pol;
            already_done_rules.push(calculate_hash(&regel1));
            already_done_rules.push(calculate_hash(&regel2));
            // neue Schlussregeln
            polplan_schlussregeln(&mat, anzahl_pole, &mut regeln, &already_done_rules);
        }
        // Polplan ist widerspruchsfrei
        break 'se_big_loop;
    }
    //println!("{}", mat);
    return (rigid, mat);
}

fn kinematik(polplan: &DCoordMat) -> DVectorf64 {
    //let f = filter_erdscheibe_ans_ende(bodies, erdscheibe);
    //println!("{}", polplan);
    let anzahl_pole = polplan.shape().0;
    // suchen einer geeigneten Spalte (in der keine nicht-bestimmten Pole sind)
    // Annahme, alle Diagonalwerte sind besetzt!
    let mut candidate = 0;
    for i in 0..anzahl_pole {
        let mut is_accepted = true;
        for j in 0..anzahl_pole {
            is_accepted = is_accepted && polplan[(j, i)].exists();
        }
        if is_accepted && !polplan[(i, i)].is_at_infinity() {
            candidate = i;
            break;
        }
    }
    let candidate = candidate;

    let mut matrix = DMatrixf64::zeros(anzahl_pole, anzahl_pole);
    let hauptpol = &polplan[(candidate, candidate)];
    for i in 0..anzahl_pole {
        if i != candidate {
            let nebenpol = &polplan[(i, candidate)];
            let hp2 = &polplan[(i, i)];
            if nebenpol.is_at_infinity() {
                matrix[(i, candidate)] = 1.0;
                matrix[(i, i)] = -1.0;
            } else {
                let (x_hp1, y_hp1) = hauptpol.get_real_coordinates();
                let pos_hp1 = S2::new(x_hp1, y_hp1);
                let (x_np, y_np) = nebenpol.get_real_coordinates();
                let pos_np = S2::new(x_np, y_np);
                let direction = pos_np - pos_hp1;
                matrix[(i, candidate)] = (direction).norm();
                // TODO wenn hauptpol1 im unendlichen ist!
                if hp2.is_at_infinity() {
                    let normal = S2::new(direction[1], -direction[0]).normalize();
                    let trans_normal = S2::new(hp2.y, -hp2.x);
                    let trans_norm = trans_normal.norm();
                    let u = normal.dot(&trans_normal) / trans_norm;
                    matrix[(i, i)] = u;
                } else {
                    let (x_hp2, y_hp2) = hp2.get_real_coordinates();
                    let pos_hp2 = S2::new(x_hp2, y_hp2);
                    let direction = pos_np - pos_hp2;
                    matrix[(i, i)] = direction.norm();
                }
            }
        } else {
            matrix[(i, i)] = 1.0;
        }
    }
    //println!("{}", matrix);
    let mut b = DVectorf64::zeros(anzahl_pole);
    b[candidate] = 1.0;
    let lu = matrix.full_piv_lu();
    lu.solve_mut(&mut b);

    return b;
}

fn displacements(
    kinematik: &DVectorf64,
    polplan: &DCoordMat,
    rigid: &Vec<RigidBody>,
    points: &Vec<Point>,
) -> Vec<Point> {
    let mapping = Mapping {
        map: sortPoints(&points),
    };
    let mut result: Vec<Point> = points
        .clone()
        .iter()
        .map(|f| Point {
            num: f.num,
            name: f.name.clone(),
            x: 0.0,
            y: 0.0,
        })
        .collect();
    for i in 0..rigid.len() - 1 {
        let rigid_points = &rigid[i].points;
        let pol = &polplan[(i, i)];
        if pol.is_at_infinity() {
            let trans = S2::new(pol.y, -pol.x).normalize();
            for j in 0..rigid_points.len() {
                let point = &points[mapping.map_from(rigid_points[j])];
                let p_num = point.num;
                let p_name = point.name.clone();
                //let point = S2::new(point.x, point.y);
                let trans_point = trans * kinematik[i];
                result[mapping.map_from(rigid_points[j])] = Point {
                    num: p_num,
                    name: p_name,
                    x: trans_point.x,
                    y: trans_point.y,
                }
            }
        } else {
            let pos = pol.get_real_coordinates();
            let pos = S2::new(pos.0, pos.1);

            for j in 0..rigid_points.len() {
                let point = &points[mapping.map_from(rigid_points[j])];
                let p_num = point.num;
                let p_name = point.name.clone();
                let point = S2::new(point.x, point.y);
                let dir = pos - point;
                let dist = dir.norm();
                let dir_normal = S2::new(dir.y, -dir.x).normalize();
                let trans_point = dir_normal * dist * kinematik[i];
                result[mapping.map_from(rigid_points[j])] = Point {
                    num: p_num,
                    name: p_name,
                    x: trans_point.x,
                    y: trans_point.y,
                }
            }
        }
    }
    //println!("{:?}",result);
    return result;
}

fn main() {
    let (v, kant) = get_tragwerk();
    let loading = get_loading();
    let _mapping = Mapping {
        map: sortPoints(&v),
    };
    // Erdscheibe hinzufügen
    let erd = RigidBody {
        points: vec![15, 16, 17, 18],
    };

    let rigid = get_rigid_bodies(&v, &kant, &erd);

    let (rigid, pole) = polplan(&v, &rigid, &erd);
    println!("{}", pole);

    let kin = kinematik(&pole).normalize();
    println!("{}", kin);

    let displace = displacements(&kin, &pole, &rigid, &v);
    println!("{:?}", displace);

    visualise(&"1.png", 720, 720, &v, &kant, &rigid, &erd, &pole, &loading);
    let rigid = get_rigid_bodies(&v, &kant, &erd);
    let rigid = filter_erdscheibe_ans_ende(&rigid, &erd);
    visualise(
        &"2.png", 720, 720, &displace, &kant, &rigid, &erd, &pole, &loading,
    );
}
