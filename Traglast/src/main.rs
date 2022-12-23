mod visualisation;

use nalgebra::{Dynamic, OMatrix, OVector, SMatrix, SVector};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::{fmt, fs::FileType};
use visualisation::visualise;

type DAdjUsize = OMatrix<usize, Dynamic, Dynamic>;
type DMatrixf64 = OMatrix<f64, Dynamic, Dynamic>;
type DVectorf64 = OVector<f64, Dynamic>;
type S2x2 = SMatrix<f64, 2, 2>;
type S2 = SVector<f64, 2>;

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

#[derive(Clone, PartialEq, Debug)]
pub struct Pol {
    pub is_at_infinity: bool,
    pub x: f64,
    pub y: f64,
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
        if self.x.is_nan() {
            write!(f, "()")
        } else {
            write!(
                f,
                "({0}, {1:.4}, {2:.4})",
                if self.is_at_infinity { "Inf" } else { "000" },
                self.x,
                self.y
            )
        }
    }
}

impl Pol {
    fn new(infty: bool, x: f64, y: f64) -> Self {
        Pol {
            is_at_infinity: infty,
            x,
            y,
        }
    }
    fn exists(&self) -> bool {
        return !self.x.is_nan() && !self.y.is_nan();
    }
    fn is_same(&self, other: &Pol) -> bool {
        if (!self.is_at_infinity && other.is_at_infinity)
            || (!self.is_at_infinity && other.is_at_infinity)
        {
            return false;
        } // einer unendlich, der andere nicht
        if self.is_at_infinity && other.is_at_infinity {
            // beide unendlich
            return (self.x * other.y - other.y * self.x).abs() < 1e16;
        }
        return (self.x - other.x).abs() < 1e-16 && (self.y - other.y).abs() < 1e-16;
    }
}
impl Pol {
    fn infer(p1: &Pol, p2: &Pol, q1: &Pol, q2: &Pol) -> Self {
        let a_stuetz = if p1.is_at_infinity && !p2.is_at_infinity {
            // Annahme, p1 oder p2 sind endlich
            (p2.x, p2.y)
        } else {
            (p1.x, p1.y)
        };
        let a_richtung = if p1.is_at_infinity {
            (p1.x, p1.y)
        } else if p2.is_at_infinity {
            (p2.x, p2.y)
        } else {
            (p1.x - p2.x, p1.y - p2.y)
        };

        let b_stuetz = if q1.is_at_infinity && !q2.is_at_infinity {
            (q2.x, q2.y)
        } else {
            (q1.x, q1.y)
        };
        let b_richtung = if q1.is_at_infinity {
            (q1.x, q1.y)
        } else if q2.is_at_infinity {
            (q2.x, q2.y)
        } else {
            (q1.x - q2.x, q1.y - q2.y)
        };

        let mut mat = S2x2::zeros();
        let mut bvec = S2::zeros();

        mat[(0, 0)] = -a_richtung.1;
        mat[(0, 1)] = a_richtung.0;
        mat[(1, 0)] = -b_richtung.1;
        mat[(1, 1)] = b_richtung.0;

        bvec[0] = a_stuetz.0 * mat[(0, 0)] + a_stuetz.1 * mat[(0, 1)];
        bvec[1] = b_stuetz.0 * mat[(1, 0)] + b_stuetz.1 * mat[(1, 1)];

        let det = mat.determinant();

        if det.abs() != 0.0 {
            let new = mat.full_piv_lu().solve(&bvec).unwrap();
            return Pol::new(false, new[0], new[1]);
        } else {
            let new_x = b_richtung.0;
            let new_y = b_richtung.1;
            //println!("{},{},{},{}, [{},{}]",p1,p2,q1,q2,new_x,new_y);
            return Pol::new(true, new_x, new_y);
        }
    }

    fn distance(&self, other: &Pol) -> f64 {
        if !self.is_at_infinity && !other.is_at_infinity {
            return ((self.x - other.x).powi(2) + (self.y - other.y).powi(2)).sqrt();
        } else {
            return f64::INFINITY;
        }
    }

    fn is_between(&self, other: &Pol, another: &Pol) -> bool {
        return ((!other.is_at_infinity && other.is_at_infinity)
            || (other.is_at_infinity && !other.is_at_infinity))
            || ((other.x..=another.x).contains(&self.x)
                || (another.x..=other.x).contains(&self.x))
                && ((other.y..=another.y).contains(&self.y)
                    || (another.y..=other.y).contains(&self.y));
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

    fn is_joined_with(&self, other: &Triangle) -> bool {
        let s_elem = [self.0, self.1, self.2];
        let o_elem = [other.0, other.1, other.2];
        let mut is_contained = [false, false, false];
        // sind 2 elemente in der anderen enthalten?
        for i in 0..3 {
            let mut is_member = false;
            for j in 0..3 {
                is_member = s_elem[i] == o_elem[j] || is_member;
            }
            is_contained[i] = is_member;
        }
        let mut count = 0;
        for i in 0..3 {
            if is_contained[i] {
                count = count + 1;
            }
        }
        return count == 2;
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
        v.push(Point {
            num: i + 8,
            name: format!("A{}", i + 1),
            x: (3.0 * i as f64),
            y: 3.0,
        });
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

    let adjSquare = &adj * &adj;

    //println!("{},{}", adj, adjSquare);

    let mut triangle_list: Vec<Triangle> = Vec::new();

    for i in 0..v.len() {
        for j in i + 1..v.len() {
            if adj[(i, j)] > 0 && adjSquare[(i, j)] > 0 {
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

    for i in 0..rigid.len() {
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
fn get_adjacency_matrix(pvec: &Vec<Point>, sortVec: &Vec<usize>, bvec: &Vec<Beam>) -> DAdjUsize {
    let mut adj_mat = DAdjUsize::zeros(pvec.len(), pvec.len());

    for i in 0..sortVec.len() {
        for j in 0..bvec.len() {
            let b = &bvec[j];
            if sortVec[i] == b.from {
                // index des von Punktes finden
                let f = b.to;
                for k in 0..sortVec.len() {
                    if sortVec[k] == f {
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
                        if i == j && k == l && k > i && mat[(k, k)].exists() {
                            // (i,j) + (k,l) Remap (i,i) + (k,k) -> (i,k) HP+HP -> (NP)
                            regeln.push(Polschlussregel::new((i, i), (k, k), (i, k)));
                        }
                        if i != j && k != l && l > j && i == k && mat[(k, l)].exists() {
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
            let i = i -rm;
            if cur != i && calculate_hash(&regeln[i]) == calculate_hash(&regeln[cur]) {
                regeln.remove(i);
                rm += 1;
            }
        }
        cur +=1;
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
    let mapping = Mapping {
        map: sortPoints(&points),
    };
    let mut rigid = filter_erdscheibe_ans_ende(bodies, erdscheibe);

    let mut mat = DCoordMat::from_element(0, 0, Pol::new(false, f64::NAN, f64::NAN));

    'se_big_loop: loop {
        let end = rigid.len() - 1;
        mat = DCoordMat::from_element(end, end, Pol::new(false, f64::NAN, f64::NAN));

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
        // Es könnte auch sein, dass die Erdscheibe einen Stab vollständig

        //println!("{:?},{}", rigid);
        let anzahl_pole = rigid.len() - 1;
        // Vector der bekannten Schlussregeln
        let mut regeln = Vec::new();
        let mut already_done_rules = Vec::new();

        polplan_schlussregeln(&mat, anzahl_pole, &mut regeln, &already_done_rules);
        //
        loop {
            // Der Polplan selbst
            println!("{:?}",regeln);
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
            //println!("{}", mat);
            println!(
                "{:?},{:?} -> [{},{}] = {}",
                regel1, regel2, newpol_i, newpol_j, new_pol
            );
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

fn kinematik(polplan: &DCoordMat) -> DMatrixf64 {
    //let f = filter_erdscheibe_ans_ende(bodies, erdscheibe);

    let mut lin_op = DMatrixf64::zeros(polplan.shape().0, polplan.shape().0);
    for unabhaengig in 0..polplan.shape().0 {
        for i in 0..polplan.shape().0 {
            if i != unabhaengig {
                let hp1 = &polplan[(unabhaengig, unabhaengig)];
                let np = &polplan[(unabhaengig, i)];
                let hp2 = &polplan[(i, i)];
                if np.is_at_infinity {
                    //println!("Nebenpol im Unendlichen!");
                    lin_op[(i, unabhaengig)] = 1.0;
                } else if hp1.is_at_infinity {
                    // HP 1 ist im Unendlichen
                    // Es gibt keine Sinnvolle übersetzung von unendlich zu finit
                    // jedenfalls nicht im Anschauungsraum...
                    lin_op[(i, unabhaengig)] = 0.0;
                } else if hp2.is_at_infinity {
                    // HP 2 ist im Unendlichen
                    // Es gibt keine Sinnvolle übersetzung von unendlich zu finit
                    // jedenfalls nicht im Anschauungsraum...
                    lin_op[(i, unabhaengig)] = 0.0;
                } else if np.is_same(hp1) || np.is_same(hp2) {
                    // Jetzt können die nur noch koplanar sein!
                    lin_op[(i, unabhaengig)] = 0.0;
                } else {
                    let sign = if np.is_between(hp1, hp2) { -1.0 } else { 1.0 };
                    let phi_m = sign * np.distance(hp1) / np.distance(hp2);
                    lin_op[(i, unabhaengig)] = phi_m;
                }
            } else {
                if polplan[(unabhaengig, unabhaengig)].is_at_infinity {
                    // Im Unendlichen ergibt ein finiter Wert keinen Sinn!
                    lin_op[(unabhaengig, unabhaengig)] = 0.0;
                } else {
                    lin_op[(unabhaengig, unabhaengig)] = 1.0;
                }
            }
        }
    }
    return lin_op;
}

fn main() {
    let (v, kant) = get_tragwerk();
    let mapping = Mapping {
        map: sortPoints(&v),
    };
    // Erdscheibe hinzufügen
    let erd = RigidBody {
        points: vec![15, 16, 17, 18],
    };

    let rigid = get_rigid_bodies(&v, &kant, &erd);

    let (rigid, pole) = polplan(&v, &rigid, &erd);

    let kin = kinematik(&pole);
    println!("{}",kin);

    visualise(&"1.png", 720, 720, &v, &kant, &rigid, &erd, &pole);
    let rigid = get_rigid_bodies(&v, &kant, &erd);
    let rigid = filter_erdscheibe_ans_ende(&rigid, &erd);
    visualise(&"2.png", 720, 720, &v, &kant, &rigid, &erd, &pole);
}
