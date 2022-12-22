mod visualisation;

use nalgebra::{Dynamic, OMatrix};
use std::fmt;
use visualisation::visualise;

type DAdjUsize = OMatrix<usize, Dynamic, Dynamic>;

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

#[derive(Debug, Clone)]
struct Polschlussregel([usize; 2], [usize; 2], [usize; 2]);

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
        return (self.x - other.x).abs() < 1e-16 && (self.y - other.y).abs() < 1e-16;
    }
}
impl Pol {
    fn infer(p1: &Pol, p2: &Pol, q1: &Pol, q2: &Pol) -> Self {
        let a_stuetz = if p1.is_at_infinity {
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

        let b_stuetz = if q1.is_at_infinity {
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

        let det = -a_richtung.0 * b_richtung.1 + a_richtung.1 * b_richtung.0;
        if det.abs() > 1e-16_f64 {
            let fx = b_stuetz.0 - a_stuetz.0;
            let fy = b_stuetz.1 - a_stuetz.1;

            let s = 1.0 / det * (-b_richtung.1 * fx - a_richtung.1 * fy);

            let new_x = b_stuetz.0 + s * b_richtung.0;
            let new_y = b_stuetz.1 + s * b_richtung.1;
            return Pol::new(false, new_x, new_y);
        } else {
            let new_x = b_richtung.0;
            let new_y = b_richtung.1;
            return Pol::new(true, new_x, new_y);
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

    kant.swap_remove(3);
    kant.swap_remove(4);
    //kant.swap_remove(7);
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
    let rgid_len = rigid.len();

    let mut cur_i = 0;
    loop {
        let mut has_changed = false;
        let mut anz_rem = 0;
        if cur_i < rigid.len() - 1 {
            for j in 1..rigid.len() - cur_i {
                if rigid[cur_i].is_joined_with(&rigid[cur_i + j - anz_rem]) {
                    let bod2 = rigid.remove(cur_i + j - anz_rem);
                    has_changed = has_changed || true;
                    anz_rem += 1;
                    rigid[cur_i] = rigid[cur_i].join_with(bod2);
                }
            }
        }
        if !has_changed {
            break;
        }
        cur_i += 1;
    }
    let mut cur_i = 0;
    loop {
        let mut has_changed = false;
        let mut anz_rem = 0;
        if cur_i < rigid.len() - 1 {
            for j in 1..rigid.len() - cur_i {
                if rigid[cur_i].is_joined_with(&rigid[cur_i + j - anz_rem]) {
                    let bod2 = rigid.remove(cur_i + j - anz_rem);
                    has_changed = has_changed || true;
                    anz_rem += 1;
                    rigid[cur_i] = rigid[cur_i].join_with(bod2);
                }
            }
        }
        if !has_changed {
            break;
        }
        cur_i += 1;
    }
    //println!("{:?}\n", rigid);

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

// Aufstellen des Polplans
fn polplan(points: &Vec<Point>, bodies: &Vec<RigidBody>, erdscheibe: &RigidBody) -> DCoordMat {
    let mapping = Mapping {
        map: sortPoints(&points),
    };
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

    // Im Unendlichen?, x, y, größe Abzüglich der erdscheibe
    let mut mat = DCoordMat::from_element(
        bodies.len() - 1,
        bodies.len() - 1,
        Pol::new(false, f64::NAN, f64::NAN),
    );

    {
        let connect_points = get_rigid_body_connectivity(points, &rigid);
        // Mittels der Erdscheibe werden die offenkundigen Hauptpole bestimmt!
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
                        if !mat[(con_bodies[j], con_bodies[j])].exists() {
                            mat[(con_bodies[j], con_bodies[j])] = Pol::new(false, po.x, po.y);
                        } else {
                            println!("Widerspruch im Polplan gefunden");
                        }
                    }
                }
            }
        }
        rigid.remove(end);
        let rigid = rigid;
        // Jetzt die Offenkundigen Nebenpole
        let connect_points = get_rigid_body_connectivity(points, &rigid);
        //println!("{:?}", rigid);
        //println!("{:?}", connect_points);
        for i in 0..connect_points.len() {
            let con_bodies = &connect_points[i].1; // alle indezes der verbundenen Rigidbodies
            let po = &points[connect_points[i].0]; // der index des punktes in point
            for j in 0..con_bodies.len() {
                // mindestens 2 Elemente müssen enthalten sein
                for k in j + 1..con_bodies.len() {
                    mat[(con_bodies[j], con_bodies[k])] = Pol::new(false, po.x, po.y);
                    mat[(con_bodies[k], con_bodies[j])] = Pol::new(false, po.x, po.y);
                }
            }
        }
    }
    //println!("{}", mat);
    let anzahl_pole = bodies.len() - 1;
    // Vector der bekannten Schlussregeln
    let mut regeln = Vec::new();

    for i in 0..anzahl_pole {
        // Iterieren durch alle möglichen
        for j in i..anzahl_pole {
            // Der i,j Pol mit
            if mat[(i, j)].exists() {
                for k in 0..anzahl_pole {
                    for l in k..anzahl_pole {
                        // mit dem k,l ten pol
                        if i == j && i == k && i != l && mat[(i, l)].exists() {
                            // (i,j) + (k,l) Remap (i,i) + (i,l) -> (l,l) HP+NP -> (HP)
                            regeln.push(Polschlussregel([i, i], [i, l], [l, l]));
                        }
                        if i == j && k == l && k > i && mat[(k, k)].exists() {
                            // (i,j) + (k,l) Remap (i,i) + (k,k) -> (i,k) HP+HP -> (NP)
                            regeln.push(Polschlussregel([i, i], [k, k], [i, k]));
                        }
                        if i != j && k != l && l > j && i == k && mat[(k, l)].exists() {
                            // (i,j) + (k,l) Remap (i,j) + (i,l) -> (j,l)
                            regeln.push(Polschlussregel([i, j], [i, l], [j, l]));
                        }
                    }
                }
            }
        }
    }

    loop {
        if regeln.len() == 0 {
            break;
        }
        let mut ig = 0;
        let mut jg = 0;
        for i in 0..regeln.len() {
            for j in i + 1..regeln.len() {
                if regeln[i].2[0] == regeln[j].2[0] && regeln[i].2[1] == regeln[j].2[1] {
                    ig = i;
                    jg = j;
                }
            }
        }

        //println!("{},{},{:?},{:?}", ig,jg,regeln[ig],regeln[jg]);
        let regel1 = regeln.remove(ig);
        let regel2 = regeln.remove(jg - 1); // da immer i < j gilt!
                                            //println!("{:?},{:?}", regel1, regel2);
                                            // Inferiere den Pol
                                            // Regeln sind ja schön, aber hier brauch ich ne Fallunterscheidung
                                            // Wenn (i,j) + (j,k) schon aufeinaner liegen, dann ist der dritte Pol auch dort!
        let newpol_i = regel1.2[0];
        let newpol_j = regel1.2[1];
        if mat[(regel1.0[0], regel1.0[1])].is_same(&mat[(regel1.1[0], regel1.1[1])]) {
            mat[(newpol_i, newpol_j)] = mat[(regel1.0[0], regel1.0[1])].clone();
            mat[(newpol_j, newpol_i)] = mat[(regel1.0[0], regel1.0[1])].clone();
        } else if mat[(regel2.0[0], regel2.0[1])].is_same(&mat[(regel2.1[0], regel2.1[1])]) {
            mat[(newpol_i, newpol_j)] = mat[(regel2.0[0], regel2.0[1])].clone();
            mat[(newpol_j, newpol_i)] = mat[(regel2.0[0], regel2.0[1])].clone();
        } else {
            let new_pol = Pol::infer(
                &mat[(regel1.0[0], regel1.0[1])],
                &mat[(regel1.1[0], regel1.1[1])],
                &mat[(regel2.0[0], regel2.0[1])],
                &mat[(regel2.1[0], regel2.1[1])],
            );

            //println!("[{},{}] = {}", newpol_i, newpol_j, new_pol);
            mat[(newpol_i, newpol_j)] = new_pol.clone();
            mat[(newpol_j, newpol_i)] = new_pol;
        }
        // neue Schlussregeln
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
                                regeln.push(Polschlussregel([i, i], [i, l], [l, l]));
                            }
                            if i == j && k == l && mat[(k, k)].exists() {
                                // (i,j) + (k,l) Remap (i,i) + (k,k) -> (i,k) HP+HP -> (NP)
                                regeln.push(Polschlussregel([i, i], [k, k], [i, k]));
                            }
                            if i != j && k != l && i == k && mat[(k, l)].exists() {
                                // (i,j) + (k,l) Remap (i,j) + (i,l) -> (j,l) NP+NP -> (NP)
                                regeln.push(Polschlussregel([i, j], [i, l], [j, l]));
                            }
                        }
                    }
                }
            }
        }
        let anz_reg = regeln.len();
        let mut anz_remov = 0;
        for re in 0..anz_reg {
            let re = re - anz_remov;
            if regeln[re].2[0] == newpol_i && regeln[re].2[1] == newpol_j
                || regeln[re].2[0] == newpol_j && regeln[re].2[1] == newpol_i
                || mat[(regeln[re].2[0], regeln[re].2[1])].exists()
            {
                regeln.remove(re);
                anz_remov += 1;
            }
        }
    }
    //println!("{}", mat);
    return mat;
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

    let pole = polplan(&v, &rigid, &erd);

    visualise(&"1.png", 720, 720, &v, &kant, &rigid, &erd, &pole);
}
