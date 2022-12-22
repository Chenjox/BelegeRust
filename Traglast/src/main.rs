use nalgebra::{Dynamic, OMatrix};
use plotters::prelude::*;
use std::fmt;

type DAdjUsize = OMatrix<usize, Dynamic, Dynamic>;

#[derive(Debug, Clone)]
struct Point {
    num: usize,
    name: String,
    x: f64,
    y: f64,
}

#[derive(Clone, PartialEq, Debug)]
struct Pol {
    is_at_infinity: bool,
    x: f64,
    y: f64,
}

impl fmt::Display for Pol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.x.is_nan() {
            write!(f, "()")
        } else {
            write!(f, "({}, {}, {})", self.is_at_infinity, self.x, self.y)
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
struct Beam {
    from: usize,
    to: usize,
    nplast: f64,
}

#[derive(Debug)]
struct Triangle(usize, usize, usize);

#[derive(Debug, Clone)]
struct RigidBody {
    points: Vec<usize>,
}

struct Mapping {
    map: Vec<usize>,
}

impl Mapping {
    // 0 -> 1
    fn map_to(&self, index: usize) -> usize {
        return self.map[index];
    }
    //
    fn map_from(&self, index: usize) -> usize {
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
        for i in &self.points {
            if other.points.contains(i) {
                count = count + 1;
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
    kant.swap_remove(5);
    kant.swap_remove(7);
    kant.swap_remove(9);

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
    let mut iter = 0;
    let max_iter = rigid.len() * (rigid.len() + 1) / 2;
    loop {
        iter = iter + 1;
        let mut has_changed = false;
        // alle Rigidbodys zusammenführen
        if let Some(bod) = rigid.pop() {
            // nehme den letzten RigidBody!
            for i in 0..rigid.len() {
                if let Some(bod2) = rigid.pop() {
                    // nehme den zweitletzten
                    if bod.is_joined_with(&bod2) {
                        rigid.insert(0, bod.join_with(bod2));
                        has_changed = true;
                    } else {
                        rigid.insert(0, bod2);
                    }
                }
            }
            if !has_changed {
                rigid.insert(0, bod);
            }
        }
        if !has_changed || iter > max_iter {
            break;
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
fn sortPoints(pvec: &Vec<Point>) -> Vec<usize> {
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
fn polplan(points: &Vec<Point>, bodies: &Vec<RigidBody>, erdscheibe: &RigidBody) {
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
                        let po = &points[mapping.map_from(connect_points[i].0)];
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
        for i in 0..connect_points.len() {
            let con_bodies = &connect_points[i].1;
            let po = &points[mapping.map_from(connect_points[i].0)];
            for j in 0..con_bodies.len() {
                // mindestens 2 Elemente müssen enthalten sein
                for k in j + 1..con_bodies.len() {
                    mat[(con_bodies[j], con_bodies[k])] = Pol::new(false, po.x, po.y);
                    mat[(con_bodies[k], con_bodies[j])] = Pol::new(false, po.x, po.y);
                }
            }
        }
        println!("{:?}", rigid);
        println!("{:?}", connect_points);
    }
    println!("{}", mat);
    // Matrix der Bedingungen
    let anzahl_pole = bodies.len() - 1;
    let mut mat_bed = DAdjUsize::from_element(anzahl_pole, anzahl_pole, 0);
    let mut best_kandidates = Vec::new();

    for i in 0..anzahl_pole {
        //Nebenpolbedingung 1 (HP HP)
        for j in i + 1..anzahl_pole {
            // (HP HP)
            if i != j && mat[(i, i)].exists() && mat[(j, j)].exists() {
                // Zwei Hauptpole sind bekannt
                mat_bed[(i, j)] += 1; // Wir können eine Aussage über den Nebenpol treffen
                mat_bed[(j, i)] += 1; // Wir können eine Aussage über den Nebenpol treffen
            }
        }
        for j in i + 1..anzahl_pole {
            // (NP NP)
            if i != j && mat[(i, j)].exists() {
                for k in j + 1..anzahl_pole {
                    if mat[(i, j)].exists() && mat[(i, k)].exists() {
                        mat_bed[(k, j)] += 1; // Wir können eine Aussage über den Nebenpol treffen
                        mat_bed[(j, k)] += 1;
                    }
                }
            }
        }
        for j in 0..anzahl_pole {
            //Hauptpolbedingungen (HPNP)
            if mat[(i, i)].exists() && i != j {
                // Der Hauptpol ist bekannt
                // Wir haben also einen Nebenpol
                if mat[(i, j)].exists() {
                    // Der Nebenpol ist bekannt
                    mat_bed[(j, j)] += 1; // Wir können eine Aussage über den Hauptpol treffen
                }
            }
        }
    }
    for i in 0..anzahl_pole {
        for j in i..anzahl_pole {
            if mat_bed[(i, j)] > 1 {
                best_kandidates.push((i, j));
            }
        }
    }
    // Nebenpolbedingungen ()
    println!("{}", mat_bed);
    println!("{:?}", best_kandidates);

    loop {
        // Finden des besten kandidaten
        if let Some((i, j)) = best_kandidates.pop() {
            // Ist es ein Haupt oder nebenpol?
            if i == j {
                // Also ein Hauptpol
                // Finden des ersten pols mit zulässigem Hauptpol und Nebenpol
                let mut kfinal = 0;
                let mut lfinal = 0;
                for k in 0..anzahl_pole - 1 {
                    if k < i {
                        if mat[(k, k)].exists() && mat[(i, k)].exists() {
                            kfinal = k;
                            break;
                        }
                    } else {
                        if mat[(k + 1, k + 1)].exists() && mat[(i, k + 1)].exists() {
                            kfinal = k + 1;
                            break;
                        }
                    }
                }
                for l in kfinal + 1..anzahl_pole {
                    if l < j {
                        if mat[(l, l)].exists() && mat[(j, l)].exists() {
                            lfinal = l;
                            break;
                        }
                    } else {
                        if mat[(l + 1, l + 1)].exists() && mat[(j, l + 1)].exists() {
                            lfinal = l + 1;
                            break;
                        }
                    }
                } // kfinal und lfinal gefunden!
                let k = kfinal;
                let l = lfinal;
                let new_pol = Pol::infer(&mat[(k, k)], &mat[(i, k)], &mat[(l, l)], &mat[(j, l)]);
                println!("{}", new_pol);
                mat[(i, j)] = new_pol;
                mat_bed[(i, j)] = 0; // hier gibt es nichts mehr
                                     // Ein neuer Hauptpol bei i und j!
                                     // aktualisieren der Bedingungsmatrix
                                     //Nebenpolbedingung 1 (HP HP)
                for k in 0..anzahl_pole {
                    // (HP HP)
                    if i != k && mat[(k, k)].exists() && !mat[(i, k)].exists() {
                        // HP HP -> NP
                        mat_bed[(i, k)] += 1; // Wir können eine Aussage über den Nebenpol treffen
                        mat_bed[(k, i)] += 1;
                    }
                    if i != k && mat[(k, i)].exists() && !mat[(k, k)].exists() {
                        // HP NP -> HP
                        mat_bed[(k, k)] += 1; // Wir können eine Aussage über den Hauptpol treffen
                    }
                }
            } else {
                // Es ist ein Nebenpol zu bestimmen
                // Erste Möglichkeit HP-HP-NP-NP
                if mat[(i, i)].exists() && mat[(j, j)].exists() {
                    // es existieren 2 Hauptpole
                    //finden der nebenpole
                    let mut kfinal = 0;
                    for k in 0..anzahl_pole {
                        if mat[(i, k)].exists() && mat[(j, k)].exists() {
                            // Ein NP Paar gefunden!
                            kfinal = k;
                            break;
                        }
                    }
                    let k = kfinal;
                    let new_pol =
                        Pol::infer(&mat[(i, i)], &mat[(j, j)], &mat[(i, k)], &mat[(j, k)]);
                    println!("{}", new_pol);
                    mat[(i, j)] = new_pol.clone();
                    mat[(j, i)] = new_pol;
                    mat_bed[(i, j)] = 0;
                    mat_bed[(j, i)] = 0;
                } else {
                    // zweite Möglichkeit NP-NP-NP-NP
                    // Der erste Index der erst zu inferierenden
                    let mut ktupel = (0, 0);
                    'kloop: for kj in 0..anzahl_pole {
                        if i != kj && mat[(i, kj)].exists() {
                            // potentiell erster NP des ersten Paars gefunden!
                            for ki in 0..anzahl_pole {
                                if ki != kj && mat[(ki, kj)].exists() {
                                    // erstes NP Paar gefunden!
                                    ktupel = (ki, kj);
                                    break 'kloop;
                                }
                            }
                        }
                    }
                    let mut ltupel = (0, 0);
                    'lloop: for li in 0..anzahl_pole {
                        if j != li && mat[(li, j)].exists() {
                            // potentiell erster NP des zweiten Paars gefunden!
                            for lj in 0..anzahl_pole {
                                if li != lj && mat[(li, lj)].exists() {
                                    // zweistes NP Paar gefunden!
                                    ktupel = (li, lj);
                                    break 'lloop;
                                }
                            }
                        }
                    }
                    let new_pol = Pol::infer(
                        &mat[(i, ktupel.1)],
                        &mat[ktupel],
                        &mat[(ltupel.0, j)],
                        &mat[ltupel],
                    );
                    println!("{}", new_pol);
                    mat[(i, j)] = new_pol.clone();
                    mat[(j, i)] = new_pol;
                    mat_bed[(i, j)] = 0;
                    mat_bed[(j, i)] = 0;
                }
                // Nebenpol wurde gefunden!
                // Überarbeiten der mat_bed
                // (HP NP)
                if mat[(i, i)].exists() && !mat[(j, j)].exists() {
                    // (ij) + (i) -> (j)
                    mat_bed[(j, j)] += 1;
                }
                if mat[(j, j)].exists() && !mat[(i, i)].exists() {
                    // (ij) + (i) -> (j)
                    mat_bed[(i, i)] += 1;
                }
                // (NP NP)
                for k in 0..anzahl_pole {
                    if k != i && k != j && mat[(i, k)].exists() && !mat[(j, k)].exists() {
                        // (ij) + (ik) -> GO(ik)
                        mat_bed[(j, k)] += 1; // Wir können eine Aussage über den Nebenpol treffen
                        mat_bed[(k, j)] += 1;
                    }
                    if k != i && k != j && mat[(j, k)].exists() && !mat[(i, k)].exists() {
                        // (ij) + (ik) -> GO(ik)
                        mat_bed[(i, k)] += 1; // Wir können eine Aussage über den Nebenpol treffen
                        mat_bed[(k, i)] += 1;
                    }
                }
            }
        } else {
            break;
        }
        let mut is_all_zero = true;
        for i in 0..anzahl_pole {
            for j in i..anzahl_pole {
                if mat_bed[(i, j)] > 0 {
                    is_all_zero = false;
                }
                if mat_bed[(i, j)] > 1 {
                    best_kandidates.push((i, j));
                }
            }
        }
        if is_all_zero {
            break;
        }
        //println!("{}", mat_bed);
    }
    println!("{}", mat);
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

    let mut rigid = get_rigid_bodies(&v, &kant, &erd);

    polplan(&v, &rigid, &erd);

    {
        let res_y = 200;
        let res_x = 300;
        let margin = 10;
        let root = BitMapBackend::new("1.png", (res_x, res_y)).into_drawing_area();
        root.fill(&WHITE).unwrap();

        for i in &v {
            root.draw(&Circle::new(
                (
                    i.x as i32 * 10 + margin,
                    (res_y as i32 - (i.y as i32 * 10 + margin)),
                ),
                5,
                Into::<ShapeStyle>::into(&GREEN).filled(),
            ))
            .unwrap();
        }
        for i in &kant {
            let from = i.from;
            let to = i.to;
            let mut line_points = Vec::new();
            {
                let po = &v[mapping.map_from(from)];
                line_points.push((
                    po.x as i32 * 10 + margin as i32,
                    (res_y as i32 - (po.y as i32 * 10 + margin)),
                ));
            }
            {
                let po = &v[mapping.map_from(to)];
                line_points.push((
                    po.x as i32 * 10 + margin as i32,
                    (res_y as i32 - (po.y as i32 * 10 + margin)),
                ));
            }
            root.draw(&Polygon::new(line_points, &BLACK)).unwrap();
        }

        // And if we want SVG backend
        // let backend = SVGBackend::new("output.svg", (800, 600));
    }
}
