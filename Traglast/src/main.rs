
use nalgebra::{Dynamic, OMatrix};
use plotters::prelude::*;

type DAdjUsize = OMatrix<usize, Dynamic, Dynamic>;

#[derive(Debug,Clone)]
struct Point {
    num: usize,
    name: String,
    x: f64,
    y: f64
}
#[derive(Debug)]
struct Beam {
    from: usize,
    to: usize,
    nplast: f64,
}

fn get_tragwerk() -> (Vec<Point>,Vec<Beam>) {
    let mut v = Vec::new();

    // Knoten (noch nicht in Koordinaten die ich brauch)
    for i in 0..7 {
        v.push(
            Point {
                num: i+1,
                name: format!("A{}",i+1),
                x: (3.0 * i as f64),
                y: 7.0
            }
        );
    }
    for i in 0..7 {
        v.push(
            Point {
                num: i+8,
                name: format!("A{}",i+1),
                x: (3.0 * i as f64),
                y: 3.0
            }
        );
    }
    for i in 0..2 {
        v.push(
            Point {
                num: i+15,
                name: format!("B{}",i+1),
                x: (3.0 * i as f64),
                y: 0.0
            }
        );
        v.push(
            Point {
                num: 18-i,
                name: format!("B{}",4-i),
                x: (3.0 * (1-i) as f64 + 15.0),
                y: 0.0
            }
        );
    }
    // Kanten 
    let mut kant = Vec::new();
    for i in 0..7 {
        kant.push(
            Beam {
                from: i+1,
                to: i+8,
                nplast: 1.0,
            }
        );
    }
    for i in 0..6 {
        kant.push(
            Beam {
                from: i+1,
                to: i+2,
                nplast: 1.0,
            }
        );
        kant.push(
            Beam {
                from: i+8,
                to: i+9,
                nplast: 1.0,
            }
        );
    }
    for i in 0..3 {
        kant.push(
            Beam {
                from: i+1,
                to: i+9,
                nplast: 1.0,
            }
        );
        kant.push(
            Beam {
                from: 7-i,
                to: 13-i,
                nplast: 1.0,
            }
        );
    }
    for i in 0..2 {
        kant.push(
            Beam {
                from: i+8,
                to: i+15,
                nplast: 1.0,
            }
        );
        kant.push(
            Beam {
                from: 14-i,
                to: 18-i,
                nplast: 1.0,
            }
        );
    }
    kant.push(
        Beam {
            from: 8,
            to: 16,
            nplast: 1.0,
        }
    );
    kant.push(
        Beam {
            from: 14,
            to: 17
            ,
            nplast: 1.0,
        }
    );

    return (v,kant);
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
                    panic!("Double Indices found: {} is already given!",test);
                }
            }
            check_vec.push(test);
        }
    }
    return sortVec;
}

// erstellen der Adjazenzmatrix
fn get_adjacency_matrix(pvec: &Vec<Point>, sortVec: &Vec<usize>, bvec: &Vec<Beam>) {
    let mut adj_mat = DAdjUsize::zeros(pvec.len(),pvec.len());

    for i in 0..sortVec.len() {
        for j in 0..bvec.len() {
            let b = &bvec[j];
            if sortVec[i] == b.from {
                // index des von Punktes finden
                let f = b.to;
                for k in 0..sortVec.len() {
                    if sortVec[k] == f {
                        adj_mat[(i,k)] = 1;
                        adj_mat[(k,i)] = 1;
                    }
                }
            }
        }
    }
    println!("{}",adj_mat);

}

fn main() {
    let (v,kant) = get_tragwerk();

    let vsort = sortPoints(&v);
    
    
    println!("{:?}",v);
    println!("{:?}",vsort);
    println!("{:?}",kant);

    get_adjacency_matrix(&v, &vsort, &kant);

    {
        let res_y = 200;
        let res_x = 300;
        let margin = 10;
        let root = BitMapBackend::new("1.png", (res_x, res_y)).into_drawing_area();
        root.fill(&WHITE).unwrap();

        for i in &v {
            root.draw(&Circle::new(
              (i.x as i32 *10 + margin,(res_y as i32 - (i.y as i32*10 + margin))),
              5,
              Into::<ShapeStyle>::into(&GREEN).filled(),
            )).unwrap();
        }
        for i in &kant {
            let from = i.from;
            let to = i.to;
            let mut line_points = Vec::new();
            for po in &v {
                if po.num == from {
                    line_points.push((po.x as i32 *10 + margin as i32,(res_y as i32 - (po.y as i32*10 + margin))));
                }
                if po.num == to {
                    line_points.push((po.x as i32 *10 + margin as i32,(res_y as i32 - (po.y as i32*10 + margin))));
                }
            }
            root.draw(&Polygon::new(
              line_points,
              &BLACK,
            )).unwrap();
        }
      
      // And if we want SVG backend
      // let backend = SVGBackend::new("output.svg", (800, 600));
    }
}