
#[derive(Debug)]
struct Point {
    num: usize,
    name: String,
    x: f64,
    y: f64
}

fn main() {
    println!("Test");
    
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
    for i in 0..4 {
        v.push(
            Point {
                num: i+15,
                name: format!("B{}",i+1),
                x: (3.0 * i as f64),
                y: 0.0
            }
        );
    }
    // Kanten 
    
    println!("{:?}",v);
}