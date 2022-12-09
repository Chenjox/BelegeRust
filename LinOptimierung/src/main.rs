mod lineare_optimierung;

use crate::lineare_optimierung::*;

fn main() {
    let zf = LineareZielfunktion {
        koeffizienten: vec![2.0, -2.0, -1.0, -1.0],
    };
    let nbs = vec![
        Nebenbedingung {
            zielwert: 1.0,
            ungleichung: true,
            koeffizienten: vec![-1.0, 2.0, 1.0, -1.0],
        },
        Nebenbedingung {
            zielwert: 1.0,
            ungleichung: true,
            koeffizienten: vec![1.0, 4.0, 0.0, -2.0],
        },
        Nebenbedingung {
            zielwert: 1.0,
            ungleichung: true,
            koeffizienten: vec![1.0, -1.0, 0.0, 1.0],
        },
    ];

    linear_optimize(&zf, &nbs);
    println!("Hello World");
}
