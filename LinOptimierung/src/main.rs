mod lineare_optimierung;

use crate::lineare_optimierung::*;

fn belegaufgabe() {
    let zf = LineareZielfunktion {
        koeffizienten: vec![-2.0, 2.0, 1.0, 1.0],
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
            zielwert: 4.0,
            ungleichung: true,
            koeffizienten: vec![1.0, -1.0, 0.0, 1.0],
        },
    ];

    let erg = linear_optimize(&zf, &nbs);
    println!("{:?}",erg);
}

fn uebungsvariable() {
    let zf = LineareZielfunktion {
        koeffizienten: vec![2.0, 3.0],
    };
    let nbs = vec![
        Nebenbedingung {
            zielwert: -5.0,
            ungleichung: true,
            koeffizienten: vec![-1.0, -1.0],
        },
        Nebenbedingung {
            zielwert: 5.0,
            ungleichung: true,
            koeffizienten: vec![-1.0, 1.0],
        },
        Nebenbedingung {
            zielwert: -2.0,
            ungleichung: true,
            koeffizienten: vec![2.0, -1.0],
        },
    ];

    let erg = linear_optimize(&zf, &nbs);
    println!("{:?}",erg);
}

fn main() {
    belegaufgabe();
}
