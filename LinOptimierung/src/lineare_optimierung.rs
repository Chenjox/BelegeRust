use nalgebra::{Dynamic, OMatrix};

type dynMatrix = OMatrix<f64, Dynamic, Dynamic>;

pub struct Nebenbedingung {
    pub zielwert: f64,
    pub ungleichung: bool, // Stets gleich oder kleiner gleich dem zielwert
    pub koeffizienten: Vec<f64>,
}

pub struct LineareZielfunktion {
    pub koeffizienten: Vec<f64>,
}

pub fn linear_optimize(zf: &LineareZielfunktion, nb: &Vec<Nebenbedingung>) -> Vec<f64> {
    // Umformen der Nebenbedingungen
    // First pass für die Anzahl der Ungleichungsnebenbedingungen und gleichungsbedingungen
    let mut ungleich: usize = 0;
    let mut gleich: usize = 0;
    for elem in nb {
        if elem.ungleichung {
            ungleich += 1;
        } else {
            gleich += 1;
        }
    }
    // ab hier ändern sich die Variablen nicht mehr
    // let mut -> let
    let unabhaengige = zf.koeffizienten.len();
    let ungleich = ungleich;
    let gleich = gleich;
    // Second pass für die Einführung von Schlupfvariablen
    // Speichern der Indizes der gleichungs und ungleichungsnb in `nb` für nachfolgende Sortierung.
    let mut gleich_indizes = vec![0; gleich.try_into().unwrap()];
    let mut ungleich_indizes = vec![0; ungleich.try_into().unwrap()];
    let mut gleich_index = 0;
    let mut ungleich_index = 0;
    for index in 0..nb.len() {
        if nb[index].ungleichung {
            ungleich_indizes[ungleich_index] = index;
            ungleich_index += 1;
        } else {
            gleich_indizes[gleich_index] = index;
            gleich_index += 1;
        }
    }
    // Kontrolle gleich_index = gleich
    // Kontrolle ungleich_index = ungleich

    // Umsortieren. in der Letzten Spalte steht die VZ RHS
    let mut nb_koeffizienten = dynMatrix::zeros((gleich + ungleich) as usize, unabhaengige + 1);
    for i in 0..nb.len() {
        for j in 0..unabhaengige {
            if i < gleich {
                nb_koeffizienten[(i, j)] = nb[gleich_indizes[i]].koeffizienten[j];
                nb_koeffizienten[(i, unabhaengige)] = nb[gleich_indizes[i]].zielwert;
            } else {
                nb_koeffizienten[(i, j)] = nb[ungleich_indizes[i - gleich]].koeffizienten[j];
                nb_koeffizienten[(i, unabhaengige)] = nb[ungleich_indizes[i]].zielwert;
            }
        }
    }
    // NB Koeffizienten stehen alle Koeffizienten der NB drin
    // zuerst Koeffizienten von Gleichungsnebenbedinungen
    // dann Koeffizienten von Ungleichungsnebenbedinungen

    // Wie viele Schlupfvariablen gibt es? für jede ungleichungsnebenbedingung eine!
    // + anzahl der koeffizienten der Zielfunktion + Zielfunktionswert
    // anzahl der Spalten des Tablaux
    let n_vars = ungleich + &zf.koeffizienten.len() + 1;
    // Bau der Matrix:
    // Erste Zeile ist Zielfunktion
    // danach n Gleichungsbedinungen
    // danach m ungleichungsnebenbedingungen mit Schlupf
    // = n + m + 1 Zeilen
    let n_zeilen = ungleich + gleich + 1;
    let mut tablaux = dynMatrix::zeros(n_zeilen, n_vars);

    // Einsortieren der Werte
    // Erste Zeile ist die Zielfunktion
    for j in 0..(unabhaengige) {
        tablaux[(0, j)] = zf.koeffizienten[j];
    }
    // Jetzt die gleichungs und ungleichsnebenbedinungen

    for i in 0..(n_zeilen - 1) {
        // zielwert in die letzte spalte
        let sign = nb_koeffizienten[(i, unabhaengige)].signum();
        tablaux[(i + 1, n_vars - 1)] = sign * nb_koeffizienten[(i, unabhaengige)];
        for j in 0..(unabhaengige) {
            // da die erste Zeile die Zielfunktion ist, muss der index i verschoben werden.
            tablaux[(i + 1, j)] = sign * nb_koeffizienten[(i, j)];
        }
    }
    for i in 0..(ungleich) {
        let sign = nb_koeffizienten[(i, unabhaengige)].signum();
        tablaux[(1 + gleich + i, unabhaengige + i)] = 1.0 * sign;
    }
    let MAX_ITER = 1000;
    let mut iter = 0;
    let mut erg_indizes = vec![usize::MAX; ungleich + &zf.koeffizienten.len()];
    // JETZT beginnt der Simplex Algorithmus //
    'simplex: loop {
        iter += 1;
        println!("{}", tablaux);
        // Auswählen der Pivotspalte
        // Größter Koeffizient der ZF ist Pivotspalte
        let mut pivotcol = 0;
        let mut max = 0.0;
        for i in 0..n_vars - 1 {
            let current = tablaux[(0, i)];
            if current > max {
                max = current;
                pivotcol = i;
            }
        }
        // Überprüfen der Positivität des Maximums
        if max <= 0.0 || iter > MAX_ITER {
            break 'simplex;
        }

        // Bestimmen der Pivotzeile
        let mut pivotrow = 1;
        let mut min = f64::INFINITY; // init mit "+unendlich"
        for i in 1..n_zeilen {
            let quotient = tablaux[(i, n_vars - 1)] / tablaux[(i, pivotcol)];
            if quotient >= 0.0 && quotient < min {
                pivotrow = i;
                min = quotient;
            }
        }
        let pivot_val = tablaux[(pivotrow, pivotcol)];

        erg_indizes[pivotcol] = pivotrow;
        // Bestimmen der neuen Ergebnisse
        // Normierung der Pivotzeile
        for j in 0..(n_vars) {
            tablaux[(pivotrow, j)] = tablaux[(pivotrow, j)] / pivot_val;
        }
        // RHS
        for i in 1..n_zeilen {
            if i != pivotrow {
                tablaux[(i, n_vars - 1)] = tablaux[(i, n_vars - 1)]
                    - tablaux[(i, pivotcol)] * tablaux[(pivotrow, n_vars - 1)];
            }
        }
        // Alle Element außer der Pivotzeilen + Spalten
        for i in 0..n_zeilen {
            for j in 0..n_vars - 1 {
                if i != pivotrow && j != pivotcol {
                    tablaux[(i, j)] =
                        tablaux[(i, j)] - tablaux[(i, pivotcol)] * tablaux[(pivotrow, j)];
                }
            }
        }
        //Zielfunktionsergebnis
        tablaux[(0, n_vars - 1)] =
            tablaux[(0, n_vars - 1)] - tablaux[(0, pivotcol)] * tablaux[(pivotrow, n_vars - 1)];
        // Alle Pivotspaltenelemente außer Pivotelement selbst
        for i in 0..n_zeilen {
            if i != pivotrow {
                tablaux[(i, pivotcol)] = 0.0;
            }
        }
    }
    //println!("{}", tablaux);
    //println!("{:?}", erg_indizes);

    let mut res = Vec::new();
    for i in 0..zf.koeffizienten.len() {
        if erg_indizes[i] != usize::MAX {
            res.push(tablaux[(erg_indizes[i], n_vars - 1)]);
        } else {
            res.push(0.0);
        }
    }

    return res;
}
