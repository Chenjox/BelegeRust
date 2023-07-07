
use plotters::prelude::*;

fn plot_moment_to_file(
    path: &str,
    einzellast: f64,
    pos_last: f64,
    laenge: f64
) {
    let root = BitMapBackend::new(path, (640, 480)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut chart = ChartBuilder::on(&root)
        .caption("Momentenverlauf", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0f64..laenge, -200f64..200f64)
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    chart
        .draw_series(LineSeries::new(
            (0..=1000).map(|x| x as f64 / 1000.0 * laenge).map(|x| {
                (
                    x,
                    //querkraft_bei_einzellast(einzellast, pos_last,laenge, x)
                    LF_zwei_einzellasten_querkraft(x,laenge,pos_last,laenge-pos_last,einzellast),
                )
            }),
            &RED,
        ))
        .unwrap()
        .label("Momentenverlauf")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .unwrap();

    root.present();
}

// Auflagerkräfte links und rechts
fn auflagerkraefte(einzellast: f64, pos_last: f64, laenge: f64) -> [f64;2] {
    let left_druck = einzellast * ((laenge - pos_last)/laenge).powi(2) * (1.0 + (2.0 * pos_last)/laenge);
    let right_druck = einzellast * (pos_last/laenge).powi(2) * (1.0 + (2.0 * (laenge -pos_last))/laenge);
    return [left_druck,right_druck];
}

// Auflagermomente links und rechts
fn auflagermomente(einzellast: f64, pos_last: f64, laenge: f64) -> [f64;2] {
    let left_mom = - einzellast * pos_last * ((laenge-pos_last)/laenge).powi(2);
    let right_mom = - einzellast * (laenge - pos_last) * (pos_last/laenge).powi(2);
    return [left_mom,right_mom];
}

// Auflagerkräfte infolge Linienlast, links und rechts
fn auflagerkraefte_lineload(linienlast: f64, laenge: f64) -> [f64;2] {
    let auflager = linienlast * laenge / 2.0;
    return [auflager,auflager];
}

fn auflagermomente_lineload(linienlast: f64, laenge: f64) -> [f64;2] {
    let auflager = - linienlast * laenge*laenge/12.0;
    return [auflager,auflager];
}

fn moment_bei_einzellast(einzellast: f64, pos_last: f64, laenge: f64, x_wert: f64) -> f64 {
    let mom = auflagermomente(einzellast, pos_last, laenge);
    let druck = auflagerkraefte(einzellast, pos_last, laenge);
    if x_wert < pos_last {
        return druck[0] * x_wert + mom[0];
    } else {
        return druck[1] * (laenge- x_wert) + mom[1];
    }
}

fn querkraft_bei_einzellast(einzellast: f64, pos_last: f64, laenge: f64, x_wert: f64) -> f64 {
    let druck = auflagerkraefte(einzellast, pos_last, laenge);
    if x_wert < pos_last {
        return druck[0];
    } else {
        return druck[0] - einzellast;
    }
}

fn moment_bei_linienlast(linienlast: f64, laenge: f64, x_wert: f64) -> f64 {
    let mom = auflagermomente_lineload(linienlast,  laenge);
    let druck = auflagerkraefte_lineload(linienlast,  laenge);
    return druck[0] * x_wert + mom[0] - linienlast * x_wert * x_wert/2.0;
}

fn LF_zwei_einzellasten_moment(x_wert: f64, laenge: f64 ,pos_last_1: f64, pos_last_2: f64, einzellast: f64)-> f64 {
    return 
    moment_bei_einzellast(einzellast, pos_last_1, laenge, x_wert)
    +moment_bei_einzellast(einzellast, pos_last_2, laenge, x_wert);
}
fn LF_zwei_einzellasten_querkraft(x_wert: f64, laenge: f64 ,pos_last_1: f64, pos_last_2: f64, einzellast: f64)-> f64 {
    return 
    querkraft_bei_einzellast(einzellast, pos_last_1, laenge, x_wert)
    +querkraft_bei_einzellast(einzellast, pos_last_2, laenge, x_wert);
}


fn nachweis_stirnplattenstoss(
    querkraft: f64, 
    moment: f64, 
    blechdicke: f64,
    blechlaenge: f64,
    schraubendurchmesser: f64, 
    anzahl_schrauben: i32
) {
// TODO ! 
}

fn nachweis_biegedrillknicken() {
    
}

fn main() {
    let laenge = 4.4; // m
    let stuetzenlast = 75.0; // Kn
    let holzlast = 0.192 * 0.7; // kN
    let einzellast = stuetzenlast + holzlast; // kN 
    let pos_last = 4.4 - 4.4/2.0 - 0.15 - 0.08; // m
    let eigengewicht_traeger = 42.5*9.81/1000.0; // In kN/m

    let emodul = 210000.0;
    let festigkeit = 235.0 /10.0; // kM/cm^2 

    let drittelspunkt = 4.4/3.0;
    let LF1_moment_stirn=     LF_zwei_einzellasten_moment(drittelspunkt,laenge, pos_last, laenge - pos_last, einzellast);
    let LF1_moment_mitte =       LF_zwei_einzellasten_moment(laenge/2.0   ,laenge, pos_last, laenge - pos_last, einzellast);
    let LF1_querkraft_mitte = LF_zwei_einzellasten_querkraft(laenge/2.0   ,laenge, pos_last, laenge - pos_last, einzellast);
    let LF1_querkraft_stirn = LF_zwei_einzellasten_querkraft(drittelspunkt, laenge, pos_last, laenge - pos_last, einzellast);
    //let LF1_querkraft = querkraft_bei_einzellast(einzellast, pos_last, laenge, drittelspunkt);

    println!("{},{},{},{}",LF1_moment_mitte,LF1_querkraft_mitte,LF1_moment_stirn,LF1_querkraft_stirn);


    let welast = LF1_moment_mitte*100.0/festigkeit * 1.1;
    println!("Benötigtes elastisches Widerstandmoment: \n W_el = {}",welast);



    plot_moment_to_file("test.png", einzellast, pos_last, laenge)
}
