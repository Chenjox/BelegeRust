use numerals::roman::Roman;
use plotters::coord::types::RangedCoordf64;
use plotters::prelude::full_palette::{
    AMBER_900, BLUEGREY_800, BLUE_900, BROWN_900, CYAN_900, DEEPORANGE_900, LIME_900, PURPLE_900,
    TEAL_900, YELLOW_900,
};
use plotters::prelude::*;

use crate::{sortPoints, Beam, DAdjUsize, DCoordMat, Mapping, Point, RigidBody};

type TPoint = (f64, f64);

// Is the turn counter clockwise?
fn turn_counter_clockwise(p1: TPoint, p2: TPoint, p3: TPoint) -> bool {
    (p3.1 - p1.1) * (p2.0 - p1.0) >= (p2.1 - p1.1) * (p3.0 - p1.0)
}

fn jarvis_march(gift: &[TPoint]) -> Option<Vec<TPoint>> {
    // There can only be a convex hull if there are more than 2 points
    if gift.len() < 3 {
        return None;
    }

    let leftmost_point = gift
        // Iterate over all points
        .iter()
        // Find the point with minimum x
        .min_by(|i, i2| i.0.total_cmp(&i2.0))
        // If there are no points in the gift, there might
        // not be a minimum. Unwrap fails (panics) the program
        // if there wasn't a minimum, but we know there always
        // is because we checked the size of the gift.
        .unwrap()
        .clone();

    let mut hull = vec![leftmost_point];

    let mut point_on_hull = leftmost_point;
    loop {
        // Search for the next point on the hull
        let mut endpoint = gift[0];
        for i in 1..gift.len() {
            if endpoint == point_on_hull
                || !turn_counter_clockwise(gift[i], hull[hull.len() - 1], endpoint)
            {
                endpoint = gift[i];
            }
        }

        point_on_hull = endpoint;

        // Stop whenever we got back to the same point
        // as we started with, and we wrapped the gift
        // completely.
        if hull[0] == endpoint {
            break;
        } else {
            hull.push(point_on_hull);
        }
    }

    Some(hull)
}

fn centroid(vec: &Vec<TPoint>) -> TPoint {
    let mut cen = (0.0, 0.0);
    for i in vec {
        cen.0 += i.0;
        cen.1 += i.1;
    }
    cen.0 = cen.0 / (vec.len() as f64);
    cen.1 = cen.1 / (vec.len() as f64);
    return cen;
}

static PALETTE: [RGBColor; 10] = [
    AMBER_900,
    BLUE_900,
    BROWN_900,
    CYAN_900,
    DEEPORANGE_900,
    LIME_900,
    PURPLE_900,
    TEAL_900,
    YELLOW_900,
    BLUEGREY_800,
];

fn draw_points<DB: DrawingBackend>(
    drawing_area: &DrawingArea<DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    points: &Vec<Point>,
) {
    for i in points {
        drawing_area
            .draw(
                &(EmptyElement::at((i.x, i.y))
                    + Circle::new((0, 0), 2, Into::<ShapeStyle>::into(&BLACK))),
            )
            .unwrap();
    }
}
fn draw_beams<DB: DrawingBackend>(
    drawing_area: &DrawingArea<DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    points: &Vec<Point>,
    beams: &Vec<Beam>,
) {
    let mapping = Mapping {
        map: sortPoints(points),
    };
    for i in beams {
        let from = i.from;
        let to = i.to;
        let mut line_points = Vec::new();
        {
            let po = &points[mapping.map_from(from)];
            line_points.push((po.x, po.y));
        }
        {
            let po = &points[mapping.map_from(to)];
            line_points.push((po.x, po.y));
        }
        drawing_area
            .draw(&Polygon::new(line_points, &BLACK))
            .unwrap();
    }
}
fn draw_rigid_bodies<DB: DrawingBackend>(
    drawing_area: &DrawingArea<DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    points: &Vec<Point>,
    rigidbodies: &Vec<RigidBody>,
    erdscheibe: &RigidBody,
) {
    let mut rigid = rigidbodies.clone();
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
    rigid.remove(end);
    //
    // hier ist keine Erdscheibe mehr drin
    let rigidbodies = rigid;
    let mapping = Mapping {
        map: sortPoints(points),
    };
    let mut palette_iter = 0;
    for body in rigidbodies {
        let mut gift = Vec::new();
        for i in &body.points {
            let p = &points[mapping.map_from(*i)];
            gift.push((p.x, p.y));
        }
        let mut line_points = Vec::new();
        if let Some(hull) = jarvis_march(&gift) {
            for p in &hull {
                line_points.push((p.0, p.1));
            }
            line_points.push(hull[0]);
        } else {
            for p in &body.points {
                let p = &points[mapping.map_from(*p)];
                line_points.push((p.x, p.y))
            }
        };
        let centroid = centroid(&line_points);
        drawing_area
            .draw(&Polygon::new(
                line_points.clone(),
                PALETTE[palette_iter].mix(0.2),
            ))
            .unwrap();
        drawing_area
            .draw(&PathElement::new(line_points, PALETTE[palette_iter]))
            .unwrap();

        drawing_area
            .draw(
                &(EmptyElement::at((centroid.0, centroid.1))
                    + Circle::new((0, 0), 3, ShapeStyle::from(&PALETTE[palette_iter]))
                    + Text::new(
                        format!("({:X})", Roman::from((palette_iter + 1) as i16)),
                        (10, 0),
                        ("sans-serif", 15.0).into_font(),
                    )),
            )
            .unwrap();
        palette_iter = (palette_iter + 1) % (10);
    }
}

fn draw_pole<DB: DrawingBackend>(
    drawing_area: &DrawingArea<DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    pole: &DCoordMat,
) {
    let mut visited_coords = Vec::new();
    let mut palette_iter = 0;
    for i in 0..pole.shape().0 {
        for j in i..pole.shape().1 {
            let (x,y) = pole[(i, j)].get_real_coordinates();
            let infty = pole[(i, j)].is_at_infinity();
            if i == j && !infty {
                drawing_area
                    .draw(
                        &(EmptyElement::at((x, y))
                            + Circle::new(
                                (0, 0),
                                5,
                                ShapeStyle::from(&PALETTE[palette_iter]).filled(),
                            )
                            + Text::new(
                                format!("({})", i + 1),
                                (10, 0),
                                ("sans-serif", 15.0).into_font(),
                            )),
                    )
                    .unwrap();
            } else if !infty {
                let mut offset = 0;

                offset += visited_coords
                    .iter()
                    .filter(|&n: &&(f64, f64)| {
                        ((n.0 - x).powi(2) + (n.1 - y).powi(2)).sqrt() < 16e-14
                    })
                    .count();

                visited_coords.push((x, y));
                drawing_area
                    .draw(
                        &(EmptyElement::at((x, y))
                            + Circle::new((0, 0), 3, ShapeStyle::from(&BLUE).filled())
                            + Text::new(
                                format!("({},{})", i + 1, j + 1),
                                (-20 as i32, -20 - 20 * offset as i32),
                                ("sans-serif", 15.0).into_font(),
                            )),
                    )
                    .unwrap();
            }
            if i == j {
                palette_iter = (palette_iter + 1) % (10);
            }
        }
    }
}

pub fn visualise(
    path: &str,
    res_x: u32,
    res_y: u32,
    points: &Vec<Point>,
    beams: &Vec<Beam>,
    rigidbodies: &Vec<RigidBody>,
    erdscheibe: &RigidBody,
    polplan: &DCoordMat,
) {
    // Ein paar diagnosti
    let mut max_x: f64 = 0.0;
    let mut min_x: f64 = 0.0;
    let mut max_y: f64 = 0.0;
    let mut min_y: f64 = 0.0;
    for i in points {
        max_x = max_x.max(i.x);
        min_x = min_x.min(i.x);
        max_y = max_y.max(i.y);
        min_y = min_y.min(i.y);
    }
    for i in polplan.iter() {
        if !i.is_at_infinity() {
            max_x = max_x.max(i.x);
            min_x = min_x.min(i.x);
            max_y = max_y.max(i.y);
            min_y = min_y.min(i.y);
        }
    }
    //
    let max = max_x.max(max_y);
    let min = min_x.min(min_y);
    let margin = 40;
    let root = BitMapBackend::new(path, (res_x, res_y))
        .into_drawing_area()
        .apply_coord_spec(Cartesian2d::<RangedCoordf64, RangedCoordf64>::new(
            min..max,
            max..min,
            (40..(res_x - margin) as i32, 40..(res_y - margin) as i32),
        ));
    root.fill(&WHITE).unwrap();

    draw_points(&root, points);
    draw_beams(&root, points, beams);
    draw_rigid_bodies(&root, points, rigidbodies, erdscheibe);
    draw_pole(&root, &polplan)

    // And if we want SVG backend
    // let backend = SVGBackend::new("output.svg", (800, 600));
}
