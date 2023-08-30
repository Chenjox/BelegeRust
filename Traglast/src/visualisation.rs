use std::f64::consts::PI;

use nalgebra::{Matrix1, Matrix2, Vector2};
use numerals::roman::Roman;
use plotters::coord::types::RangedCoordf64;
use plotters::prelude::full_palette::{
  AMBER_900, BLUEGREY_800, BLUE_900, BROWN_900, CYAN_900, DEEPORANGE_900, LIME_900, PURPLE_900,
  TEAL_900, YELLOW_900,
};
use plotters::prelude::*;

use crate::{Beam2DMatrix, Load2DFlattenMatrix, Point2DMatrix};

type S2x2 = Matrix2<f64>;
type S2 = Vector2<f64>;
type TPoint = (f64, f64);

//
fn get_rot_matrix(angle: f64) -> S2x2 {
  return S2x2::new(angle.cos(), -angle.sin(), angle.sin(), angle.cos());
}

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

pub fn visualise(
  path: &str,
  res_x: u32,
  res_y: u32,
  points: &Point2DMatrix,
  beams: &Beam2DMatrix,
  loads: &Load2DFlattenMatrix,
) {
  // Ein paar diagnosti
  let mut max_x: f64 = f64::MIN;
  let mut min_x: f64 = f64::MAX;
  let mut max_y: f64 = f64::MIN;
  let mut min_y: f64 = f64::MAX;
  let num_points = points.ncols();
  for i in 0..num_points {
    max_x = max_x.max(points[(0, i)]);
    min_x = min_x.min(points[(0, i)]);
    max_y = max_y.max(points[(1, i)]);
    min_y = min_y.min(points[(1, i)]);
  }
  //for i in polplan.iter() {
  //    if !i.is_at_infinity() {
  //        max_x = max_x.max(i.x);
  //        min_x = min_x.min(i.x);
  //        max_y = max_y.max(i.y);
  //        min_y = min_y.min(i.y);
  //    }
  //}
  //
  let max = max_x.max(max_y);
  let min = min_x.min(min_y);
  let margin = 40;
  let root = BitMapBackend::new(path, (res_x, res_y))
    .into_drawing_area()
    .apply_coord_spec(Cartesian2d::<RangedCoordf64, RangedCoordf64>::new(
      min..max,
      max..min,
      //((0..res_x as i32),(0..res_y as i32))
      (
        (margin as i32)..(res_x as i32 - margin),
        (margin as i32)..(res_y as i32 - margin),
      ),
    ));
  root.fill(&WHITE).unwrap();

  //draw_box(max_x, min_x, max_y, min_y, &root);
  draw_coordinate_system(&root);
  draw_points(&root, &points);
  draw_beams(&root, points, beams);
  draw_loads(&root, &loads, &points)
  // And if we want SVG backend
  // let backend = SVGBackend::new("output.svg", (800, 600));
}

fn draw_loads<DB: DrawingBackend>(
  drawing_area: &DrawingArea<DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
  loads: &Load2DFlattenMatrix,
  points: &Point2DMatrix,
) {
  //println!("{:?},{:?}",points.shape(),loads.shape());

  for (i, point) in points.column_iter().enumerate() {
    let coords = point;
    let loadvec = { S2::new(loads[2 * i], loads[2 * i + 1]) }.normalize();
    let begin_load = coords - loadvec * 1.1;
    let end_load = coords - loadvec * 0.1;
    let arrow1 = end_load - get_rot_matrix(PI / 4.0) * loadvec * 0.2;
    let arrow2 = end_load - get_rot_matrix(-PI / 4.0) * loadvec * 0.2;
    drawing_area
      .draw(&PathElement::new(
        vec![
          (begin_load.x, begin_load.y),
          (end_load.x, end_load.y),
          (arrow1.x, arrow1.y),
          (arrow2.x, arrow2.y),
          (end_load.x, end_load.y),
        ],
        &RED,
      ))
      .unwrap();
  }
}

fn draw_box<DB: DrawingBackend>(
  max_x: f64,
  min_x: f64,
  max_y: f64,
  min_y: f64,
  drawing_area: &DrawingArea<DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
) {
  drawing_area
    .draw(&PathElement::new(
      vec![
        (min_x, min_y),
        (min_x, max_y),
        (max_x, max_y),
        (max_x, min_y),
        (min_x, min_y),
      ],
      &RED,
    ))
    .unwrap();
}

fn draw_coordinate_system<DB: DrawingBackend>(
  drawing_area: &DrawingArea<DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
) {
  drawing_area
    .draw(&Polygon::new(vec![(0., 0.), (1., 0.)], &BLACK))
    .unwrap();
  drawing_area
    .draw(&Polygon::new(vec![(0., 0.), (0., 1.)], &BLACK))
    .unwrap();
}

fn draw_points<DB: DrawingBackend>(
  drawing_area: &DrawingArea<DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
  points: &Point2DMatrix,
) {
  let mut num = 0;
  for i in points.column_iter() {
    drawing_area
      .draw(
        &(EmptyElement::at((i.x, i.y))
          + Circle::new((0, 0), 2, Into::<ShapeStyle>::into(&BLACK))
          + Text::new(
            format!("{}", num),
            (10, -20),
            ("sans-serif", 12.0).into_font(),
          )),
      )
      .unwrap();
    num += 1;
  }
}

fn draw_beams<DB: DrawingBackend>(
  drawing_area: &DrawingArea<DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
  points: &Point2DMatrix,
  beams: &Beam2DMatrix,
) {
  for i in beams.column_iter() {
    let from = i.x;
    let to = i.y;
    let mut line_points = Vec::new();
    {
      let po = points.column(from);
      line_points.push((po.x, po.y));
    }
    {
      let po = points.column(to);
      line_points.push((po.x, po.y));
    }
    drawing_area
      .draw(&Polygon::new(line_points, &BLACK))
      .unwrap();
  }
}
