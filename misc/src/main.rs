use nalgebra::{Dynamic, OMatrix, SMatrix, SVector};
use plotters::{coord::types::RangedCoordusize, prelude::*, style::colors};
use rand::prelude::*;

type DynMatrix = OMatrix<f64, Dynamic, Dynamic>;

fn fill_random<R: Rng + ?Sized>(rng: &mut R, mat: &mut DynMatrix) {
  for i in 0..mat.shape().0 {
    for j in 0..mat.shape().1 {
      mat[(i, j)] = rng.gen();
    }
  }
}

fn get_interpolate_color(color1: &RGBColor, color2: &RGBColor, mixing: f64) -> RGBAColor {
  let rgb1 = color1.rgb();
  let rgb2 = color2.rgb();
  let mut rgb3 = (0, 0, 0);

  rgb3.0 = ((rgb1.0 as f64) * (1.0 - mixing) + (rgb2.0 as f64) * (mixing)) as u8;
  rgb3.1 = ((rgb1.1 as f64) * (1.0 - mixing) + (rgb2.1 as f64) * (mixing)) as u8;
  rgb3.2 = ((rgb1.2 as f64) * (1.0 - mixing) + (rgb2.2 as f64) * (mixing)) as u8;
  let alpha = color1.alpha() * (1.0 - mixing) + color2.alpha() * 1.0;
  return RGBAColor(rgb3.0, rgb3.1, rgb3.2, alpha);
}

const OUT_FILE_NAME: &'static str = "misc/earthquake.gif";
const OUT_FILE_NAME_2: &'static str = "misc/earthquake_magni.png";
pub fn main() {
  let size = 75;
  let root_gif = BitMapBackend::gif(OUT_FILE_NAME, (size, size), 100)
    .unwrap()
    .into_drawing_area()
    .apply_coord_spec(Cartesian2d::<RangedCoordusize, RangedCoordusize>::new(
      0..(size as usize),
      0..(size as usize),
      (0..size as i32, 0..size as i32),
    ));
  let size = size as usize;
  // Infinite Square grid
  let mut rng = thread_rng();
  let mut grid = DynMatrix::zeros(size, size);
  let mut update_grid = DynMatrix::zeros(size, size);

  fill_random(&mut rng, &mut grid);

  let alpha = 0.09;

  let mut t = 0.0;
  let tstart = 100.0;
  let tdelta = 3.0;
  let mut iterations: usize = 0;
  let max_iteration = 1_000_000;
  let mut avalanche = 0;
  let mut history: Vec<[f64; 2]> = Vec::new();
  loop {
    if t > tstart + tdelta || iterations > max_iteration {
      break;
    }
    let max = grid.max();
    if max >= 1.0 {
      // Maximum größer 1!s
      for i in 0..size {
        for j in 0..size {
          if grid[(i, j)] >= 1.0 {
            update_grid[(i, j)] = -grid[(i, j)];
            for ii in [-1_isize, 1] {
              match update_grid.get_mut(((i as isize + ii) as usize, j as usize)) {
                Some(t) => *t += grid[(i, j)] * alpha,
                None => {}
              }
              match update_grid.get_mut((i as usize, (j as isize + ii) as usize)) {
                Some(t) => *t += grid[(i, j)] * alpha,
                None => {}
              }
            }
          }
        }
      } // Update Grid
      grid = grid + &update_grid;
      avalanche = avalanche + 1;
      // Reset Matrix
      update_grid.fill(0.0);
    } else {
      iterations += 1;
      for i in 0..size {
        for j in 0..size {
          grid[(i, j)] = grid[(i, j)] + (1.0 - max);
        }
      }
      t = t + (1.0 - max);
      history.push([t, avalanche as f64]);
      if iterations % 100 == 0 {
        println!("{},t = {}, a = {}", iterations, t, avalanche);
      }
      avalanche = 0;
      if t > tstart {
        root_gif.fill(&WHITE).unwrap();
        for i in 0..size {
          for j in 0..size {
            root_gif
              .draw(&Circle::new(
                (i, j),
                1,
                Into::<ShapeStyle>::into(get_interpolate_color(&BLUE, &RED, grid[(i, j)])).filled(),
              ))
              .unwrap();
          }
        }
        root_gif.present().unwrap();
      }
    }
    //println!("{}",grid);
    // visual
  }
  let history = history;
  let history_avg = history
    .windows(10)
    .map(|w| {
      w.iter()
        .fold([0.0f64; 2], |a, e| [a[0] + e[0], a[1] + e[1]])
    })
    .map(|w| [w[0] / 10.0, w[1] / 10.0])
    .collect::<Vec<[f64; 2]>>();
  let delta_time = history
    .iter()
    .map(|f| f[0])
    .collect::<Vec<f64>>()
    .windows(2)
    .map(|f| f[1] - f[0])
    .collect::<Vec<f64>>();
  let len = history.len();

  let max_t = history.last().unwrap()[0];
  let max_avalanche = history.iter().map(|e| e[1]).reduce(f64::max).unwrap();
  let root = BitMapBackend::new(OUT_FILE_NAME_2, (2 * len as u32, 768)).into_drawing_area();
  root.fill(&WHITE).unwrap();
  let mut chart = ChartBuilder::on(&root)
    .margin(5)
    .build_cartesian_2d(0f64..max_t, 0f64..(max_avalanche + 5.0))
    .unwrap();

  chart
    .draw_series(LineSeries::new(history.iter().map(|f| (f[0], f[1])), RED))
    .unwrap();
  chart
    .draw_series(LineSeries::new(
      history_avg.iter().map(|f| (f[0], f[1])),
      BLUE,
    ))
    .unwrap();
  root.present().unwrap();
}
