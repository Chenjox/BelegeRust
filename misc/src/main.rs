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
    let mut rgb3 = (0,0,0);

    rgb3.0 = ((rgb1.0 as f64) *(1.0 - mixing) + (rgb2.0 as f64) *(mixing)) as u8;
    rgb3.1 = ((rgb1.1 as f64) *(1.0 - mixing) + (rgb2.1 as f64) *(mixing)) as u8;
    rgb3.2 = ((rgb1.2 as f64) *(1.0 - mixing) + (rgb2.2 as f64) *(mixing)) as u8;
    let alpha = color1.alpha() *(1.0 - mixing) + color2.alpha() * 1.0;
    return RGBAColor(rgb3.0, rgb3.1, rgb3.2, alpha)
}

const OUT_FILE_NAME: &'static str = "misc/earthquake.gif";
pub fn main() {
    let size = 100;
    let root_gif = BitMapBackend::gif(OUT_FILE_NAME, (100, 100), 100)
        .unwrap()
        .into_drawing_area()
        .apply_coord_spec(Cartesian2d::<RangedCoordusize, RangedCoordusize>::new(
            0..size,
            0..size,
            (0..100, 0..100),
        ));
    // Infinite Square grid
    let mut rng = thread_rng();
    let mut grid = DynMatrix::zeros(size, size);
    let mut update_grid = DynMatrix::zeros(size, size);

    fill_random(&mut rng, &mut grid);

    let alpha = 0.24;

    let mut t = 0.0;
    let tstart = 500.0;
    let tdelta = 5.0;
    let mut avalanche = 0;
    let mut history: Vec<[f64; 2]> = Vec::new();
    loop {
        if t > tstart + tdelta {
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
            avalanche = avalanche +1;
            // Reset Matrix
            update_grid.fill(0.0);
        } else {
            for i in 0..size {
                for j in 0..size {
                    grid[(i, j)] = grid[(i, j)] + (1.0 - max);
                }
            }
            t = t + (1.0 - max);
            history.push([t, avalanche as f64]);
            println!("t = {}, a = {}", t, avalanche);
            avalanche = 0;
            if t > tstart {
                root_gif.fill(&WHITE).unwrap();
                for i in 0..size {
                    for j in 0..size {
                        root_gif
                            .draw(&Circle::new(
                                (i, j),
                                1,
                                Into::<ShapeStyle>::into(get_interpolate_color(&BLUE,&RED, grid[(i, j)])).filled(),
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
}
