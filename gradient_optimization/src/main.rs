use num_dual::*;
use plotters::prelude::*;
struct History {
    pub current: Vec<(f64, f64)>,
    pub next: Vec<(f64, f64)>,
}

impl History {
    fn push(&mut self, current: (f64, f64), next: (f64, f64)) {
        self.current.push(current);
        self.next.push(next);
    }
}

trait GradientOptimization2D {
    /// The target function in question.
    /// Must be differentialable at all points inside the search domain
    fn target_function<D: DualNum<f64>>(&self, x1: D, x2: D) -> D;

    /// The contraint function(s)
    /// Must return a boolean which indicates a violation of the constraints.
    fn violates_constraint_functions<D: DualNum<f64>>(&self, x1: D, x2: D) -> bool;

    /// the start point for the optimization.
    /// must not violate the contraint function
    fn start_point(&self) -> (f64, f64);

    /// A ordered list of step sizes, starting with greatest.
    /// all values must be strictly positive.
    fn steps_sizes(&self) -> &Vec<f64>;

    /// The accepted error between iterations
    fn eps(&self) -> f64;

    fn gradient_minimize(&self, hist: &mut History) -> (f64, f64) {
        let start_point = self.start_point();
        let mut current_point = StaticVec::new_vec([start_point.0, start_point.1]);
        let current_dual_point = current_point.map(DualVec64::<2>::from).derive();
        let mut current = self.target_function(current_dual_point[0], current_dual_point[1]);

        let mut stepsize_index = 0;
        let max_stepsize_index = self.steps_sizes().len();
        let max_iter = 1000;
        let mut iter = 0;
        loop {
            // Ein Iterationsschritt
            iter = iter + 1;
            if iter > max_iter {
                break;
            }
            let new_point = current_point
                - current.eps * self.steps_sizes()[stepsize_index] / current.eps.norm();
            let new_dual_point = new_point.map(DualVec64::<2>::from).derive();
            let new = self.target_function(new_dual_point[0], new_dual_point[1]);
            hist.push(
                (current_point[0], current_point[1]),
                (new_point[0], new_point[1]),
            );
            if (new.re - current.re).abs() < self.eps() {
                break;
            }
            // Nebenbedingungen
            if self.violates_constraint_functions(new_dual_point[0], new_dual_point[1]) {
                stepsize_index = if stepsize_index + 1 == max_stepsize_index {
                    stepsize_index
                } else {
                    stepsize_index + 1
                }; // Verkleinerung der Schrittweite, sofern möglich.
                continue; // Neuer Iterationsschritt
            }
            if new.re < current.re {
                // Zielfunktion wurde verbessert
                current_point = new_point;
                current = new;
                stepsize_index = if stepsize_index == 0 {
                    stepsize_index
                } else {
                    stepsize_index - 1
                }; // Vergrößerung der Schrittweite, sofern möglich.
                continue; // Neuer Iterationsschritt
            } else {
                stepsize_index = if stepsize_index + 1 == max_stepsize_index {
                    stepsize_index
                } else {
                    stepsize_index + 1
                }; // Verkleinerung der Schrittweite, sofern möglich.
                continue; // Neuer Iterationsschritt
            }
        }

        //println!("{}",current_point);
        return (current_point[0], current_point[1]);
    }
}

struct BelegFunction {
    pub start_point: (f64, f64),
    pub search_steps: Vec<f64>,
    pub eps: f64,
}

impl GradientOptimization2D for BelegFunction {
    fn target_function<D: DualNum<f64>>(&self, x1: D, x2: D) -> D {
        (x1 * x2).mul(100.0) + x1 / x2 + (x1 * x2).powi(-1)
    }
    fn violates_constraint_functions<D: DualNum<f64>>(&self, x1: D, x2: D) -> bool {
        !(x1.re() >= 0.25 && x1.re() <= 1.0 && x2.re() >= 0.25 && x2.re() <= 1.0)
    }
    fn start_point(&self) -> (f64, f64) {
        self.start_point
    }
    fn steps_sizes(&self) -> &Vec<f64> {
        &self.search_steps
    }
    fn eps(&self) -> f64 {
        self.eps
    }
}

fn visualize(history: History) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new("sol.png", (1080, 1080)).into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.margin(10, 10, 10, 10);
    // After this point, we should be able to construct a chart context
    let mut chart = ChartBuilder::on(&root)
        // Set the caption of the chart
        .caption("Iterationsverlauf", ("sans-serif", 40).into_font())
        // Set the size of the label region
        .x_label_area_size(20)
        .y_label_area_size(40)
        // Finally attach a coordinate on the drawing area and make a chart context
        .build_cartesian_2d(0f64..1.25f64, 0f64..1.25f64)?;

    // Then we can draw a mesh
    chart
        .configure_mesh()
        // We can customize the maximum number of labels allowed for each axis
        .x_labels(5)
        .y_labels(5)
        // We can also change the format of the label text
        .y_label_formatter(&|x| format!("{:.3}", x))
        .draw()?;

    chart.draw_series(LineSeries::new(history.current, &RED))?;
    chart.draw_series(LineSeries::new([(0.25,0.25),(0.25,1.0),(1.0,1.0),(1.0,0.25),(0.25,0.25)], &BLACK))?;
    // Similarly, we can draw point series
    chart.draw_series(PointSeries::of_element(
        history.next,
        2,
        &RED,
        &|c, s, st| {
            return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
            + Circle::new((0,0),s,st.filled()) // At this point, the new pixel coordinate is established
            ;
        },
    ))?;

    root.present()?;
    Ok(())
}

fn main() {
    let opt = BelegFunction {
        start_point: (0.5, 0.5),
        search_steps: vec![0.01, 0.005, 0.0025, 0.001],
        eps: 0.01,
    };
    let mut hist = History {
        current: Vec::new(),
        next: Vec::new(),
    };
    let res = opt.gradient_minimize(&mut hist);
    println!("{},{}", res.0, res.1);
    visualize(hist).unwrap();
}
