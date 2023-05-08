
struct GebrochenRational {
    p_vector: Vec<f64>,
    q_vector: Vec<f64>
}

impl GebrochenRational {
    fn new(p_vector: Vec<f64>, q_vector: Vec<f64>) -> Self {
        if p_vector.len() -2 != q_vector.len() {
            panic!("Lenghts of given Vectors do not match! {} != {}", p_vector.len(), q_vector.len());
        }
        Self { p_vector , q_vector }
    }
}

fn main() {
    let g = GebrochenRational::new(
        vec![
            10337706.5,
            48744.172,
            455.205728,
            1.13950676,
            0.00483245447,
            3.11054119e-6,
        ], vec![
            0.00490648418,
            2.71147593e-5,
            5.42912842e-8,
            6.94958499e-12,
        ]);
    println!("Hello World");
}