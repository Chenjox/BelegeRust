use ndarray::Array2;


struct GebrochenRational {
    p_vector: Vec<f64>,
    q_vector: Vec<f64>,
    n: usize
}

impl GebrochenRational {
    fn new(p_vector: Vec<f64>, q_vector: Vec<f64>) -> Self {
        if p_vector.len() -1 != q_vector.len() {
            panic!("Lenghts of given Vectors do not match! {} != {}", p_vector.len(), q_vector.len());
        }
        let mut q_vector = q_vector;
        let n = p_vector.len();
        Self { p_vector , q_vector, n }
    }

    fn get_linear_approximation(&self) -> (Array2<f64>,Array2<f64>) {
        let mut p = self.p_vector.clone();
        let mut q = self.q_vector.clone();
        q.push(0.0);
        let n = self.n;
        let mut A : Array2<f64> = Array2::zeros([n-1,n-1]);
        let mut B : Array2<f64> = Array2::zeros([n-1,n-1]);

        println!("{},{}",p.len(),q.len());
        // FIXME: Funktioniert noch nicht.
        for i in 0..(n-1) {
            let mut r = vec![0.0; self.p_vector.len()];
            // Indexschlacht die erste (s_0 und s_1)
            // Zerlegen in Linearen 
            let s_1 = p[n-i-1]/q[n-2-i];
            let s_0 = if i != n-2 {(q[n-3-i]-s_1*p[n-i-2])/(q[n-2-i])} else { 0.0 };

            // Indexschlacht die zweite
            if i < n-2 {
                for g in 1..(n-i-2) {
                    let g = n-i-2 -g;
                    r[g] = p[g] - s_0 * q[g] - s_1 * q[g-1];
                }
                r[0] = p[0] - s_0 * q[0];
            }

            
            println!("{:?}", r);

            p = q;
            q = r;


            A[[i,i]] = s_0;
            B[[i,i]] = s_1;

            // Indexschlacht die zweite
            
        }
        return (A,B);
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
            1.0,
            0.00490648418,
            2.71147593e-5,
            5.42912842e-8,
            6.94958499e-12,
        ]);
    g.get_linear_approximation();
    println!("Hello World");
}