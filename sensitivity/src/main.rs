

fn rosenbrock(input_vec: &Vec<f64>) -> f64 {
    let mut running_sum = 0.0;
    for i in 0..input_vec.len()-1 {
        running_sum += 100.0*(input_vec[i+1] - input_vec[i].powi(2)).powi(2) + (1.0 - input_vec[i]).powi(2);
    }
    return running_sum;
}

fn rosenbrock_gradient(input_vec: &Vec<f64>) -> Vec<f64> {
    let mut result = vec![0.0;input_vec.len()];
    for i in 0..input_vec.len(){
        if i == 0 {
            result[i] = 400.0*input_vec[0].powi(3) - 400.0*input_vec[0]*input_vec[1]+ 2.0*input_vec[0] - 2.0;
        }else if i < input_vec.len()-1 {
            result[i] = 400.0*input_vec[i].powi(3) - 400.0*input_vec[i]*input_vec[i+1]+ 202.0*input_vec[i] - 200.0*input_vec[i-1].powi(2) - 2.0;
        }else {
            result[i] = 200.0*input_vec[i] - 200.0*input_vec[i-1].powi(2);
        }
    }
    return result;
}

fn gaussian_integration() {
    let min = -2.046;
    let max = 2.046;
    let scale = (max - min)*0.5;
    let shift = (max + min)*0.5;

    let n = 5;
    let m = 3;
    let mvec = [-(3.0f64/5.0).sqrt(), 0.0, (3.0f64/5.0).sqrt()];
    let wvec = [5.0/9.0, 8.0/9.0, 5.0/9.0];
    
}

fn main() {
    println!("Hello World");
}