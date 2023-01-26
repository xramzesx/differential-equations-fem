use crate::utils::solver::gaussian_elimination;
use crate::utils::calc::integrate;

//// MAIN FUNCTIONS ////

pub fn generate_solution (
    n: i32, 
    from: f64, 
    to: f64
) -> impl Fn(f64) -> f64 {

    //// GENERATE BASIS FUNCTIONS ////

    let mut e = Vec::new();
    let mut de = Vec::new();

    let h = (to - from) / (n as f64);


    for i in 0..n {
        e.push(get_e_i(
            from + i as f64 * h, 
            h
        ));
    }

    for i in 0..n {
        de.push(get_de_i(
            from + i as f64 * h, 
            h
        ));
    }

    //// SOLVE EQUATIONS ////
    
    println!("Generating equation matrix");
    let mut matrix = generate_matrix(n, &e, &de, from, to);
    
    println!("Solving equations");
    let result = gaussian_elimination(&mut matrix);
    
    move |x| { result.iter().enumerate().fold(0.0, |sum, (i, u) | sum + u * e[i](x) )}
}

pub fn generate_matrix(
    n: i32, 
    e: &Vec<impl Fn(f64) -> f64>,
    de: &Vec<impl Fn(f64) -> f64>,
    from: f64, 
    to: f64
) -> Vec<Vec<f64>> {

    //// Initialize n x n + 1 matrix 
    let mut result = vec![vec![0.0; (n + 1) as usize]; n as usize];

    //// FILL MATRIX ////
    
    for i in 0..(n as usize) {
        result[i][i] = B(&n, &from, &to, &e, &de, i, i);
    }

    for i in 1..(n as usize) {
        result[i - 1][i] = B(&n, &from, &to, &e, &de, i - 1, i);
        result[i][i - 1] = result[i - 1][i];
    }

    for i in 0..(n as usize) {
        result[i][n as usize] = L(&n, &from, &to, &e, &de, i);
    }

    result
}
pub fn get_e_i ( 
    xi : f64, 
    h : f64 
) -> impl Fn(f64) -> f64 {
    return move |x| {
        0.0_f64.max(1.0 - ( (x - xi) / h ).abs())
    };
}

pub fn get_de_i(
    xi : f64,
    h : f64
) -> impl Fn(f64) -> f64 {
    return move |x| {
        if xi - h < x && x <= xi {
            1f64 / h 
        } else if xi < x && x < xi + h {
            -1f64 / h 
        } else {
            0.0
        }
    }
}

//// INTEGRAL FUNCTIONS ////
pub const POINTS: i32 = 100;

fn B ( 
    n : &i32, 
    from: &f64,
    to : &f64,
    e: &Vec<impl Fn(f64) -> f64>,
    de: &Vec<impl Fn(f64) -> f64>,
    u: usize,
    v: usize
) -> f64 {
    let h = (to - from) / (*n as f64);

    let left_boundary = from
        .max(from + (u as f64 - 1f64) * h)
        .max(from + (v as f64 - 1f64) * h);
    let right_boundary = to
        .min(from + (u as f64 + 1f64) * h)
        .min(from + (v as f64 + 1f64) * h);

    - e[v](0.0) * e[u](0.0) + integrate( 
        move |x| de[v](x) * de[u](x), 
        0f64.max(left_boundary), 
        2f64.min(right_boundary), 
        POINTS
    )
}

fn L (
    n : &i32, 
    from: &f64,
    to : &f64,
    e: &Vec<impl Fn(f64) -> f64>,
    de: &Vec<impl Fn(f64) -> f64>,
    v : usize
) -> f64 {

    let h = (to - from) / (*n as f64);

    let left_boundary = from
        .max(from + (v as f64 - 1f64) * h);
    let right_boundary = to
        .min(from + (v as f64 + 1f64) * h);
    
    let first_integral = integrate( 
        |x:f64| 100.0 * x / (x + 1.0) * get_e_i(from + v as f64 * h, h)(x), 
        0f64.max(left_boundary).min(right_boundary), 
        1f64.min(right_boundary).max(left_boundary), 
        POINTS 
    ); 

    let second_integral = 50.0 * integrate(
        move |x| e[v](x),
        1f64.max(left_boundary),
        2f64.min(right_boundary), 
        POINTS
    );

    let free_term = - 20.0 * e[v](0.0);

    first_integral + second_integral + free_term
}
//// TESTS ///

#[cfg(test)]
pub mod tests {
    use crate::utils::fem::generate_matrix;
    use crate::utils::fem::generate_solution;
    use crate::get_e_i;
    use crate::integrate;
    
    #[test] 
    fn test_integrate() {
        let n:i32 = 5;
        let h:f64 = 2f64 / (n as f64);

        let f = |x:f64| 100.0 * x / (x + 1.0) * get_e_i(0.0, h)(x);

        println!("F(.2): {}", f(0.2));

        println!("IntF : {}", integrate(f, 0f64, 1f64, 1000))
    }
}