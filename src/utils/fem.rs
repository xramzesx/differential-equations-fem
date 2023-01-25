use crate::utils::solver::gaussian_elimination;
use crate::utils::calc::derivative;
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

    let mut matrix = generate_matrix(n, &e, &de, from, to);
    
    for row in &matrix {
        println!("> {:?}", row);
    }
    println!("");
    let result = gaussian_elimination(&mut matrix);
    
    for row in &matrix {
        println!("> {:?}", row);
    }

    println!("{:?}", result);

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
        result[i][i] = B(&e[i],&de[i], &e[i], &de[i]);
    }


    for i in 1..(n as usize) {
        result[i - 1][i] = B(&e[i-1], &de[i-1], &e[i], &de[i]);
        result[i][i - 1] = result[i - 1][i];
    }

    // for i in 0..(n as usize) {
    //     for j in 0..(n as usize) {
    //         result[i][j] = B(&e[i], &e[j]);
    //     }
    // }

    for i in 0..(n as usize) {
        result[i][n as usize] = L(&e[i], &de[i]);
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

    // return move |x| {
    //     if xi - h < x && x <= xi {
    //         (x - xi) / h + 1.0
    //     } else if xi < x && x < xi + h {
    //         (-x + xi) / h + 1.0
    //     } else {
    //         0.0
    //     }
    // }
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
// TODO: get this from config 
pub const POINTS: i32 = 10000;

fn B ( 
    // u: &Basis,
    // v: &Basis
    u : &dyn Fn(f64) -> f64,  
    du : &dyn Fn(f64) -> f64,  
    v : &dyn Fn(f64) -> f64,
    dv : &dyn Fn(f64) -> f64
) -> f64 {
    // - v(0.0) * u(0.0) + integrate( move |x| derivative(v)(x) * derivative(u)(x), 0.0, 2.0, POINTS)
    // - (v.f)(0.0) * (u.f)(0.0) + integrate( move |x| (v.df)(x) * (u.df)(x), 0.0, 2.0, POINTS)
    - v(0.0) * u(0.0) + integrate( move |x| dv(x) * du(x), 0.0, 2.0, POINTS)
}

fn L (
    v : &dyn Fn(f64) -> f64,
    dv : &dyn Fn(f64) -> f64
    // v : &Basis
) -> f64 {
    100.0 * integrate( 
        move |x| x / (x + 1.0) * v(x), 
        0.0, 
        1.0, 
        POINTS 
    ) + 50.0 * integrate(
        move |x| v(x),
        1.0,
        2.0, 
        POINTS
    ) - 20.0 * v(0.0)
}
//// TESTS ///

#[cfg(test)]
pub mod tests {
    use crate::utils::fem::generate_matrix;
    use crate::utils::fem::generate_solution;
    #[test]
    fn test_get_e_i() {
        
    }

    #[test]
    fn test_generate_matrix () {
        // let matrix = generate_matrix(3, 0.0, 2.0);
        // let matrix = generate_solution(3, 0.0, 2.0);
        
        // for row in matrix {
        //     println!("{:?}", row);
        // }
    }
}