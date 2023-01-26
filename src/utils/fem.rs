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

    // println!("> {:?}", matrix);

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
        // result[i][i] = B(&e[i],&de[i], &e[i], &de[i]);
        result[i][i] = B(&n, &from, &to, &e, &de, i, i);
    }


    for i in 1..(n as usize) {
        result[i - 1][i] = B(&n, &from, &to, &e, &de, i - 1, i);
        // result[i - 1][i] = B(&e[i-1], &de[i-1], &e[i], &de[i]);
        result[i][i - 1] = result[i - 1][i];
    }

    // for i in 0..(n as usize) {
    //     for j in 0..(n as usize) {
    //         result[i][j] = B(&n, &from, &to, &e, &de, i, j);

    //     }
    // }

    for i in 0..(n as usize) {
        // result[i][n as usize] = L(&e[i], &de[i]);
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
pub const POINTS: i32 = 1000;

fn B ( 
    n : &i32, 
    from: &f64,
    to : &f64,
    e: &Vec<impl Fn(f64) -> f64>,
    de: &Vec<impl Fn(f64) -> f64>,
    u: usize,
    v: usize
) -> f64 {
    // - v(0.0) * u(0.0) + integrate( move |x| derivative(v)(x) * derivative(u)(x), 0.0, 2.0, POINTS)
    // - (v.f)(0.0) * (u.f)(0.0) + integrate( move |x| (v.df)(x) * (u.df)(x), 0.0, 2.0, POINTS)

    // let leftBoundary = 

    let h = (to - from) / (*n as f64);

    let leftBoundary = from
        .max(from + (u as f64 - 1f64) * h)
        .max(from + (v as f64 - 1f64) * h);
    let rightBoundary = to
        .min(from + (u as f64 + 1f64) * h)
        .min(from + (v as f64 + 1f64) * h);

    - e[v](0.0) * e[u](0.0) + integrate( 
        move |x| de[v](x) * de[u](x), 
        0f64.max(leftBoundary), 
        2f64.min(rightBoundary), 
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

    let leftBoundary = from
        .max(from + (v as f64 - 1f64) * h);
    let rightBoundary = to
        .min(from + (v as f64 + 1f64) * h);

    // println!("lr: {} | {} : {}", leftBoundary, rightBoundary, e[v](0.2));
    
    // let t = move |x| 100.0 * x / (x + 1.0) * (*e[v])(x);
    
    // let f = |x:f64| 100.0 * x / (x + 1.0) * get_e_i(from + h * (v as f64), h)(x);
    let f = |x:f64| 100.0 * x / (x + 1.0) * get_e_i(0.0, h)(x);
    
    // println!("h: {}", h);
    // println!("IntF : {}", integrate(f, 0f64, 1f64, POINTS));
    // println!("t: {}", f(0.2));

    let firstIntegral = integrate( 
        |x:f64| 100.0 * x / (x + 1.0) * get_e_i(from + v as f64 * h, h)(x), 
        // move |x| x / (x + 1.0) * e[v](x), 
        0f64.max(leftBoundary).min(rightBoundary), 
        1f64.min(rightBoundary).max(leftBoundary), 
        POINTS 
    ); 

    let secondIntegral = 50.0 * integrate(
        move |x| e[v](x),
        1f64.max(leftBoundary),
        2f64.min(rightBoundary), 
        POINTS
    );

    let freeTerm = - 20.0 * e[v](0.0);

    // println!("[{}][{}]({},{}): {} {} {}", v,from + h * (v as f64),leftBoundary, rightBoundary, firstIntegral, secondIntegral, freeTerm);

    firstIntegral + secondIntegral + freeTerm
}
//// TESTS ///

#[cfg(test)]
pub mod tests {
    use crate::utils::fem::generate_matrix;
    use crate::utils::fem::generate_solution;
    use crate::get_e_i;
    use crate::integrate;
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

    
    #[test] 
    fn test_integrate() {
        let n:i32 = 5;
        let h:f64 = 2f64 / (n as f64);

        let f = |x:f64| 100.0 * x / (x + 1.0) * get_e_i(0.0, h)(x);

        println!("F(.2): {}", f(0.2));

        println!("IntF : {}", integrate(f, 0f64, 1f64, 1000))
    }
}