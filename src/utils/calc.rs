use std::f64::consts::PI;

//// MAIN FUNCTIONS ////

pub fn derivative ( 
    f : fn(f64) -> f64  
) -> impl Fn( f64 ) -> f64 {
    const H: f64 = 0.00001;    
    return move |x| (f(x + H) - f(x - H)) / ( 2.0 * H )
}

pub fn integrate ( 
    f : fn(f64) -> f64, 
    from: f64, 
    to: f64,
    points : i32
) -> f64 {
    let (x, w) = get_roots(points, from, to);
    let mut sum = 0.0;
    for i in 0..x.len() {
        sum += w[i as usize] * f(x[i as usize]);
    }
    return sum;
} 

//// UTILITIES ////

pub fn get_roots ( n : i32, x1: f64, x2: f64 ) -> (Vec<f64>, Vec<f64>) {
    const EPS: f64 = 1.0e-14;
    let mut x: Vec<f64> = vec![0.0; n as usize];
    let mut w: Vec<f64> = vec![0.0; n as usize];
    
    let xm: f64 = 0.5 * (x2 + x1); 
    let xl: f64 = 0.5 * (x2 - x1); 
    let mut pp: f64 = 1.0;

    let mid: i32 = (n + 1) / 2;
    
    for i in 0..mid {
        let mut z: f64 = (PI * (f64::from(i) + 0.75) / (n as f64 + 0.5)).cos();
        let mut z1:f64 = z + EPS * 2.0;
        
        while (z -z1).abs() > EPS {
            let mut p1: f64 = 1.0;
            let mut p2: f64 = 0.0;
            let mut tmp:f64;
            
            for j in 0..n {
                tmp = p2;
                p2 = p1;
                p1 = ((2.0 * f64::from(j) + 1.0) * z * p2 - j as f64 * tmp) / (f64::from(j) + 1.0);
            }

            pp = f64::from(n) * ( z * p1 - p2 ) / ( z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        }

        // Scale the root to the desired interval,
        x[i as usize]=xm-xl*z;
        
        // and put in its symmetric counterpart.
        x[(n - 1 - i) as usize] = xm + xl * z; 
        
        // Compute the weight
        w[i as usize]= 2.0 * xl / ( ( 1.0 -z * z) * pp * pp); 
        
        // and its symmetric counterpart.
        w[(n-1-i) as usize]=w[i as usize];
  
    }

    return (x, w);
}

//// TESTS ////

#[cfg(test)]
pub mod tests {
    use get_roots;

    #[test]
    fn test_get_roots() {
        let ex: Vec<f64> = vec![0.42264973081037427, 1.5773502691896257 ];
        let ew: Vec<f64> = vec![1.0000000000000002, 1.0000000000000002 ];        
        
        let (x, w) = get_roots(2, 0.0, 2.0);

        for i in 0..x.len() {
            assert_eq!(ex[i as usize], x[i as usize]);
        }
    }
}
