use std::f64::consts::PI;
use std::hash::{Hash, Hasher};

//// MAIN FUNCTIONS ////

pub fn derivative ( 
    f : impl Fn(f64) -> f64  
) -> impl Fn( f64 ) -> f64 {  
    const H: f64 = 0.0000000001;    
    return move |x| (f(x + H) - f(x - H)) / ( 2.0 * H )
}

pub fn integrate ( 
    f : impl Fn(f64) -> f64, 
    from: f64, 
    to: f64,
    points : i32
) -> f64 {
    let Roots{x, w} = get_roots(points, -1.0, 1.0);
    let p = (to - from) / 2f64;
    let q = (to + from) / 2f64;

    let mut sum = 0.0;
    
    for i in 0..x.len() {
        sum += w[i as usize] * f(p * x[i as usize] + q);
    }
    return sum * p;
} 

//// UTILITIES ////
static mut ROOT_NS : Vec<i32> = Vec::new();
static mut ROOT_XS : Vec<Vec<f64>> = Vec::new();
static mut ROOT_WS : Vec<Vec<f64>> = Vec::new();

unsafe fn get_cached_roots (n1: &i32) -> (Vec<f64>, Vec<f64>) {
    let mut x1 = vec![0.0; *n1 as usize];
    let mut w1 = vec![0.0; *n1 as usize];

    for i in 0..ROOT_NS.len() {
        let n = ROOT_NS[i as usize].clone();
        let x = ROOT_XS[i as usize].clone();
        let w = ROOT_XS[i as usize].clone();
    
        if n == *n1 {
            return (x, w)
        }
    }
    
    (x1, w1)
}

unsafe fn insert_cached_roots (n1: &i32, x1: &Vec<f64>, w1 :&Vec<f64> ) {
    let n = n1.clone();
    let x = x1.clone();
    let w = w1.clone();

    ROOT_NS.push(n);
    ROOT_XS.push(x);
    ROOT_WS.push(w);
}

#[derive(Debug, Clone)]
pub struct Roots {
    x: Vec<f64>, 
    w: Vec<f64>,
}

impl Hash for Roots {
    fn hash<H: Hasher>(&self, state: &mut H) {
        for item in &self.x {
            item.to_bits().hash(state);
        }
        for item in &self.w {
            item.to_bits().hash(state);
        }
    }
}

pub fn get_roots ( n : i32, x1: f64, x2: f64 ) -> Roots {    
    //// TODO: make it works ////////////////////////////////////////
    // This part of the code is unsafe and was created for caching 
    // repeated results. Unfortunately, there were some bugs 
    // related to storing smth as a global, mutable variable in Rust.
    //
    // Although this part of the code sometimes works, it is 
    // kept for potential future contributors who may want 
    // to reuse this project.
    ////////////////////////////////////////////////////////////////
    
    // let (mut x, mut w) : (Vec<f64>,Vec<f64>) = unsafe {
    //     get_cached_roots(&n)
    // };
    
    // let root_contains = unsafe { ROOT_NS.contains(&n) };
    // if root_contains {
    //     return Roots{ x: x, w : w};
    // }

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
    // TODO: make it works
    // unsafe {
    //     insert_cached_roots (&n, &x, &w);
    // }
    return Roots{x : x, w : w};
}





//// TESTS ////

#[cfg(test)]
pub mod tests {
    use crate::utils::calc::get_roots;
    use crate::utils::calc::Roots;

    #[test]
    fn test_get_roots() {
        let ex: Vec<f64> = vec![0.42264973081037427, 1.5773502691896257 ];
        let ew: Vec<f64> = vec![1.0000000000000002, 1.0000000000000002 ];        
        
        let Roots{x, w} = get_roots(2, 0.0, 2.0);

        for i in 0..x.len() {
            assert_eq!(ex[i as usize], x[i as usize]);
            assert_eq!(ew[i as usize], w[i as usize]);
        }
    }

}

