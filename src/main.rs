mod utils;
use utils::calc::*;
// use derivative::*;


fn test (i : f64) -> impl Fn(f64) -> f64 {
    return move |x| { x + i};
}


fn u(x: f64) -> f64 {
    return x*x;
}

fn du(x: f64) -> f64{
    return derivative(u)(x);
}

fn main () {
    let mut x:f64 = 9.0;
    let t = test(1.0)(x);
    println!("{}",t);
    println!("{}",x);
    x+=1.0;

    println!("{}",t);
    println!("{}",x);
    x = 1.0;
    println!("{} , {}", u(x), du(x));

}