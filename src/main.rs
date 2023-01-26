mod utils;
use utils::fem::*;
use plotters::prelude::*;
use crate::utils::calc::derivative;
use crate::utils::calc::integrate;

fn main () -> Result<(), Box<dyn std::error::Error>>  {
    
    let n = 50;
    let from :f64 = 0.0;
    let to   :f64 = 2.0;

    let u = generate_solution(n, from, to);
    // let u = derivative( |x| x * x );
    let mut e = Vec::new();
    let h = (to - from) / (n as f64);

    for i in 0..n {
        e.push( get_e_i(from + i as f64 * h, h) );
    }

    let mut max_v : f32 = 2.0;
    let mut min_v : f32 = -2.0;

    for x in -200..=200 {
        max_v = max_v.max( u( x as f64 / 100.0) as f32 );
        min_v = min_v.min( u( x as f64 / 100.0) as f32 );
    }


    println!("{}, {}", max_v, min_v );

    let root = BitMapBackend::new("output.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Równianie transportu ciepła", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0f32..2f32, min_v..(max_v + 1.0))?;

    chart.configure_mesh().draw()?;

    // for i in 0..(e.len() as usize) {
    //     chart
    //         .draw_series(LineSeries::new(
    //             (-100..=100).map(|x| x as f32 / 50.0).map(|x| (x, e[i](x as f64) as f32)),
    //             &RED,
    //         ));
    //     }

    chart
        .draw_series(LineSeries::new(
            (-100..=100).map(|x| x as f32 / 50.0).map(|x| (x, u(x as f64) as f32)),
            &RED,
        ))?
        .label("u(x)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    // chart
    //     .draw_series(LineSeries::new(
    //         (-100..=100).map(|x| x as f32 / 50.0).map(|x| (x, du(x as f64) as f32)),
    //         &RED,
    //     ))?
    //     .label("u(x)")
    //     .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;

    
    
    
    
    for i in 0..20 {
        println!("u: {}", u(i as f64 / 10.0));
    }


    println!("integral : {}",  integrate(|x| 4f64, 0f64, 1f64, 50));
    Ok(())
}
