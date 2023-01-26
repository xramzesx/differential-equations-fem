mod utils;
use utils::fem::*;
use plotters::prelude::*;

fn main () -> Result<(), Box<dyn std::error::Error>>  {
    
    if std::env::args().len() < 2 {
        panic!("Please enter at least one parameter (n)");
    }
    
    let n = std::env::args()
        .nth(1)
        .and_then(|n| n.parse::<i32>().ok())
        .unwrap_or(-1);

    let file = std::env::args()
        .nth(2)
        .unwrap_or("output.png".to_owned());

    if n < 2 {
        panic!("n must be greater than 2!");
    }

    let from :f64 = 0.0;
    let to   :f64 = 2.0;

    let u = generate_solution(n, from, to);

    println!("Generating plot");

    let mut max_v : f32 = 2.0;
    let mut min_v : f32 = -2.0;

    for x in -200..=200 {
        max_v = max_v.max( u( x as f64 / 100.0) as f32 );
        min_v = min_v.min( u( x as f64 / 100.0) as f32 );
    }


    let root = BitMapBackend::new(&file, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Równianie transportu ciepła", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0f32..2f32, min_v..(max_v + 1.0))?;

    chart.configure_mesh().draw()?;

    chart
        .draw_series(LineSeries::new(
            (-100..=100).map(|x| x as f32 / 50.0).map(|x| (x, u(x as f64) as f32)),
            &RED,
        ))?
        .label("u(x)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;
    print!("Done!");
    Ok(())
}
