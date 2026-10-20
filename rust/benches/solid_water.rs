use std::{hint::black_box, time::Instant};
use thermopack_rs::{SolidWater, SolidWaterState};
fn fields(s: SolidWaterState) -> [f64; 12] {
    let g = s.gibbs;
    [
        g.value,
        g.dp,
        g.dpp,
        g.dt,
        g.dtt,
        g.dtp,
        s.volume,
        s.entropy,
        s.enthalpy,
        s.internal_energy,
        s.helmholtz_energy,
        s.baseline_heat_capacity_p,
    ]
}
fn main() {
    let mut model = SolidWater::new();
    model.initialize(100., -10000.).unwrap();
    let prepared = model.prepare(250.).unwrap();
    let check = std::env::args().any(|s| s == "--check");
    let n = if cfg!(debug_assertions) { 100 } else { 500000 };
    for mode in ["fresh", "prepared"] {
        let evaluate = |p| {
            if mode == "fresh" {
                black_box(&model).state(black_box(250.), p).unwrap()
            } else {
                black_box(&prepared).state(p).unwrap()
            }
        };
        if check {
            for i in 0..100 {
                let values = fields(evaluate(1e5 + i as f64 * 1e6));
                println!(
                    "check_{mode}_{i},{}",
                    values
                        .iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                );
            }
            continue;
        }
        let start = Instant::now();
        let mut sum = 0.;
        for i in 0..n {
            let values = black_box(fields(evaluate(black_box(1e5 + (i % 100) as f64 * 1e6))));
            sum += values.iter().sum::<f64>();
        }
        println!(
            "{mode},{:.6},{sum:.17e}",
            start.elapsed().as_nanos() as f64 / n as f64
        );
    }
}
