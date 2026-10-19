use std::{hint::black_box, time::Instant};
use thermopack_rs::{Mbwr, MbwrKind, Phase};
fn main() {
    let check = std::env::args().any(|s| s == "--check");
    let iterations = if cfg!(debug_assertions) { 100 } else { 20000 };
    for (family, kind) in [("19", MbwrKind::Bender19), ("32", MbwrKind::Mbwr32)] {
        let model = Mbwr::new(kind, "C1").unwrap();
        let parameters = model.source_density_parameters().unwrap();
        for (workload, t, phase) in [
            ("liquid", 150., Phase::Liquid),
            ("vapor", 150., Phase::Vapor),
            ("supercritical", 300., Phase::Vapor),
        ] {
            let prepared = model.prepare_baseline_density(t, parameters).unwrap();
            for mode in ["fresh", "prepared"] {
                let name = format!("{family}_{workload}_{mode}");
                let evaluate = |p| {
                    if mode == "fresh" {
                        black_box(&model)
                            .baseline_density(black_box(t), p, phase, black_box(parameters))
                            .unwrap()
                    } else {
                        black_box(&prepared).density(p, phase).unwrap()
                    }
                };
                if check {
                    for i in 0..100 {
                        let r = evaluate(5e5 + i as f64 * 5e3);
                        println!(
                            "check_{name}_{i},{:.17e},{},{:.17e},{:.17e}",
                            r.density,
                            if r.phase == Phase::Liquid { 1 } else { 2 },
                            r.pressure,
                            r.pressure_density_derivative
                        );
                    }
                    continue;
                }
                let mut sum = 0.;
                let start = Instant::now();
                for i in 0..iterations {
                    let r = black_box(evaluate(black_box(5e5 + (i % 100) as f64 * 5e3)));
                    sum += r.density
                        + (if r.phase == Phase::Liquid { 1. } else { 2. })
                        + r.pressure * 1e-9
                        + r.pressure_density_derivative * 1e-5;
                }
                println!(
                    "{name},{:.6},{sum:.17e}",
                    start.elapsed().as_nanos() as f64 / iterations as f64
                );
            }
        }
    }
}
