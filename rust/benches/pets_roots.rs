use std::{hint::black_box, time::Instant};
use thermopack_rs::{Pets, Phase};

fn main() {
    let model = Pets::new().unwrap();
    let check = std::env::args().any(|s| s == "--check");
    let iterations = if cfg!(debug_assertions) { 100 } else { 20000 };
    for (name, t, phase) in [
        ("liquid", 120., Phase::Liquid),
        ("vapor", 120., Phase::Vapor),
        ("supercritical", 300., Phase::Vapor),
    ] {
        let prepared = model.prepare_equilibrium(t).unwrap();
        for mode in ["fresh", "prepared", "baseline"] {
            let evaluate = |p| match mode {
                "fresh" => black_box(&model).volume(black_box(t), p, phase).unwrap(),
                "prepared" => black_box(&prepared).volume(p, phase).unwrap(),
                _ => {
                    black_box(&model)
                        .baseline_volume_iteration(black_box(t), p, phase)
                        .unwrap()
                        .volume
                }
            };
            if check {
                for i in 0..100 {
                    let p = 5e5 + i as f64 * 5e3;
                    let v = evaluate(p);
                    let state = model.pressure_tv(t, v, 1.).unwrap();
                    println!(
                        "check_{name}_{mode}_{i},{v:.17e},{:.17e},{:.17e},{:.17e}",
                        state.value, state.dt, state.dv
                    );
                }
                continue;
            }
            let mut sum = 0.;
            let start = Instant::now();
            for i in 0..iterations {
                sum += evaluate(black_box(5e5 + (i % 100) as f64 * 5e3));
            }
            println!(
                "{name}_{mode},{:.6},{sum:.17e}",
                start.elapsed().as_nanos() as f64 / iterations as f64
            );
        }
    }
}
