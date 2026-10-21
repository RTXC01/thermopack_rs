use std::{hint::black_box, time::Instant};
use thermopack_rs::spinodal_properties::*;
use thermopack_rs::{Cubic, CubicKind, Phase};
fn main() {
    let check = std::env::args().any(|a| a == "--check");
    for (profile, (kind, name)) in [
        (CubicKind::PengRobinson, "C1"),
        (CubicKind::SoaveRedlichKwong, "N2"),
    ]
    .into_iter()
    .enumerate()
    {
        let model = Cubic::new(kind, &[name]).unwrap();
        let tc = model.components()[0].critical_temperature;
        let start = Instant::now();
        let prepared = model
            .prepare_pure_spinodal_search(&[1.0], 1e6, None)
            .unwrap();
        let setup = start.elapsed().as_nanos() as f64;
        for mode in ["fresh", "prepared"] {
            let solve = |i: usize| {
                let target = black_box(tc * (0.94 + 0.04 * (i % 25) as f64 / 24.0));
                let point = if mode == "fresh" {
                    baseline_pure_spinodal_intersection(
                        black_box(&model),
                        &[1.0],
                        SpinodalProperty::Temperature,
                        target,
                        1e6,
                        Some(Phase::Liquid),
                        None,
                    )
                } else {
                    black_box(&prepared).locate(
                        SpinodalProperty::Temperature,
                        target,
                        Some(Phase::Liquid),
                    )
                }
                .unwrap()
                .unwrap();
                [point.temperature, point.volume, point.pressure]
            };
            if check {
                for i in 0..25 {
                    let p = solve(i);
                    println!(
                        "check_{}_{mode}_{i},{:.17e},{:.17e},{:.17e}",
                        profile + 1,
                        p[0],
                        p[1],
                        p[2]
                    );
                }
                continue;
            }
            let iterations = if cfg!(debug_assertions) { 1 } else { 100 };
            let mut sums = [0.0; 3];
            let start = Instant::now();
            for i in 0..iterations {
                let values = black_box(solve(i));
                for j in 0..3 {
                    sums[j] += values[j];
                }
            }
            println!(
                "native_{}_{mode},{:.6},{:.17e},{:.17e},{:.17e}",
                profile + 1,
                start.elapsed().as_nanos() as f64 / iterations as f64,
                sums[0],
                sums[1],
                sums[2]
            );
        }
        if !check {
            println!("setup_{},{setup:.6}", profile + 1);
        }
    }
}
