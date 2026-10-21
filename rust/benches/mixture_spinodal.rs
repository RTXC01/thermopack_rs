use std::{hint::black_box, time::Instant};
use thermopack_rs::mixture_spinodal::*;
use thermopack_rs::{Cubic, CubicKind};
fn main() {
    let check = std::env::args().any(|a| a == "--check");
    for (profile, (kind, names, z)) in [
        (CubicKind::PengRobinson, vec!["C1", "C2"], vec![0.5, 0.5]),
        (
            CubicKind::SoaveRedlichKwong,
            vec!["C1", "C2"],
            vec![0.5, 0.5],
        ),
        (
            CubicKind::PengRobinson,
            vec!["C1", "C2", "C3"],
            vec![0.2, 0.3, 0.5],
        ),
    ]
    .into_iter()
    .enumerate()
    {
        let model = Cubic::new(kind, &names).unwrap();
        let solve = || {
            baseline_trace_spinodal(
                black_box(&model),
                black_box(&z),
                SpinodalTraceOptions {
                    log_volume_step: 0.04,
                    ..Default::default()
                },
            )
            .unwrap()
        };
        if check {
            for (i, p) in solve().points.iter().enumerate() {
                println!(
                    "check_{}_{i},{:.17e},{:.17e},{:.17e}",
                    profile + 1,
                    p.temperature,
                    p.volume,
                    p.pressure
                );
            }
            continue;
        }
        let mut sums = [0.0; 4];
        let iterations = if cfg!(debug_assertions) { 1 } else { 10 };
        let start = Instant::now();
        for _ in 0..iterations {
            let curve = black_box(solve());
            sums[0] += curve.points.len() as f64;
            for p in curve.points {
                sums[1] += p.temperature;
                sums[2] += p.volume;
                sums[3] += p.pressure;
            }
        }
        println!(
            "native_{},{:.6},{:.17e},{:.17e},{:.17e},{:.17e}",
            profile + 1,
            start.elapsed().as_nanos() as f64 / iterations as f64,
            sums[0],
            sums[1],
            sums[2],
            sums[3]
        );
    }
}
