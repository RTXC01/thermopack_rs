use std::{hint::black_box, time::Instant};
use thermopack_rs::envelope::{EnvelopeEnd, EnvelopeOptions};
use thermopack_rs::{Cubic, CubicKind, CubicOptions, MixingRule, SaturationKind};

fn main() {
    let check = std::env::args().any(|s| s == "--check");
    for profile in 1..=4 {
        let names: &[&str] = match profile {
            3 => &["C1", "C2", "C3"],
            4 => &["C1", "C3"],
            _ => &["C1", "C2"],
        };
        let z = if profile == 3 {
            vec![0.2, 0.3, 0.5]
        } else {
            vec![0.5, 0.5]
        };
        let n = z.len();
        let model = Cubic::with_options(
            if profile == 2 {
                CubicKind::SoaveRedlichKwong
            } else {
                CubicKind::PengRobinson
            },
            names,
            CubicOptions {
                mixing_rule: if profile == 4 {
                    MixingRule::Unifac
                } else {
                    MixingRule::VdW
                },
                ..Default::default()
            },
        )
        .unwrap();
        let solve = || {
            let curve = black_box(&model)
                .baseline_phase_envelope(black_box(&z), EnvelopeOptions::default())
                .unwrap();
            assert_eq!(curve.end, EnvelopeEnd::InitialPressure);
            curve
        };
        if check {
            for (i, p) in solve().points.iter().enumerate() {
                let values: Vec<_> = [
                    p.temperature,
                    p.pressure,
                    f64::from(p.kind == SaturationKind::Dew),
                    p.liquid_volume,
                    p.vapor_volume,
                ]
                .into_iter()
                .chain(p.liquid.iter().copied())
                .chain(p.vapor.iter().copied())
                .collect();
                println!(
                    "check_{profile}_{i},{}",
                    values
                        .iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                );
            }
            continue;
        }
        let iterations = if cfg!(debug_assertions) { 1 } else { 10 };
        let mut sums = vec![0.; 6 + 2 * n];
        let start = Instant::now();
        for _ in 0..iterations {
            let curve = black_box(solve());
            sums[0] += curve.points.len() as f64;
            for p in curve.points {
                sums[1] += p.temperature;
                sums[2] += p.pressure;
                sums[3] += p.liquid_volume;
                sums[4] += p.vapor_volume;
                sums[5] += f64::from(p.kind == SaturationKind::Dew);
                for i in 0..n {
                    sums[6 + i] += p.liquid[i];
                    sums[6 + n + i] += p.vapor[i];
                }
            }
        }
        println!(
            "native_{profile},{:.6},{}",
            start.elapsed().as_nanos() as f64 / iterations as f64,
            sums.iter()
                .map(|v| format!("{v:.17e}"))
                .collect::<Vec<_>>()
                .join(",")
        );
    }
}
