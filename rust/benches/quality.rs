use std::{hint::black_box, time::Instant};
use thermopack_rs::envelope::EnvelopeConstraint;
use thermopack_rs::quality::{QualityGuess, QualityState};
use thermopack_rs::{Cubic, CubicKind, CubicOptions, MixingRule};
fn fields(p: QualityState) -> Vec<f64> {
    let mut v = vec![
        p.temperature,
        p.pressure,
        p.vapor_fraction,
        p.liquid_volume,
        p.vapor_volume,
    ];
    v.extend(p.liquid);
    v.extend(p.vapor);
    v.push(p.dpdt);
    v
}
fn main() {
    let check = std::env::args().any(|v| v == "--check");
    let iterations = if cfg!(debug_assertions) { 3 } else { 1000 };
    for (profile, (eos, names, config)) in [
        ("PR", "C1,C2", "classic"),
        ("SRK", "C1,C2", "classic"),
        ("PR", "C1,C3", "translated"),
        ("PR", "C1,C3", "unifac"),
    ]
    .into_iter()
    .enumerate()
    {
        let model = Cubic::with_options(
            if eos == "PR" {
                CubicKind::PengRobinson
            } else {
                CubicKind::SoaveRedlichKwong
            },
            &names.split(',').collect::<Vec<_>>(),
            CubicOptions {
                volume_translation: config == "translated",
                mixing_rule: if config == "unifac" {
                    MixingRule::Unifac
                } else {
                    MixingRule::VdW
                },
                ..Default::default()
            },
        )
        .unwrap();
        let mut inputs = Vec::new();
        for row in include_str!("../../tests/fixtures/quality.csv")
            .lines()
            .skip(1)
        {
            let case = row.split(';').next().unwrap();
            let f: Vec<_> = case.split('|').collect();
            if f[1] != eos || f[2] != names || f[3] != config {
                continue;
            }
            let beta = f[5].parse::<f64>().unwrap();
            let o: Vec<_> = f[6].split(':').collect();
            if beta <= 0. || beta >= 1. || o[0] == "K" {
                continue;
            }
            let value = o[1].parse().unwrap();
            let constraint = if o[0] == "T" {
                EnvelopeConstraint::Temperature(value)
            } else {
                EnvelopeConstraint::Pressure(value)
            };
            let g: Vec<f64> = f[7].split(',').map(|v| v.parse().unwrap()).collect();
            inputs.push((
                beta,
                constraint,
                QualityGuess {
                    temperature: g[0],
                    pressure: g[1],
                    ratios: g[2..].to_vec(),
                },
            ));
        }
        assert!(!inputs.is_empty());
        let solve = |i: usize| {
            let (beta, constraint, guess) = black_box(&inputs[i % inputs.len()]);
            fields(
                black_box(&model)
                    .baseline_quality_state_with_guess(
                        black_box(&[0.5, 0.5]),
                        *beta,
                        *constraint,
                        guess,
                    )
                    .unwrap(),
            )
        };
        if check {
            for i in 0..inputs.len() {
                println!(
                    "check_{}_{i},{}",
                    profile + 1,
                    solve(i)
                        .iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                );
            }
        } else {
            let start = Instant::now();
            let mut sums = [0.; 10];
            for i in 0..iterations {
                for (sum, v) in sums.iter_mut().zip(black_box(solve(i))) {
                    *sum += v;
                }
            }
            println!(
                "native_{},{:.6},{}",
                profile + 1,
                start.elapsed().as_nanos() as f64 / iterations as f64,
                sums.iter()
                    .map(|v| format!("{v:.17e}"))
                    .collect::<Vec<_>>()
                    .join(",")
            );
        }
    }
}
