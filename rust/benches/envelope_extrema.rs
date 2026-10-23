use std::{hint::black_box, time::Instant};
use thermopack_rs::envelope::EnvelopeExtremum;
use thermopack_rs::{
    Cubic, CubicKind, CubicOptions, MixingRule, MixtureSaturation, SaturationKind,
};

fn main() {
    let check = std::env::args().any(|s| s == "--check");
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
            if eos == "SRK" {
                CubicKind::SoaveRedlichKwong
            } else {
                CubicKind::PengRobinson
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
        for row in include_str!("../../tests/fixtures/envelope_extrema.csv")
            .lines()
            .skip(1)
        {
            let case = row.split(';').next().unwrap();
            let f: Vec<_> = case.split('|').collect();
            if f[1] != eos
                || f[2] != names
                || f[3] != config
                || f[8] != "100000.0,15000000.0,80.0,0.0,1.0,0.25,1000"
            {
                continue;
            }
            let g: Vec<f64> = f[7].split(',').map(|v| v.parse().unwrap()).collect();
            let beta: f64 = f[6].parse().unwrap();
            let mut x: Vec<_> = g[2..]
                .iter()
                .map(|k| 0.5 / (1. - beta + beta * k))
                .collect();
            let mut y: Vec<_> = x.iter().zip(&g[2..]).map(|(x, k)| x * k).collect();
            let sx = x.iter().sum::<f64>();
            let sy = y.iter().sum::<f64>();
            x.iter_mut().for_each(|v| *v /= sx);
            y.iter_mut().for_each(|v| *v /= sy);
            inputs.push((
                if f[5] == "P" {
                    EnvelopeExtremum::Pressure
                } else {
                    EnvelopeExtremum::Temperature
                },
                MixtureSaturation {
                    kind: if beta == 0. {
                        SaturationKind::Bubble
                    } else {
                        SaturationKind::Dew
                    },
                    temperature: g[0],
                    pressure: g[1],
                    liquid: x,
                    vapor: y,
                    liquid_volume: 1.,
                    vapor_volume: 1.,
                    dpdt: 0.,
                    iterations: 0,
                    fugacity_residual: 0.,
                },
            ));
        }
        assert_eq!(inputs.len(), 2);
        let solve = |i: usize| {
            let (kind, guess) = black_box(&inputs[i % inputs.len()]);
            let result = black_box(&model)
                .refine_envelope_extremum(black_box(&[0.5, 0.5]), guess, *kind)
                .unwrap();
            let p = result.state;
            let mut v = vec![
                p.temperature,
                p.pressure,
                f64::from(p.kind == SaturationKind::Dew),
                p.liquid_volume,
                p.vapor_volume,
            ];
            v.extend(p.liquid);
            v.extend(p.vapor);
            v.push(result.stationarity_residual);
            v
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
            let mut sums = [0.; 10];
            let start = Instant::now();
            for i in 0..iterations {
                for (s, v) in sums.iter_mut().zip(black_box(solve(i))) {
                    *s += v;
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
