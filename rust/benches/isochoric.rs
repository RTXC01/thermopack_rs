use std::{hint::black_box, time::Instant};
use thermopack_rs::{Cubic, CubicKind, CubicOptions, FlashPhase, FlashResult};

fn canonical(f: &FlashResult) -> Vec<f64> {
    let two = f.phase == FlashPhase::TwoPhase;
    let mut values = vec![
        f.temperature,
        f.pressure,
        f64::from(two),
        if two { f.vapor_fraction } else { 0. },
        if two { f.liquid_fraction } else { 1. },
    ];
    values.extend(if two {
        f.liquid.clone()
    } else {
        f.feed.clone()
    });
    values.extend(if two {
        f.vapor.clone()
    } else {
        vec![0.; f.feed.len()]
    });
    values
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let tcpr = args.iter().any(|a| a == "--tcpr");
    let check = args.iter().any(|a| a == "--check");
    let profile = args.iter().any(|a| a == "--profile");
    let model = Cubic::with_options(
        CubicKind::PengRobinson,
        &["C1", "C2", "C3"],
        CubicOptions {
            parameter_reference: if tcpr { "tcPR" } else { "DEFAULT" },
            ..Default::default()
        },
    )
    .unwrap();
    let z = [0.5, 0.3, 0.2];
    // The comparison driver supplies the exact same source-derived inputs to
    // both programs. A standalone Cargo invocation can prepare native inputs.
    let inputs: Vec<Vec<f64>> = if let Some(i) = args.iter().position(|a| a == "--inputs") {
        std::fs::read_to_string(&args[i + 1])
            .unwrap()
            .lines()
            .map(|line| line.split(',').map(|v| v.trim().parse().unwrap()).collect())
            .collect()
    } else {
        (0..100)
            .map(|i| {
                let t = 180. + i as f64 * 0.8;
                let f = model.tp_flash(t, 1e6, &z).unwrap();
                let p = model.equilibrium_properties(&f).unwrap();
                vec![t, 1e6, p.internal_energy, p.entropy, p.volume]
            })
            .collect()
    };
    assert_eq!(inputs.len(), 100);
    for entropy in [false, true] {
        let solve = |row: &[f64]| {
            if entropy {
                let mut options = model.volume_flash_options();
                options.temperature_guess = Some(row[0] * 0.9);
                options.pressure_guess = Some(row[1] * 1.2);
                black_box(&model).sv_flash_with_options(
                    black_box(&z),
                    black_box(row[3]),
                    black_box(row[4]),
                    options,
                )
            } else {
                black_box(&model).uv_flash(black_box(&z), black_box(row[2]), black_box(row[4]))
            }
            .unwrap()
        };
        let name = if entropy { "sv" } else { "uv" };
        for (i, row) in inputs.iter().enumerate() {
            let f = solve(row);
            assert!((f.temperature - row[0]).abs() < 1e-5);
            let p = model.equilibrium_properties(&f).unwrap();
            let value = if entropy {
                p.entropy
            } else {
                p.internal_energy
            };
            let scale = thermopack_rs::R * if entropy { 1. } else { f.temperature };
            assert!((value - row[if entropy { 3 } else { 2 }]).abs() < 1e-9 * scale);
            assert!((p.volume - row[4]).abs() < 1e-9 * row[4]);
            if check {
                println!(
                    "check_{name}_{i},{}",
                    canonical(&f)
                        .iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                );
            }
            if profile {
                let start = Instant::now();
                for _ in 0..5 {
                    black_box(solve(row));
                }
                println!(
                    "profile_{name}_{i},{:.3},{}",
                    start.elapsed().as_nanos() as f64 / 5.,
                    f.iterations
                );
            }
        }
        if check || profile {
            continue;
        }
        let iterations = if cfg!(debug_assertions) { 20 } else { 2000 };
        let start = Instant::now();
        let mut checksum = [0.; 3];
        for i in 0..iterations {
            let f = solve(&inputs[i % inputs.len()]);
            checksum[0] += f.temperature;
            checksum[1] += f.pressure;
            checksum[2] += if f.phase == FlashPhase::TwoPhase {
                f.vapor_fraction
            } else {
                0.
            };
        }
        println!(
            "{name},{:.6},{:.17e},{:.17e},{:.17e}",
            start.elapsed().as_nanos() as f64 / iterations as f64,
            checksum[0],
            checksum[1],
            checksum[2]
        );
    }
}
