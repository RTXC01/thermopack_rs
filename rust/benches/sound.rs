use std::{hint::black_box, time::Instant};
use thermopack_rs::{Cubic, CubicKind, CubicOptions};

fn model(names: &[&str], tcpr: bool) -> Cubic {
    Cubic::with_options(
        CubicKind::PengRobinson,
        names,
        CubicOptions {
            parameter_reference: if tcpr { "tcPR" } else { "DEFAULT" },
            ..Default::default()
        },
    )
    .unwrap()
}

fn main() {
    let tcpr = std::env::args().any(|arg| arg == "--tcpr");
    let model_tv = model(&["CO2", "C1", "N2"], tcpr);
    let z = [0.7, 0.2, 0.1];
    let prepared = model_tv.prepare(300., &z).unwrap();
    let iterations = 500_000;
    for fresh in [false, true] {
        let start = Instant::now();
        let mut checksum = 0.;
        for i in 0..iterations {
            let volume = black_box(1e-4 + (i % 1000) as f64 * 1e-6);
            checksum += if fresh {
                black_box(&model_tv).baseline_speed_of_sound_tv(
                    black_box(300.),
                    volume,
                    black_box(&z),
                )
            } else {
                black_box(&prepared).baseline_speed_of_sound(volume)
            }
            .unwrap();
        }
        println!(
            "tv_{},{:.6},{checksum:.17e}",
            if fresh { "with_mixing" } else { "prepared" },
            start.elapsed().as_nanos() as f64 / iterations as f64
        );
    }
    let model_eq = model(&["C1", "C2", "C3"], tcpr);
    let flashes: Vec<_> = (0..100)
        .map(|i| {
            model_eq
                .baseline_tp_flash(180. + i as f64 * 0.8, 1e6, &[0.5, 0.3, 0.2])
                .unwrap()
        })
        .collect();
    let iterations = 10_000;
    let start = Instant::now();
    let mut checksum = 0.;
    for i in 0..iterations {
        checksum += black_box(&model_eq)
            .baseline_equilibrium_speed_of_sound(black_box(&flashes[i % 100]))
            .unwrap();
    }
    println!(
        "equilibrium,{:.6},{checksum:.17e}",
        start.elapsed().as_nanos() as f64 / iterations as f64
    );
}
