use std::{hint::black_box, time::Instant};
use thermopack_rs::pure_saturation_curve::{
    PureSaturationEnd, PureSaturationOptions, PureSaturationStart,
};
use thermopack_rs::{Cubic, CubicKind, CubicOptions};

fn main() {
    let check = std::env::args().any(|s| s == "--check");
    let envelope = std::env::args().any(|s| s == "--envelope");
    for profile in 1..=4 {
        let model = match profile {
            1 => Cubic::new(CubicKind::PengRobinson, &["C1"]).unwrap(),
            2 => Cubic::new(CubicKind::SoaveRedlichKwong, &["CO2"]).unwrap(),
            3 => Cubic::with_options(
                CubicKind::PengRobinson,
                &["H2O"],
                CubicOptions {
                    volume_translation: true,
                    ..Default::default()
                },
            )
            .unwrap(),
            _ => Cubic::quantum(&["H2"]).unwrap(),
        };
        let c = &model.components()[0];
        let options = PureSaturationOptions {
            start: if profile == 2 {
                PureSaturationStart::Temperature(0.65 * c.critical_temperature)
            } else {
                PureSaturationStart::Pressure(0.05 * c.critical_pressure)
            },
            logarithmic: profile == 2 || profile == 4,
            ..Default::default()
        };
        if envelope {
            run_envelope(profile, &model, options, check);
            continue;
        }
        let solve = || {
            let curve = black_box(&model)
                .baseline_pure_saturation_curve(0, black_box(options))
                .unwrap();
            assert_eq!(curve.end, PureSaturationEnd::SourceCriticalParameters);
            curve
        };
        if check {
            for (i, p) in solve().points.iter().enumerate() {
                println!(
                    "check_{profile}_{i},{:.17e},{:.17e},{:.17e},{:.17e}",
                    p.temperature, p.pressure, p.liquid_volume, p.vapor_volume
                );
            }
            continue;
        }
        let iterations = if cfg!(debug_assertions) { 1 } else { 20 };
        let mut sums = [0.; 5];
        let start = Instant::now();
        for _ in 0..iterations {
            let curve = black_box(solve());
            sums[0] += curve.points.len() as f64;
            for p in curve.points {
                sums[1] += p.temperature;
                sums[2] += p.pressure;
                sums[3] += p.liquid_volume;
                sums[4] += p.vapor_volume;
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

fn run_envelope(profile: usize, model: &Cubic, options: PureSaturationOptions, check: bool) {
    use thermopack_rs::SaturationKind;
    use thermopack_rs::envelope::{EnvelopeConstraint, EnvelopeEnd, EnvelopeOptions};
    let solve = || {
        let curve = black_box(model)
            .baseline_phase_envelope(
                black_box(&[1.]),
                EnvelopeOptions {
                    start: match options.start {
                        PureSaturationStart::Temperature(t) => EnvelopeConstraint::Temperature(t),
                        PureSaturationStart::Pressure(p) => EnvelopeConstraint::Pressure(p),
                    },
                    maximum_points: 100,
                    ..Default::default()
                },
            )
            .unwrap();
        assert_eq!(curve.end, EnvelopeEnd::SourceCriticalParameters);
        assert!(curve.points.iter().all(|p| p.kind == SaturationKind::Bubble
            && p.liquid == [1.]
            && p.vapor == [1.]
            && p.fugacity_residual.is_finite()));
        curve
    };
    if check {
        for (i, p) in solve().points.iter().enumerate() {
            println!(
                "check_{profile}_{i},{:.17e},{:.17e},{:.17e},{:.17e}",
                p.temperature, p.pressure, p.liquid_volume, p.vapor_volume
            );
        }
        return;
    }
    let iterations = if cfg!(debug_assertions) { 1 } else { 20 };
    let start = Instant::now();
    let mut sums = [0.; 5];
    for _ in 0..iterations {
        let curve = black_box(solve());
        sums[0] += curve.points.len() as f64;
        for p in curve.points {
            sums[1] += p.temperature;
            sums[2] += p.pressure;
            sums[3] += p.liquid_volume;
            sums[4] += p.vapor_volume;
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
