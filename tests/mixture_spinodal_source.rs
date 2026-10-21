//! Exact-source coordinates, raw pressure conventions and physical closure.
use thermopack_rs::mixture_spinodal::*;
use thermopack_rs::{Cubic, CubicKind, CubicOptions, Phase};

#[test]
fn full_source_mixture_points_and_curves_match() {
    let mut label = String::new();
    let mut points = Vec::new();
    let mut physical = Vec::new();
    let mut count = 0;
    let mut calls = 0;
    let mut warnings = 0;
    for row in include_str!("fixtures/mixture_spinodal.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = row.split(';').collect();
        let c: Vec<_> = f[0].split('|').collect();
        let index: usize = f[1].parse().unwrap();
        if label != f[0] {
            if !label.is_empty() {
                assert_eq!(count, points.len(), "{label}");
            }
            label = f[0].to_owned();
            count = 0;
            calls += 1;
            let kind = match c[1] {
                "SRK" => CubicKind::SoaveRedlichKwong,
                "RK" => CubicKind::RedlichKwong,
                "VDW" => CubicKind::VanDerWaals,
                _ => CubicKind::PengRobinson,
            };
            let names: Vec<_> = c[2].split(',').collect();
            let z: Vec<f64> = c[3].split(',').map(|x| x.parse().unwrap()).collect();
            let model = Cubic::with_options(
                kind,
                &names,
                CubicOptions {
                    volume_translation: c[4] == "translated",
                    ..Default::default()
                },
            )
            .unwrap();
            let p = c[6].parse().unwrap();
            points = if c[5] == "BASE_CURVE" {
                let options = SpinodalTraceOptions {
                    initial_pressure: p,
                    log_volume_step: c[7].parse().unwrap(),
                    ..Default::default()
                };
                let curve = baseline_trace_spinodal(&model, &z, options).unwrap();
                assert_ne!(curve.end, SpinodalTraceEnd::PointLimit);
                let normal = trace_spinodal(&model, &z, options).unwrap();
                physical = normal.points;
                assert_eq!(curve.points.len(), physical.len());
                curve.points
            } else {
                let point = baseline_mixture_spinodal_point(
                    &model,
                    p,
                    if c[5] == "BASE_L" {
                        Phase::Liquid
                    } else {
                        Phase::Vapor
                    },
                    &z,
                )
                .unwrap();
                physical = vec![point];
                vec![point]
            };
            for point in &physical {
                let state = model.prepare(point.temperature, &z).unwrap();
                assert!(
                    state
                        .mixture_stability(point.volume)
                        .unwrap()
                        .eigenvalue
                        .abs()
                        < 3e-10
                );
                assert!(
                    (state.pressure(point.volume).unwrap().value - point.pressure).abs() < 1e-6
                );
            }
        }
        assert_eq!(index, count);
        let point = points[index];
        let expected: Vec<f64> = f[2..6].iter().map(|x| x.parse().unwrap()).collect();
        for (j, a) in [point.temperature, point.volume, point.pressure]
            .into_iter()
            .enumerate()
        {
            assert!(
                (a - expected[j]).abs() < [2e-5, 1e-10, 2.0][j] + 2e-6 * expected[j].abs(),
                "{label} {index} {j}: {a} {}",
                expected[j]
            );
        }
        assert!(
            (physical[index].pressure - expected[3]).abs() < 2.0 + 2e-6 * expected[3].abs(),
            "{label} {index} physical pressure"
        );
        if (expected[2] - expected[3]).abs() > 1.0 + 1e-7 * expected[2].abs().max(expected[3].abs())
        {
            warnings += 1;
        }
        count += 1;
    }
    assert_eq!(count, points.len());
    assert_eq!(calls, 130);
    assert_eq!(warnings, 22);
}
