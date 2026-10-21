use thermopack_rs::mixture_spinodal::{SpinodalTraceOptions, trace_spinodal};
use thermopack_rs::spinodal_properties::*;
use thermopack_rs::{Cubic, CubicKind, CubicOptions, Phase};

fn kind(s: &str) -> CubicKind {
    match s {
        "SRK" => CubicKind::SoaveRedlichKwong,
        "RK" => CubicKind::RedlichKwong,
        "VDW" => CubicKind::VanDerWaals,
        _ => CubicKind::PengRobinson,
    }
}
fn property(code: &str) -> SpinodalProperty {
    match code {
        "1" => SpinodalProperty::Entropy,
        "2" => SpinodalProperty::Volume,
        "3" => SpinodalProperty::Enthalpy,
        "4" => SpinodalProperty::Temperature,
        "5" => SpinodalProperty::Pressure,
        _ => SpinodalProperty::InternalEnergy,
    }
}

#[test]
fn pure_source_grids_include_both_stored_critical_entries() {
    let mut label = "";
    let mut points = Vec::new();
    let mut count = 0;
    for row in include_str!("fixtures/pure_spinodal_grid.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = row.split(';').collect();
        let c: Vec<_> = f[0].split('|').collect();
        if label != f[0] {
            label = f[0];
            let model = Cubic::with_options(
                kind(c[1]),
                &[c[2]],
                CubicOptions {
                    volume_translation: c[3] == "translated",
                    ..Default::default()
                },
            )
            .unwrap();
            let spec: Vec<_> = c[4].split(':').collect();
            points = baseline_pure_spinodal_curve(
                &model,
                &[1.0],
                c[5].parse().unwrap(),
                None,
                Some(spec[1].parse().unwrap()),
            )
            .unwrap()
            .points;
            assert_eq!(points.len(), 200);
            // All noncritical vapor roots must remain away from the repulsive
            // pole even when the positive-slope mechanical bracket is narrow.
            for p in points
                .iter()
                .enumerate()
                .filter_map(|(i, p)| (![99, 100].contains(&i)).then_some(p))
            {
                let state = model.prepare(p.temperature, &[1.0]).unwrap();
                assert!(state.mixture_stability(p.volume).unwrap().eigenvalue.abs() < 1e-7);
                assert!(p.pressure.abs() < 1e10);
            }
        }
        let i: usize = f[1].parse().unwrap();
        let p = points[i];
        let expected: Vec<f64> = f[2..5].iter().map(|v| v.parse().unwrap()).collect();
        for (j, native) in [p.temperature, p.volume, p.pressure]
            .into_iter()
            .enumerate()
        {
            let base = [2e-5, 1e-11, 2.0][j] + 2e-6 * expected[j].abs();
            let mut tolerance = base;
            if [99, 100].contains(&i) && j == 1 {
                tolerance +=
                    thermopack_rs::R * expected[0] / expected[2] * (64.0 * f64::EPSILON).cbrt();
            }
            assert!(
                (native - expected[j]).abs() <= tolerance,
                "{label} {i} {j}: {native} {}",
                expected[j]
            );
        }
        count += 1;
    }
    assert_eq!(count, 4000);
}

#[test]
fn source_property_returns_match_fresh_and_cached_searches() {
    let mut count = 0;
    for row in include_str!("fixtures/spinodal_property.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = row.split(';').collect();
        let c: Vec<_> = f[0].split('|').collect();
        let q: Vec<_> = c[4].split(':').collect();
        let model = Cubic::with_options(
            kind(c[1]),
            &[c[2]],
            CubicOptions {
                volume_translation: c[3] == "translated",
                ..Default::default()
            },
        )
        .unwrap();
        let phase = match q[2] {
            "1" => Some(Phase::Liquid),
            "2" => Some(Phase::Vapor),
            _ => None,
        };
        let prop = property(q[1]);
        let target = q[3].parse().unwrap();
        let p = c[5].parse().unwrap();
        let native =
            baseline_pure_spinodal_intersection(&model, &[1.0], prop, target, p, phase, None)
                .unwrap();
        let prepared = model.prepare_pure_spinodal_search(&[1.0], p, None).unwrap();
        let cached = prepared.locate(prop, target, phase).unwrap();
        assert_eq!(native.is_some(), f[1] == "1", "{}", f[0]);
        assert_eq!(cached.is_some(), native.is_some());
        if let Some(actual) = native {
            let expected: Vec<f64> = f[2..5].iter().map(|v| v.parse().unwrap()).collect();
            let cached = cached.unwrap();
            for (i, value) in [actual.temperature, actual.volume, actual.pressure]
                .into_iter()
                .enumerate()
            {
                assert!(
                    (value - expected[i]).abs() < [2e-5, 1e-11, 2.0][i] + 2e-6 * expected[i].abs(),
                    "{} {i}: {value} {}",
                    f[0],
                    expected[i]
                );
            }
            assert!((cached.temperature - actual.temperature).abs() < 1e-7);
            assert!((cached.volume / actual.volume - 1.0).abs() < 1e-7);
        }
        count += 1;
    }
    assert_eq!(count, 276);
}

#[test]
fn baseline_cubic_volume_and_invalid_search_inputs_are_explicit() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        let model = Cubic::new(kind, &["C1"]).unwrap();
        for t in [150.0, 250.0, 400.0] {
            for p in [1e5, 1e6, 1e7] {
                for phase in [Phase::Liquid, Phase::Vapor] {
                    let state = model.prepare(t, &[1.0]).unwrap();
                    let v = state.baseline_volume(p, phase).unwrap();
                    assert!((state.pressure(v).unwrap().value - p).abs() < 1e-7 * p);
                    assert!((v / state.volume(p, phase).unwrap() - 1.0).abs() < 1e-6);
                }
            }
        }
    }
    let model = Cubic::new(CubicKind::PengRobinson, &["C1"]).unwrap();
    let search = model
        .prepare_pure_spinodal_search(&[1.0], 1e6, None)
        .unwrap();
    assert!(
        search
            .locate(SpinodalProperty::Temperature, f64::NAN, None)
            .is_err()
    );
    assert!(
        search
            .locate(
                SpinodalProperty::Temperature,
                180.0,
                Some(Phase::MinimumGibbs)
            )
            .is_err()
    );
    assert!(
        model
            .prepare_pure_spinodal_search(&[0.0], 1e6, None)
            .is_err()
    );
    assert!(baseline_pure_spinodal_curve(&model, &[1.0], 0.0, None, None).is_err());
    assert!(
        model
            .prepare(200.0, &[1.0])
            .unwrap()
            .baseline_volume(1e6, Phase::MinimumGibbs)
            .is_err()
    );
}

#[test]
fn physical_intersections_close_all_six_specifications_for_pure_and_mixture_curves() {
    for names in [vec!["C1"], vec!["C1", "C2"]] {
        let model = Cubic::new(CubicKind::PengRobinson, &names).unwrap();
        let z = vec![1.0 / names.len() as f64; names.len()];
        let curve = trace_spinodal(
            &model,
            &z,
            SpinodalTraceOptions {
                log_volume_step: 0.04,
                ..Default::default()
            },
        )
        .unwrap();
        for prop in [
            SpinodalProperty::Entropy,
            SpinodalProperty::Volume,
            SpinodalProperty::Enthalpy,
            SpinodalProperty::Temperature,
            SpinodalProperty::Pressure,
            SpinodalProperty::InternalEnergy,
        ] {
            let target = 0.5
                * (prop.evaluate(&model, &z, curve.points[5]).unwrap()
                    + prop.evaluate(&model, &z, curve.points[6]).unwrap());
            let roots = spinodal_intersections(&model, &z, &curve, prop, target).unwrap();
            assert!(!roots.is_empty(), "{names:?} {prop:?}");
            for p in roots {
                assert!(
                    model
                        .prepare(p.temperature, &z)
                        .unwrap()
                        .mixture_stability(p.volume)
                        .unwrap()
                        .eigenvalue
                        .abs()
                        < 3e-10
                );
                let value = prop.evaluate(&model, &z, p).unwrap();
                assert!(
                    (value - target).abs() < 2e-9 * target.abs().max(1e-12),
                    "{prop:?}: {value} {target}"
                );
            }
        }
        assert!(
            spinodal_intersections(&model, &z, &curve, SpinodalProperty::Temperature, -1.0)
                .unwrap()
                .is_empty()
        );
        assert!(
            spinodal_intersections(&model, &z, &curve, SpinodalProperty::Temperature, f64::NAN)
                .is_err()
        );
    }
}
