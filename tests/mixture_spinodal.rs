use thermopack_rs::mixture_spinodal::*;
use thermopack_rs::spinodal::spinodal_point;
use thermopack_rs::{Cubic, CubicKind, Phase};

#[test]
fn local_mixture_limits_and_curves_close() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        let model = Cubic::new(kind, &["C1", "C2"]).unwrap();
        let z = [0.5, 0.5];
        for p in [1e5, 1e6, 3e6] {
            for phase in [Phase::Liquid, Phase::Vapor] {
                let point = spinodal_point(&model, p, phase, &z).unwrap();
                let state = model.prepare(point.temperature, &z).unwrap();
                assert!(
                    state
                        .mixture_stability(point.volume)
                        .unwrap()
                        .eigenvalue
                        .abs()
                        < 3e-10
                );
                assert!((state.pressure(point.volume).unwrap().value - p).abs() < 3e-10 * p);
            }
        }
        let curve = trace_spinodal(&model, &z, Default::default()).unwrap();
        assert_ne!(curve.end, SpinodalTraceEnd::PointLimit);
        for p in &curve.points {
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
        }
        assert!(curve.points.windows(2).all(|w| w[1].volume > w[0].volume));
    }
}

#[test]
fn explicit_seeds_endpoints_and_source_branch_defect_are_separate() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
    let z = [0.5, 0.5];
    let liquid = spinodal_point(&model, 1e6, Phase::Liquid, &z).unwrap();
    let vapor = spinodal_point(&model, 1e6, Phase::Vapor, &z).unwrap();
    let source_vapor = baseline_mixture_spinodal_point(&model, 1e6, Phase::Vapor, &z).unwrap();
    assert!((source_vapor.volume / liquid.volume - 1.0).abs() < 1e-7);
    assert!(vapor.volume > 5.0 * liquid.volume);
    assert!((vapor.temperature - 154.022266526).abs() < 1e-6);
    let corrected = spinodal_point_with_guess(
        &model,
        1e6,
        &z,
        vapor.temperature * 1.01,
        vapor.volume * 1.02,
    )
    .unwrap();
    assert!((corrected.volume / vapor.volume - 1.0).abs() < 1e-8);
    let inverse =
        spinodal_volume_at_temperature(&model, vapor.temperature, &z, vapor.volume * 1.02).unwrap();
    let inverse_t =
        spinodal_temperature_at_volume(&model, vapor.volume, &z, vapor.temperature * 1.02).unwrap();
    assert!((inverse.volume / vapor.volume - 1.0).abs() < 1e-8);
    assert!((inverse_t.temperature - vapor.temperature).abs() < 1e-6);
    let state = model.prepare(liquid.temperature, &z).unwrap();
    assert!(state.pressure(liquid.volume).unwrap().dv < -1e8);
    let low = model
        .prepare(liquid.temperature - 0.01, &z)
        .unwrap()
        .mixture_stability(liquid.volume)
        .unwrap();
    let high = model
        .prepare(liquid.temperature + 0.01, &z)
        .unwrap()
        .mixture_stability(liquid.volume)
        .unwrap();
    assert!(low.eigenvalue < 0.0 && high.eigenvalue > 0.0);
    let physical = trace_spinodal(&model, &z, Default::default()).unwrap();
    let baseline = baseline_trace_spinodal(&model, &z, Default::default()).unwrap();
    assert_eq!(physical.end, SpinodalTraceEnd::MinimumTemperature);
    assert_eq!(physical.points.len(), baseline.points.len());
    let a = physical.points.last().unwrap();
    let b = baseline.points.last().unwrap();
    assert!((a.temperature - b.temperature).abs() < 1e-6);
    assert!((a.pressure - b.pressure).abs() > 1.0);
    let pressure_end = trace_spinodal(
        &model,
        &z,
        SpinodalTraceOptions {
            minimum_vapor_temperature: Some(20.0),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(pressure_end.end, SpinodalTraceEnd::InitialPressure);
    assert!((pressure_end.points.last().unwrap().pressure - 1e5).abs() < 1e-3);
    let truncated = trace_spinodal(
        &model,
        &z,
        SpinodalTraceOptions {
            maximum_points: 2,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(truncated.end, SpinodalTraceEnd::PointLimit);
    assert_eq!(truncated.points.len(), 2);
}

#[test]
fn embedded_pure_translation_and_invalid_inputs() {
    let pure = Cubic::new(CubicKind::PengRobinson, &["C2"]).unwrap();
    let mut embedded = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
    embedded.enable_volume_translation("DEFAULT").unwrap();
    let mut shifted = pure.clone();
    shifted.enable_volume_translation("DEFAULT").unwrap();
    for phase in [Phase::Liquid, Phase::Vapor] {
        let a = spinodal_point(&shifted, 1e6, phase, &[1.0]).unwrap();
        let b = spinodal_point(&embedded, 1e6, phase, &[0.0, 1.0]).unwrap();
        assert!((a.temperature - b.temperature).abs() < 1e-7);
        assert!((a.volume / b.volume - 1.0).abs() < 1e-10);
    }
    let curve = trace_spinodal(&pure, &[1.0], Default::default()).unwrap();
    assert_ne!(curve.end, SpinodalTraceEnd::PointLimit);
    assert!(baseline_trace_spinodal(&pure, &[1.0], Default::default()).is_err());
    for z in [
        vec![1.0],
        vec![0.0, 0.0],
        vec![-0.5, 1.5],
        vec![f64::NAN, 1.0],
    ] {
        assert!(spinodal_point(&embedded, 1e6, Phase::Liquid, &z).is_err());
    }
    assert!(spinodal_point(&pure, 1e6, Phase::MinimumGibbs, &[1.0]).is_err());
    for options in [
        SpinodalTraceOptions {
            maximum_points: 1,
            ..Default::default()
        },
        SpinodalTraceOptions {
            log_volume_step: 0.0,
            ..Default::default()
        },
        SpinodalTraceOptions {
            log_volume_step: f64::NAN,
            ..Default::default()
        },
        SpinodalTraceOptions {
            minimum_vapor_temperature: Some(-1.0),
            ..Default::default()
        },
    ] {
        assert!(trace_spinodal(&pure, &[1.0], options).is_err());
    }
}
