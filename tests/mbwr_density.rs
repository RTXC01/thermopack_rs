use thermopack_rs::{Mbwr, MbwrDensityParameters, MbwrDensityPath, MbwrKind, Phase, R};

#[test]
fn owned_preparation_and_custom_ideal_density_have_source_phase_flags() {
    let parameters = MbwrDensityParameters::new(190., 4.6e6, 10000., 0.01).unwrap();
    for kind in [MbwrKind::Bender19, MbwrKind::Mbwr32] {
        let prepared = {
            let model =
                Mbwr::with_parameters(kind, &vec![0.; kind.coefficient_count()], 0.1).unwrap();
            model.prepare_baseline_density(300., parameters).unwrap()
        };
        let raw = prepared.source_parameters();
        assert_eq!(raw.len(), if kind == MbwrKind::Bender19 { 9 } else { 16 });
        assert_eq!(raw[0], 0.1);
        assert_eq!(raw[1], 1000. * R * 300.);
        assert!(raw[2..].iter().all(|&v| v == 0.));
        for p in [0.01, 100., 1e5, 1e8] {
            for phase in [Phase::Liquid, Phase::Vapor] {
                let result = prepared.density(p, phase).unwrap();
                assert_eq!(result.phase, phase);
                assert_eq!(result.path, MbwrDensityPath::RequestedPhaseNewton);
                assert!(result.pressure_converged && result.mechanically_stable);
                assert!((result.density - p / (R * 300.)).abs() < 1e-12 + 1e-12 * result.density);
                assert!((result.pressure_density_derivative - R * 300.).abs() < 1e-10);
            }
        }
    }
}

#[test]
fn tabulated_initialization_and_model_convenience_match_prepared_queries() {
    for kind in [MbwrKind::Bender19, MbwrKind::Mbwr32] {
        let model = Mbwr::new(kind, "C1").unwrap();
        let data = model.source_density_parameters().unwrap();
        let critical = model.source_critical_point().unwrap();
        assert_eq!(data.critical_temperature, critical.temperature);
        assert_eq!(data.critical_pressure, critical.pressure);
        assert_eq!(data.critical_density, 1. / critical.volume);
        assert_eq!(data.acentric_factor, 0.01131);
        let prepared = model.prepare_baseline_density(150., data).unwrap();
        for phase in [Phase::Liquid, Phase::Vapor] {
            let a = model.baseline_density(150., 1e6, phase, data).unwrap();
            let b = prepared.density(1e6, phase).unwrap();
            assert_eq!(a.density, b.density);
            assert_eq!(a.phase, b.phase);
            assert_eq!(a.path, b.path);
            assert!(a.pressure_converged && a.mechanically_stable);
        }
    }
    let no_table = Mbwr::new(MbwrKind::Mbwr32, "IC4").unwrap();
    assert!(no_table.source_density_parameters().is_err());
    let no_acf = Mbwr::from_record("mbwrcx12").unwrap();
    assert!(no_acf.source_density_parameters().is_err());
}

#[test]
fn invalid_domains_and_bypassed_metadata_constructor_are_rejected() {
    let model = Mbwr::new(MbwrKind::Mbwr32, "C1").unwrap();
    let data = model.source_density_parameters().unwrap();
    for value in [0., -1., f64::INFINITY, f64::NAN] {
        assert!(MbwrDensityParameters::new(value, 4e6, 10000., 0.).is_err());
        assert!(MbwrDensityParameters::new(200., value, 10000., 0.).is_err());
        assert!(MbwrDensityParameters::new(200., 4e6, value, 0.).is_err());
        assert!(model.prepare_baseline_density(value, data).is_err());
        assert!(
            model
                .baseline_density(150., value, Phase::Vapor, data)
                .is_err()
        );
    }
    assert!(MbwrDensityParameters::new(200., 4e6, 10000., f64::NAN).is_err());
    assert!(
        model
            .prepare_baseline_density(
                150.,
                MbwrDensityParameters {
                    critical_pressure: 0.,
                    ..data
                }
            )
            .is_err()
    );
    assert!(
        model
            .baseline_density(150., 1e6, Phase::MinimumGibbs, data)
            .is_err()
    );
    let custom = Mbwr::with_parameters(MbwrKind::Bender19, &[0.; 19], 0.1).unwrap();
    assert!(custom.source_density_parameters().is_err());
}
