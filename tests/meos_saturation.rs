use thermopack_rs::{LennardJonesKind, Multiparameter};

fn models() -> Vec<Multiparameter> {
    let mut m: Vec<_> = ["C3", "R134A", "N-H2", "O-H2", "P-H2"]
        .iter()
        .map(|c| Multiparameter::new(c).unwrap())
        .collect();
    for kind in [LennardJonesKind::Full, LennardJonesKind::TruncatedShifted] {
        for c in ["NE", "AR", "KR", "XE"] {
            m.push(Multiparameter::lennard_jones(c, kind).unwrap());
        }
    }
    m.push(Multiparameter::reduced_lennard_jones_truncated_shifted());
    m
}

fn case_model(case: &str) -> Multiparameter {
    let fields: Vec<_> = case.split('|').collect();
    let config: Vec<_> = fields[1].split(':').collect();
    match config[0] {
        "LJ_MEOS" => Multiparameter::lennard_jones(config[1], LennardJonesKind::Full).unwrap(),
        "LJTS_MEOS" => {
            Multiparameter::lennard_jones(config[1], LennardJonesKind::TruncatedShifted).unwrap()
        }
        "LJTS_REDUCED" => Multiparameter::reduced_lennard_jones_truncated_shifted(),
        "NIST_MEOS" => Multiparameter::new(config[1]).unwrap(),
        _ => panic!("unrecognized model {case}"),
    }
}

fn check_fixture(fixture: &str, expected_count: usize) {
    let mut count = 0;
    for row in fixture.lines().skip(1) {
        let fields: Vec<_> = row.split(';').collect();
        let case: Vec<_> = fields[0].split('|').collect();
        let model = case_model(fields[0]);
        let input = case[2].parse().unwrap();
        let point = match case[5] {
            "SAT_P" => model.saturation_pressure(input),
            "SAT_T" => model.saturation_temperature(input),
            _ => panic!("unknown saturation mode"),
        }
        .unwrap_or_else(|e| panic!("{}: {e}", fields[0]));
        let actual = [
            point.temperature,
            point.pressure,
            point.liquid_volume,
            point.vapor_volume,
            point.dpdt,
        ];
        for (i, value) in actual.into_iter().enumerate() {
            let expected: f64 = fields[i + 1].parse().unwrap();
            let floor = [2e-7, 1e-2, 1e-12, 1e-12, 1e-4][i];
            assert!(
                (value - expected).abs() < floor + 2e-6 * value.abs().max(expected.abs()),
                "{} field {i}: {value} != {expected}",
                fields[0]
            );
        }
        count += 1;
    }
    assert_eq!(count, expected_count);
}

#[test]
fn qualified_public_source_saturation_coordinates_match() {
    check_fixture(include_str!("fixtures/meos_saturation.csv"), 164);
}

#[test]
fn source_kernel_qualified_recoveries_remain_reproducible() {
    // These are native solutions checked against source kernels, rather than
    // coordinates successfully returned by the source saturation solver.
    check_fixture(include_str!("fixtures/meos_saturation_recovered.csv"), 68);
}

#[test]
fn analytic_clapeyron_slope_matches_curve_differences() {
    for model in models() {
        for tr in [0.8, 0.99] {
            let t = tr * model.tc();
            let h = 1e-5 * t;
            let point = model.saturation_pressure(t).unwrap();
            let hi = model.saturation_pressure(t + h).unwrap();
            let lo = model.saturation_pressure(t - h).unwrap();
            let finite = (hi.pressure - lo.pressure) / (2. * h);
            assert!(
                (point.dpdt - finite).abs() < 2e-8 + 3e-6 * finite.abs(),
                "{} {tr}: {} != {finite}",
                model.record().ident,
                point.dpdt
            );
        }
    }
}

#[test]
fn invalid_inputs_and_unresolved_critical_endpoint_are_rejected() {
    for model in models() {
        for t in [
            0.,
            -1.,
            f64::NAN,
            f64::INFINITY,
            model.tc(),
            1.1 * model.tc(),
        ] {
            assert!(model.saturation_pressure(t).is_err());
        }
        let pc = model
            .prepare(model.tc())
            .unwrap()
            .pressure(1. / model.rc())
            .unwrap()
            .value;
        for p in [0., -1., f64::NAN, f64::INFINITY, pc, 1.1 * pc] {
            assert!(model.saturation_temperature(p).is_err());
        }
        let point = model.saturation_pressure(0.8 * model.tc()).unwrap();
        for t in [0., -1., f64::NAN, model.tc()] {
            assert!(
                model
                    .saturation_temperature_with_guess(point.pressure, t)
                    .is_err()
            );
        }
    }
}

#[test]
fn inverse_temperature_guesses_on_both_sides_recover_the_same_curve() {
    for model in models() {
        for tr in [0.5, 0.9, 0.999] {
            let expected = model.saturation_pressure(tr * model.tc()).unwrap();
            for guess in [0.55, 0.95] {
                let recovered = model
                    .saturation_temperature_with_guess(expected.pressure, guess * model.tc())
                    .unwrap_or_else(|e| {
                        panic!("{} target {tr}, guess {guess}: {e}", model.record().ident)
                    });
                assert!(
                    (recovered.temperature - expected.temperature).abs()
                        < 1e-7 * expected.temperature
                );
            }
        }
    }
}

#[test]
fn all_model_coexistence_curves_close_pressure_mu_latent_heat_and_inverse() {
    for model in models() {
        for tr in [0.45, 0.5, 0.6, 0.75, 0.9, 0.95, 0.99, 0.999, 0.9999] {
            let t = tr * model.tc();
            let sat = model.saturation_pressure(t).unwrap_or_else(|e| {
                panic!(
                    "{} {:?} {tr}: {e}",
                    model.record().ident,
                    model.lennard_jones_parameters()
                )
            });
            let s = model.prepare(t).unwrap();
            let l = s.state(sat.liquid_volume).unwrap();
            let v = s.state(sat.vapor_volume).unwrap();
            assert!(sat.liquid_volume < 1. / model.rc() && sat.vapor_volume > 1. / model.rc());
            assert!(l.pressure.dv < 0. && v.pressure.dv < 0.);
            assert!(
                (l.pressure.value - sat.pressure).abs()
                    < 1e-6 * sat.pressure + 1e-7 * model.gas_constant() * t * model.rc()
            );
            assert!((v.pressure.value - sat.pressure).abs() < 1e-8 * sat.pressure);
            assert!((l.ln_fugacity.unwrap() - v.ln_fugacity.unwrap()).abs() < 1e-6);
            let latent = v.enthalpy - l.enthalpy;
            assert!((latent - t * (v.entropy - l.entropy)).abs() < 1e-6 * latent.abs());
            assert!(
                (sat.dpdt - latent / (t * (sat.vapor_volume - sat.liquid_volume))).abs()
                    < 1e-7 * sat.dpdt
            );
            let recovered = model
                .saturation_temperature(sat.pressure)
                .unwrap_or_else(|e| panic!("inverse {} {tr}: {e}", model.record().ident));
            assert!((recovered.temperature - t).abs() < 1e-7 * t);
        }
    }
}
