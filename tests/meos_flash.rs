use thermopack_rs::{CaloricFlashOptions, FlashPhase, LennardJonesKind, Multiparameter};

fn models() -> Vec<Multiparameter> {
    let mut models: Vec<_> = ["C3", "R134A", "N-H2", "O-H2", "P-H2"]
        .iter()
        .map(|c| Multiparameter::new(c).unwrap())
        .collect();
    for kind in [LennardJonesKind::Full, LennardJonesKind::TruncatedShifted] {
        for c in ["NE", "AR", "KR", "XE"] {
            models.push(Multiparameter::lennard_jones(c, kind).unwrap());
        }
    }
    models.push(Multiparameter::reduced_lennard_jones_truncated_shifted());
    models
}

fn case_model(case: &str) -> Multiparameter {
    let config: Vec<_> = case.split('|').nth(1).unwrap().split(':').collect();
    match config[0] {
        "LJ_MEOS" => Multiparameter::lennard_jones(config[1], LennardJonesKind::Full).unwrap(),
        "LJTS_MEOS" => {
            Multiparameter::lennard_jones(config[1], LennardJonesKind::TruncatedShifted).unwrap()
        }
        "LJTS_REDUCED" => Multiparameter::reduced_lennard_jones_truncated_shifted(),
        _ => Multiparameter::new(config[1]).unwrap(),
    }
}

#[test]
fn coexistence_tp_returns_a_stable_endpoint_without_quality_claim() {
    let mut count = 0;
    for row in include_str!("fixtures/meos_flash_coexistence_tp.csv")
        .lines()
        .skip(1)
    {
        let fields: Vec<_> = row.split(';').collect();
        let model = case_model(fields[0]);
        let t = fields[1].parse().unwrap();
        let p = fields[2].parse().unwrap();
        let f = model.tp_flash(t, p).unwrap();
        let q = model.equilibrium_properties(&f).unwrap();
        assert_ne!(f.phase, FlashPhase::TwoPhase);
        let prep = model.prepare(t).unwrap();
        let pressure = prep.pressure(q.volume).unwrap();
        assert!(pressure.dv < 0.);
        assert!((pressure.value - p).abs() < 1e-4 + 2e-8 * p);
        let sat = model.saturation_pressure(t).unwrap();
        let member = [sat.liquid_volume, sat.vapor_volume]
            .iter()
            .any(|v| (q.volume / v - 1.).abs() < 2e-6);
        assert!(
            member,
            "{}: {} is not a coexistence endpoint",
            fields[0], q.volume
        );
        count += 1;
    }
    assert_eq!(count, 46);
}

fn options(model: &Multiparameter, guess: f64) -> CaloricFlashOptions {
    CaloricFlashOptions {
        temperature_guess: Some(guess),
        minimum_temperature: 0.3 * model.tc(),
        maximum_temperature: 3. * model.tc(),
        ..Default::default()
    }
}

#[test]
fn homogeneous_caloric_roundtrips_and_closure_across_model_families() {
    for model in models() {
        let t = 1.2 * model.tc();
        let pc = if model.pc() > 0. {
            model.pc()
        } else {
            model
                .prepare(model.tc())
                .unwrap()
                .pressure(1. / model.rc())
                .unwrap()
                .value
        };
        let p = 1.1 * pc;
        let flash = model.tp_flash(t, p).unwrap();
        let state = model.equilibrium_properties(&flash).unwrap();
        assert!(
            (state.enthalpy - state.internal_energy - p * state.volume).abs()
                < 1e-8 * model.gas_constant() * t
        );
        for enthalpy in [false, true] {
            let result = if enthalpy {
                model.ph_flash_with_options(p, state.enthalpy, options(&model, t * 1.2))
            } else {
                model.ps_flash_with_options(p, state.entropy, options(&model, t * 1.2))
            }
            .unwrap_or_else(|e| panic!("{} H={enthalpy}: {e}", model.record().ident));
            assert!(
                (result.temperature - t).abs() < 2e-8 * t,
                "{}: {} != {t}",
                model.record().ident,
                result.temperature
            );
            assert_eq!(result.pressure, p);
            assert_eq!(result.liquid, [1.]);
            assert_eq!(result.vapor, [1.]);
        }
    }
}

#[test]
fn latent_splits_bypass_homogeneous_temperature_bounds() {
    for model in models() {
        let t = 0.8 * model.tc();
        let seed = model.saturated_flash(t, 0.3).unwrap();
        let q = model.equilibrium_properties(&seed).unwrap();
        let bounds = CaloricFlashOptions {
            minimum_temperature: 2. * model.tc(),
            maximum_temperature: 3. * model.tc(),
            ..Default::default()
        };
        for result in [
            model.ph_flash_with_options(seed.pressure, q.enthalpy, bounds),
            model.ps_flash_with_options(seed.pressure, q.entropy, bounds),
        ] {
            let f = result.unwrap_or_else(|e| panic!("{}: {e}", model.record().ident));
            assert_eq!(f.phase, FlashPhase::TwoPhase);
            assert_eq!(f.pressure, seed.pressure);
            assert!((f.temperature - t).abs() < 2e-8 * t);
            assert!((f.vapor_fraction - 0.3).abs() < 2e-8);
        }
    }
}

#[test]
fn reduced_default_bounds_include_normal_reduced_temperatures() {
    let model = Multiparameter::reduced_lennard_jones_truncated_shifted();
    let f = model.tp_flash(1.4, 0.2).unwrap();
    let q = model.equilibrium_properties(&f).unwrap();
    assert!((model.ph_flash(0.2, q.enthalpy).unwrap().temperature - 1.4).abs() < 1e-8);
    assert!((model.ps_flash(0.2, q.entropy).unwrap().temperature - 1.4).abs() < 1e-8);
}

#[test]
fn recomputed_properties_ignore_cached_volumes_and_reject_invalid_splits() {
    let model = Multiparameter::new("C3").unwrap();
    let mut f = model.tp_flash(450., 1e6).unwrap();
    let q = model.equilibrium_properties(&f).unwrap();
    f.liquid_volume = Some(-1.);
    f.vapor_volume = Some(f64::NAN);
    assert_eq!(model.equilibrium_properties(&f).unwrap().volume, q.volume);
    f.feed = vec![0.5];
    assert!(model.equilibrium_properties(&f).is_err());
    let mut f = model.saturated_flash(300., 0.5).unwrap();
    f.liquid_fraction = 0.3;
    assert!(model.equilibrium_properties(&f).is_err());
    f.liquid_fraction = 0.5;
    f.pressure *= 1.1;
    assert!(model.equilibrium_properties(&f).is_err());
    for t in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(model.tp_flash(t, 1e5).is_err());
    }
    for p in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(model.ph_flash(p, 100.).is_err());
    }
    for beta in [-0.1, 1.1, f64::NAN] {
        assert!(model.saturated_flash(300., beta).is_err());
    }
}

fn check_fixture(fixture: &str, expected_count: usize) {
    let mut count = 0;
    for row in fixture.lines().skip(1) {
        let fields: Vec<_> = row.split(';').collect();
        let case: Vec<_> = fields[0].split('|').collect();
        let config: Vec<_> = case[1].split(':').collect();
        let model = match config[0] {
            "LJ_MEOS" => Multiparameter::lennard_jones(config[1], LennardJonesKind::Full).unwrap(),
            "LJTS_MEOS" => {
                Multiparameter::lennard_jones(config[1], LennardJonesKind::TruncatedShifted)
                    .unwrap()
            }
            "LJTS_REDUCED" => Multiparameter::reduced_lennard_jones_truncated_shifted(),
            _ => Multiparameter::new(config[1]).unwrap(),
        };
        let control: Vec<_> = case[5].split(':').collect();
        let t = case[2].parse().unwrap();
        let p = case[3].parse().unwrap();
        let options = CaloricFlashOptions {
            temperature_guess: if control[3] == "none" {
                None
            } else {
                Some(control[3].parse().unwrap())
            },
            minimum_temperature: control[4].parse().unwrap(),
            maximum_temperature: control[5].parse().unwrap(),
            ..Default::default()
        };
        let f = match control[1] {
            "TP" => model.tp_flash(t, p),
            "PH" => model.ph_flash_with_options(p, control[2].parse().unwrap(), options),
            "PS" => model.ps_flash_with_options(p, control[2].parse().unwrap(), options),
            _ => panic!("unknown flash {}", fields[0]),
        }
        .unwrap_or_else(|e| panic!("{}: {e}", fields[0]));
        let q = model
            .equilibrium_properties(&f)
            .unwrap_or_else(|e| panic!("{}: {e}", fields[0]));
        let two = f.phase == FlashPhase::TwoPhase
            && f.vapor_fraction > 1e-8
            && f.vapor_fraction < 1. - 1e-8;
        let actual = [
            f.temperature,
            f.pressure,
            f64::from(two),
            if two { f.vapor_fraction } else { 0. },
            if two {
                f.liquid_volume.unwrap()
            } else {
                q.volume
            },
            if two {
                f.vapor_volume.unwrap()
            } else {
                q.volume
            },
            q.enthalpy,
            q.entropy,
        ];
        for (i, a) in actual.into_iter().enumerate() {
            let b: f64 = fields[i + 1].parse().unwrap();
            let floor = [2e-7, 1e-4, 0., 2e-7, 1e-12, 1e-12, 2e-5, 2e-7][i];
            assert!(
                (a - b).abs() <= floor + 2e-6 * a.abs().max(b.abs()),
                "{} field {i}: {a} != {b}",
                fields[0]
            );
        }
        count += 1;
    }
    assert_eq!(count, expected_count);
}

#[test]
fn qualified_source_flash_properties_match() {
    check_fixture(include_str!("fixtures/meos_flash.csv"), 1708);
}

#[test]
fn source_kernel_qualified_caloric_recoveries_match() {
    check_fixture(include_str!("fixtures/meos_flash_recovered.csv"), 100);
}
