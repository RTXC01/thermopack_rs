use thermopack_rs::{LennardJonesKind, Multiparameter};

fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1e-12),
        "{a:e} != {b:e}"
    );
}

fn models() -> Vec<Multiparameter> {
    let mut result: Vec<_> = ["C3", "R134A", "N-H2", "O-H2", "P-H2"]
        .into_iter()
        .map(|name| Multiparameter::new(name).unwrap())
        .collect();
    for kind in [LennardJonesKind::Full, LennardJonesKind::TruncatedShifted] {
        for gas in ["NE", "AR", "KR", "XE"] {
            result.push(Multiparameter::lennard_jones(gas, kind).unwrap());
        }
    }
    result.push(Multiparameter::reduced_lennard_jones_truncated_shifted());
    result
}

#[test]
fn fixed_volume_cv_matches_energy_and_entropy_differences_across_variants() {
    for model in models() {
        for tr in [0.5, 0.8, 0.95, 0.99, 0.999] {
            let t = tr * model.tc();
            let sat = model.saturation_pressure(t).unwrap();
            for beta in [0.2, 0.5, 0.8] {
                let volume = (1. - beta) * sat.liquid_volume + beta * sat.vapor_volume;
                let state = model.equilibrium_tv(t, volume).unwrap();
                assert_eq!(state.flash.phase, thermopack_rs::FlashPhase::TwoPhase);
                let h = (2e-6 * t).min(0.001 * (model.tc() - t));
                let up = model.equilibrium_tv(t + h, volume).unwrap();
                let lo = model.equilibrium_tv(t - h, volume).unwrap();
                close(
                    state.heat_capacity_v,
                    (up.properties.internal_energy - lo.properties.internal_energy) / (2. * h),
                    3e-4,
                );
                close(
                    state.heat_capacity_v / t,
                    (up.properties.entropy - lo.properties.entropy) / (2. * h),
                    3e-4,
                );
            }
        }
    }
}

#[test]
fn coexistence_endpoints_and_input_boundaries() {
    let model = Multiparameter::new("R134A").unwrap();
    let t = 0.78 * model.tc();
    let sat = model.saturation_pressure(t).unwrap();
    for (volume, beta) in [(sat.liquid_volume, 0.), (sat.vapor_volume, 1.)] {
        let state = model.equilibrium_tv(t, volume).unwrap();
        assert_eq!(state.flash.phase, thermopack_rs::FlashPhase::TwoPhase);
        close(state.flash.vapor_fraction, beta, 1e-12);
    }
    assert!(model.equilibrium_tv(f64::NAN, sat.liquid_volume).is_err());
    assert!(model.equilibrium_tv(t, f64::INFINITY).is_err());
    assert!(model.uv_flash(f64::NAN, sat.liquid_volume).is_err());
    assert!(model.sv_flash(1., 0.).is_err());
}

#[test]
fn dense_cold_targets_can_start_above_the_pressure_limit() {
    for model in [
        Multiparameter::new("C3").unwrap(),
        Multiparameter::new("R134A").unwrap(),
    ] {
        let t = 0.5 * model.tc();
        let volume = model.saturation_pressure(t).unwrap().liquid_volume;
        let target = model.equilibrium_tv(t, volume).unwrap();
        let initial = model.equilibrium_tv(model.tc(), volume).unwrap();
        let options = thermopack_rs::VolumeFlashOptions {
            minimum_temperature: 2.,
            ..model.volume_flash_options()
        };
        assert!(initial.flash.pressure > options.maximum_pressure);
        let uv = model
            .uv_flash_with_options(target.properties.internal_energy, volume, options)
            .unwrap();
        let sv = model
            .sv_flash_with_options(target.properties.entropy, volume, options)
            .unwrap();
        close(uv.temperature, t, 2e-8);
        close(sv.temperature, t, 2e-8);
        assert!(uv.pressure <= options.maximum_pressure);
        let rejected = thermopack_rs::VolumeFlashOptions {
            maximum_pressure: 0.9 * target.flash.pressure,
            ..options
        };
        assert!(
            model
                .uv_flash_with_options(target.properties.internal_energy, volume, rejected)
                .is_err()
        );
        assert!(
            model
                .sv_flash_with_options(target.properties.entropy, volume, rejected)
                .is_err()
        );
    }
}

#[test]
fn uv_component_triple_temperature_repairs_seeds_without_clipping_the_root() {
    for kind in [LennardJonesKind::Full, LennardJonesKind::TruncatedShifted] {
        let model = Multiparameter::lennard_jones("KR", kind).unwrap();
        let t = 0.5 * model.tc();
        let triple = thermopack_rs::component("KR")
            .unwrap()
            .source_triple_temperature();
        assert!(t < triple && t > model.volume_flash_options().minimum_temperature);
        let sat = model.saturation_pressure(t).unwrap();
        let volume = 0.5 * (sat.liquid_volume + sat.vapor_volume);
        let target = model
            .equilibrium_tv(t, volume)
            .unwrap()
            .properties
            .internal_energy;
        let cold = model.uv_flash(target, volume).unwrap();
        let guessed = model
            .uv_flash_with_options(
                target,
                volume,
                thermopack_rs::VolumeFlashOptions {
                    temperature_guess: Some(0.93 * t),
                    ..model.volume_flash_options()
                },
            )
            .unwrap();
        close(cold.temperature, t, 2e-8);
        close(guessed.temperature, t, 2e-8);
    }
}

#[test]
fn uv_sv_roundtrip_homogeneous_and_latent_states() {
    for model in models() {
        let (t, volume) = (1.2 * model.tc(), 3. / model.rc());
        let state = model.equilibrium_tv(t, volume).unwrap();
        let mut options = model.volume_flash_options();
        options.minimum_temperature = options.minimum_temperature.min(0.5 * t);
        let uv = model
            .uv_flash_with_options(state.properties.internal_energy, volume, options)
            .unwrap();
        let sv = model
            .sv_flash_with_options(state.properties.entropy, volume, options)
            .unwrap();
        close(uv.temperature, t, 2e-7);
        close(sv.temperature, t, 2e-7);
        let mut pressure_seeded = options;
        pressure_seeded.pressure_guess = Some(0.5 * state.flash.pressure);
        let uv_seeded = model
            .uv_flash_with_options(state.properties.internal_energy, volume, pressure_seeded)
            .unwrap();
        close(uv_seeded.temperature, uv.temperature, 1e-12);
        close(uv_seeded.pressure, uv.pressure, 1e-12);

        {
            let sat_t = 0.82 * model.tc();
            let sat = model.saturation_pressure(sat_t).unwrap();
            let beta = 0.41;
            let volume = (1. - beta) * sat.liquid_volume + beta * sat.vapor_volume;
            let latent = model.equilibrium_tv(sat_t, volume).unwrap();
            let mut options = model.volume_flash_options();
            options.minimum_temperature = options.minimum_temperature.min(0.5 * sat_t);
            let recovered = model
                .uv_flash_with_options(latent.properties.internal_energy, volume, options)
                .unwrap();
            close(recovered.temperature, sat_t, 3e-7);
            close(recovered.vapor_fraction, beta, 2e-6);
            let entropy_recovered = model
                .sv_flash_with_options(latent.properties.entropy, volume, options)
                .unwrap();
            close(entropy_recovered.temperature, sat_t, 3e-7);
            close(entropy_recovered.vapor_fraction, beta, 2e-6);
        }
    }
}

fn fixture_flash(case: &str) -> (Multiparameter, thermopack_rs::FlashResult) {
    let f: Vec<_> = case.split('|').collect();
    assert_eq!(f.len(), 6);
    let config: Vec<_> = f[1].split(':').collect();
    let model = match config[0] {
        "LJ_MEOS" => Multiparameter::lennard_jones(config[1], LennardJonesKind::Full).unwrap(),
        "LJTS_MEOS" => {
            Multiparameter::lennard_jones(config[1], LennardJonesKind::TruncatedShifted).unwrap()
        }
        "LJTS_REDUCED" => Multiparameter::reduced_lennard_jones_truncated_shifted(),
        _ => Multiparameter::new(config[1]).unwrap(),
    };
    let c: Vec<_> = f[5].split(':').collect();
    assert_eq!(c.len(), 9);
    assert_eq!(c[0], "ISO");
    let guess = |s: &str| {
        if s == "none" {
            None
        } else {
            Some(s.parse().unwrap())
        }
    };
    let options = thermopack_rs::VolumeFlashOptions {
        temperature_guess: guess(c[3]),
        pressure_guess: guess(c[4]),
        minimum_temperature: c[5].parse().unwrap(),
        maximum_temperature: c[6].parse().unwrap(),
        minimum_pressure: c[7].parse().unwrap(),
        maximum_pressure: c[8].parse().unwrap(),
        ..model.volume_flash_options()
    };
    let target: f64 = c[2].parse().unwrap();
    let volume: f64 = f[3].parse().unwrap();
    let flash = match c[1] {
        "UV" => model.uv_flash_with_options(target, volume, options),
        "SV" => model.sv_flash_with_options(target, volume, options),
        _ => panic!("unknown fixture kind"),
    }
    .unwrap_or_else(|e| panic!("{case}: {e}"));
    assert!(
        (options.minimum_temperature..=options.maximum_temperature).contains(&flash.temperature)
    );
    assert!((options.minimum_pressure..=options.maximum_pressure).contains(&flash.pressure));
    let q = model.equilibrium_properties(&flash).unwrap();
    assert!((q.volume / volume - 1.).abs() < 2e-8, "{case}: specified V");
    let (actual, scale) = if c[1] == "UV" {
        (q.internal_energy, model.gas_constant() * flash.temperature)
    } else {
        (q.entropy, model.gas_constant())
    };
    assert!(
        (actual - target).abs() <= 3e-10 * scale,
        "{case}: specified caloric target"
    );
    match flash.phase {
        thermopack_rs::FlashPhase::Single => {
            assert_eq!((flash.vapor_fraction, flash.liquid_fraction), (-1., -1.));
        }
        thermopack_rs::FlashPhase::Liquid => {
            assert_eq!((flash.vapor_fraction, flash.liquid_fraction), (0., 1.));
        }
        thermopack_rs::FlashPhase::Vapor => {
            assert_eq!((flash.vapor_fraction, flash.liquid_fraction), (1., 0.));
        }
        thermopack_rs::FlashPhase::TwoPhase => {
            assert!((0. ..=1.).contains(&flash.vapor_fraction));
            assert!((0. ..=1.).contains(&flash.liquid_fraction));
            assert!((flash.vapor_fraction + flash.liquid_fraction - 1.).abs() <= 1e-12);
        }
    }
    (model, flash)
}

fn source_close(case: &str, index: usize, actual: f64, expected: f64, reduced: bool) {
    let floors = if reduced {
        [2e-8, 1e-10, 0., 3e-6, 3e-6, 1e-8, 1e-8, 1e-8, 0.]
    } else {
        [2e-6, 1e-3, 0., 3e-6, 3e-6, 1e-12, 1e-4, 1e-6, 0.]
    };
    let relative = if index == 2 || index == 8 { 0. } else { 3e-7 };
    assert!(
        actual.is_finite()
            && expected.is_finite()
            && (actual - expected).abs()
                <= floors[index] + relative * actual.abs().max(expected.abs()),
        "{case} field {index}: {actual} != {expected}"
    );
}

fn replay_source(fixture: &str, count: usize, declared_phase_gaps: bool) {
    let mut rows = 0;
    for row in fixture.lines().skip(1) {
        let f: Vec<_> = row.split(';').collect();
        assert_eq!(f.len(), 10);
        let (model, flash) = fixture_flash(f[0]);
        let q = model.equilibrium_properties(&flash).unwrap();
        let two = flash.phase == thermopack_rs::FlashPhase::TwoPhase;
        let actual = [
            flash.temperature,
            flash.pressure,
            f64::from(two),
            if two { flash.vapor_fraction } else { 0. },
            if two { flash.liquid_fraction } else { 1. },
            q.volume,
            q.internal_energy,
            q.entropy,
            flash.phase as i32 as f64,
        ];
        let expected: Vec<f64> = f[1..].iter().map(|s| s.parse().unwrap()).collect();
        let case: Vec<_> = f[0].split('|').collect();
        let seed_t: f64 = case[2].parse().unwrap();
        let seed_v: f64 = case[3].parse().unwrap();
        let endpoint = seed_t < model.tc()
            && model.saturation_pressure(seed_t).is_ok_and(|sat| {
                (seed_v / sat.liquid_volume - 1.).abs() < 2e-6
                    || (seed_v / sat.vapor_volume - 1.).abs() < 2e-6
            });
        if declared_phase_gaps {
            assert!(endpoint, "{}: phase gap outside endpoint grid", f[0]);
        }
        // Endpoint flag fidelity is audited in the source executable's own
        // arithmetic environment. Replay always checks physical coordinates;
        // a different floating-point build may choose the opposite side.
        for i in [0, 1, 5, 6, 7] {
            source_close(f[0], i, actual[i], expected[i], model.gas_constant() == 1.);
        }
        if !endpoint {
            for i in [2, 3, 4, 8] {
                source_close(f[0], i, actual[i], expected[i], model.gas_constant() == 1.);
            }
        }
        rows += 1;
    }
    assert_eq!(rows, count);
}

#[test]
fn original_source_uv_sv_coordinates_replay_without_fortran() {
    replay_source(
        include_str!("fixtures/multiparameter_isochoric.csv"),
        181,
        false,
    );
}

#[test]
fn declared_source_endpoint_flags_retain_coordinate_and_specification_checks() {
    replay_source(
        include_str!("fixtures/multiparameter_isochoric_phase_gaps.csv"),
        59,
        true,
    );
}

#[test]
fn failed_source_calls_have_original_kernel_qualified_native_recoveries() {
    let mut rows = 0;
    for row in include_str!("fixtures/multiparameter_isochoric_recovered.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = row.split(';').collect();
        assert_eq!(f.len(), 6);
        let (model, flash) = fixture_flash(f[0]);
        let q = model.equilibrium_properties(&flash).unwrap();
        for ((index, actual), expected) in [
            (0, flash.temperature),
            (1, flash.pressure),
            (5, q.volume),
            (6, q.internal_energy),
            (7, q.entropy),
        ]
        .into_iter()
        .zip(&f[1..])
        {
            source_close(
                f[0],
                index,
                actual,
                expected.parse().unwrap(),
                model.gas_constant() == 1.,
            );
        }
        rows += 1;
    }
    assert_eq!(rows, 320);
}
