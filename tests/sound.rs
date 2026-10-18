mod support;
use thermopack_rs::{Cubic, CubicKind, CubicOptions, FlashPhase, FlashResult, Phase};

fn pure_split(t: f64, p: f64, z: &[f64], quality: f64) -> FlashResult {
    FlashResult {
        temperature: t,
        pressure: p,
        feed: z.to_vec(),
        phase: FlashPhase::TwoPhase,
        liquid: z.to_vec(),
        vapor: z.to_vec(),
        vapor_fraction: quality,
        liquid_fraction: 1. - quality,
        liquid_volume: None,
        vapor_volume: None,
        iterations: 0,
        fugacity_residual: 0.,
    }
}

#[test]
fn tv_and_equilibrium_sound_match_the_exact_source() {
    for (fixture, expected_count, tv) in [
        (include_str!("fixtures/tv_sound.csv"), 1258, true),
        (include_str!("fixtures/equilibrium_sound.csv"), 817, false),
    ] {
        let mut count = 0;
        for row in fixture.lines().skip(1) {
            let columns: Vec<_> = row.split(';').collect();
            let fields: Vec<_> = columns[0].split('|').collect();
            let model = support::fixture_model(&fields);
            let t = fields[2].parse().unwrap();
            let value = fields[3].parse().unwrap();
            let z: Vec<f64> = fields[4].split(',').map(|s| s.parse().unwrap()).collect();
            let actual = if tv {
                let amount: f64 = fields[5]
                    .split('@')
                    .find(|op| op.starts_with("amount:"))
                    .map(|op| op.split(':').nth(1).unwrap().parse().unwrap())
                    .unwrap_or(1.);
                let moles: Vec<_> = z.iter().map(|z| z * amount).collect();
                vec![
                    model
                        .baseline_speed_of_sound_tv(t, value * amount, &moles)
                        .unwrap(),
                ]
            } else {
                let config: Vec<_> = fields[5].split('@').collect();
                let flash =
                    if let Some(quality) = config.iter().find(|op| op.starts_with("quality:")) {
                        pure_split(
                            t,
                            value,
                            &z,
                            quality.split(':').nth(1).unwrap().parse().unwrap(),
                        )
                    } else {
                        model.baseline_tp_flash(t, value, &z).unwrap()
                    };
                let mut roots = [Phase::Liquid, Phase::Vapor];
                if let Some(data) = config.iter().find(|op| op.starts_with("roots:")) {
                    for (root, label) in roots.iter_mut().zip(data.split(':').skip(1)) {
                        *root = if label == "L" {
                            Phase::Liquid
                        } else {
                            Phase::Vapor
                        };
                    }
                }
                vec![
                    flash.phase as i32 as f64,
                    model
                        .baseline_equilibrium_speed_of_sound_with_phases(&flash, roots)
                        .unwrap(),
                ]
            };
            assert_eq!(actual.len(), columns.len() - 1);
            for (i, (a, b)) in actual.iter().zip(&columns[1..]).enumerate() {
                let b: f64 = b.parse().unwrap();
                let tolerance = if !tv && i == 0 {
                    0.
                } else {
                    1e-5 + 2e-7 * a.abs().max(b.abs())
                };
                assert!(
                    (a - b).abs() <= tolerance,
                    "{} field {i}: {a}, {b}",
                    columns[0]
                );
            }
            count += 1;
        }
        assert_eq!(count, expected_count);
    }
}

fn isentropic_check(model: &Cubic, flash: &FlashResult) {
    let properties = model.equilibrium_properties(flash).unwrap();
    let speed = model.equilibrium_speed_of_sound(flash).unwrap();
    let p = flash.pressure;
    let dp = p * 1e-4;
    let low = model
        .ps_flash(p - dp, &flash.feed, properties.entropy)
        .unwrap();
    let high = model
        .ps_flash(p + dp, &flash.feed, properties.entropy)
        .unwrap();
    let vl = model.equilibrium_properties(&low).unwrap().volume;
    let vh = model.equilibrium_properties(&high).unwrap().volume;
    let mass: f64 = flash
        .feed
        .iter()
        .zip(model.components())
        .map(|(z, c)| z * c.molar_mass * 1e-3)
        .sum();
    let difference = (-properties.volume.powi(2) * 2. * dp / ((vh - vl) * mass)).sqrt();
    assert!(
        (speed - difference).abs() < 1e-3 + 3e-4 * speed,
        "{:?} T={} P={p}: {speed}, {difference}",
        model.kind(),
        flash.temperature
    );
    if let Some(response) = model.equilibrium_response(flash).unwrap() {
        let predicted = (-properties.volume.powi(2)
            / (mass * response.isentropic_volume_dp(flash.temperature)))
        .sqrt();
        assert!((speed - predicted).abs() < 1e-7 + 1e-9 * speed);
    }
}

#[test]
fn physical_sound_follows_constant_entropy_in_mixtures_and_pure_coexistence() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        let model = Cubic::with_options(
            kind,
            &["C1", "C2", "C3"],
            CubicOptions {
                volume_translation: kind != CubicKind::VanDerWaals,
                ..Default::default()
            },
        )
        .unwrap();
        for t in [180., 250., 350.] {
            isentropic_check(&model, &model.tp_flash(t, 1e6, &[0.5, 0.3, 0.2]).unwrap());
        }
        let pure = Cubic::new(kind, &["C1"]).unwrap();
        let sat = pure.saturation_pressure(150.).unwrap();
        for quality in [0.01, 0.5, 0.99] {
            isentropic_check(&pure, &pure_split(150., sat.pressure, &[1.], quality));
        }
    }
}

#[test]
fn tv_stiffness_does_not_require_positive_pressure() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1"]).unwrap();
    let state = model.prepare(100., &[1.]).unwrap();
    let volume = 1e-4;
    assert!(state.pressure(volume).unwrap().value < 0.);
    assert_eq!(state.baseline_speed_of_sound(volume).unwrap(), 0.);
    assert!(state.speed_of_sound(volume).is_err());
    assert!(state.residual(volume).is_err());
    let stretched = 5e-5;
    assert!(state.pressure(stretched).unwrap().value < 0.);
    let physical = state.speed_of_sound(stretched).unwrap();
    assert!(physical > 0. && physical.is_finite());
    assert!((physical - state.baseline_speed_of_sound(stretched).unwrap()).abs() < 1e-7);
    assert!(state.baseline_speed_of_sound(state.covolume()).is_err());
    let sparse = Cubic::new(CubicKind::PengRobinson, &["C1", "C2", "N2"]).unwrap();
    let p = model.saturation_pressure(150.).unwrap().pressure;
    let a = model
        .equilibrium_speed_of_sound(&pure_split(150., p, &[1.], 0.5))
        .unwrap();
    let b = sparse
        .equilibrium_speed_of_sound(&pure_split(150., p, &[1., 0., 0.], 0.5))
        .unwrap();
    assert!((a - b).abs() < 1e-7);
    let mut bad = pure_split(150., p, &[1.], 0.5);
    bad.vapor_fraction = f64::NAN;
    assert!(model.equilibrium_speed_of_sound(&bad).is_err());
}

#[test]
fn tv_total_volume_and_mole_numbers_scale_together() {
    let model = Cubic::quantum(&["HE", "H2", "NE"]).unwrap();
    let z = [0.2, 0.4, 0.4];
    let state = model.prepare(80., &z).unwrap();
    let v = state.volume(1e6, Phase::Vapor).unwrap();
    let physical = state.speed_of_sound(v).unwrap();
    let baseline = state.baseline_speed_of_sound(v).unwrap();
    for amount in [0.01, 0.3, 3.5, 100.] {
        let moles: Vec<_> = z.iter().map(|z| z * amount).collect();
        let a = model.speed_of_sound_tv(80., v * amount, &moles).unwrap();
        let b = model
            .baseline_speed_of_sound_tv(80., v * amount, &moles)
            .unwrap();
        assert!((a - physical).abs() < 1e-7);
        assert!((b - baseline).abs() < 1e-7);
    }
    for moles in [[0., 0., 0.], [-0.1, 0.5, 0.6], [f64::NAN, 0.5, 0.5]] {
        assert!(model.speed_of_sound_tv(80., v, &moles).is_err());
    }
    assert!(model.speed_of_sound_tv(80., v, &[1.]).is_err());
    assert!(model.speed_of_sound_tv(80., 0., &z).is_err());
}
