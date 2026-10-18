mod support;
use support::fixture_model;
use thermopack_rs::{Cubic, CubicKind, CubicOptions, FlashPhase, Phase};

fn compare_fixture(text: &str) -> usize {
    let mut count = 0;
    for row in text.lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let model = if fields[5] == "TP" {
            let kind = match fields[0] {
                "PR" => CubicKind::PengRobinson,
                "SRK" => CubicKind::SoaveRedlichKwong,
                "RK" => CubicKind::RedlichKwong,
                "VDW" => CubicKind::VanDerWaals,
                _ => panic!("kind"),
            };
            Cubic::new(kind, &fields[1].split(',').collect::<Vec<_>>()).unwrap()
        } else {
            fixture_model(&fields)
        };
        let z: Vec<f64> = fields[4].split(',').map(|s| s.parse().unwrap()).collect();
        let flash = model
            .baseline_tp_flash(fields[2].parse().unwrap(), fields[3].parse().unwrap(), &z)
            .unwrap_or_else(|e| panic!("{}: {e}", columns[0]));
        let mut actual = vec![
            flash.phase as i32 as f64,
            flash.vapor_fraction,
            flash.liquid_fraction,
        ];
        actual.extend(flash.liquid.iter().chain(&flash.vapor));
        assert_eq!(actual.len(), columns.len() - 2);
        for (i, (a, b)) in actual.iter().zip(&columns[2..]).enumerate() {
            let b: f64 = b.parse().unwrap();
            let tolerance = if i == 0 {
                0.
            } else {
                2e-7 + 2e-7 * a.abs().max(b.abs())
            };
            assert!(
                (a - b).abs() <= tolerance,
                "{} field {i}: {a}, {b}",
                columns[0]
            );
        }
        count += 1;
    }
    count
}

#[test]
fn source_iteration_order_matches_raw_classical_and_translated_splits() {
    assert_eq!(compare_fixture(include_str!("fixtures/tp_flash.csv")), 1043);
    assert_eq!(
        compare_fixture(include_str!("fixtures/translated_tp.csv")),
        1130
    );
}

fn quantum_model(kind: CubicKind, components: &[&str]) -> Cubic {
    Cubic::with_options(
        kind,
        components,
        CubicOptions {
            parameter_reference: "QuantumCubic",
            ..CubicOptions::default()
        },
    )
    .unwrap()
}

#[test]
fn early_ss_return_preserves_the_sources_stationary_metastable_split() {
    for (components, z) in [
        (&["H2", "NE"][..], vec![0.6, 0.4]),
        (&["HE", "H2", "NE"][..], vec![0., 0.5, 0.5]),
    ] {
        let model = quantum_model(CubicKind::PengRobinson, components);
        let flash = model.baseline_tp_flash(3., 1e3, &z).unwrap();
        assert_eq!(flash.phase, FlashPhase::TwoPhase);
        assert!(flash.fugacity_residual < 1e-8);
        assert!((flash.vapor_fraction - z[z.len() - 2]).abs() < 1e-12);
        assert!((flash.liquid[z.len() - 1] - 1.).abs() < 1e-12);
        assert!((flash.vapor[z.len() - 2] - 1.).abs() < 1e-12);
        let g = |x: &[f64], phase: Phase| {
            let state = model.prepare(3., x).unwrap();
            let volume = state.volume(1e3, phase).unwrap();
            let mut phi = vec![0.; x.len()];
            state.ln_fugacity_tp_into(1e3, volume, &mut phi).unwrap();
            x.iter()
                .zip(phi)
                .filter(|(x, _)| **x > 0.)
                .map(|(x, f)| x * (x.ln() + f))
                .sum::<f64>()
        };
        let feed = g(&z, Phase::Liquid).min(g(&z, Phase::Vapor));
        let split = flash.liquid_fraction * g(&flash.liquid, Phase::Liquid)
            + flash.vapor_fraction * g(&flash.vapor, Phase::Vapor);
        assert!(split > feed + 0.5);
    }
}

#[test]
fn gibbs_change_stopping_retains_vanishing_phase_without_claiming_chemical_closure() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
    ] {
        let model = quantum_model(kind, &["C1", "CO2", "N2"]);
        let flash = model.baseline_tp_flash(80., 1e5, &[0.5, 0.3, 0.2]).unwrap();
        assert_eq!(flash.phase, FlashPhase::TwoPhase);
        assert!(flash.vapor_fraction > 0. && flash.vapor_fraction < 1e-20);
        assert!(flash.fugacity_residual > 1.);
        for (x, z) in flash.liquid.iter().zip([0.5, 0.3, 0.2]) {
            assert!((x - z).abs() < 1e-12);
        }
        // Raw ghost-vapor compositions still differ from the source by up to
        // 3e-5. The complete diagnostic retains those three failures.
    }
}
