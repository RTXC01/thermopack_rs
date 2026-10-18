mod support;
use support::fixture_model;
use thermopack_rs::{Cubic, CubicKind, Phase, RachfordRiceIteration};

#[test]
fn substitution_matches_source_values_and_all_flag_combinations() {
    let mut counts = [0; 4];
    for row in include_str!("fixtures/rr_iteration.csv").lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let model = fixture_model(&fields);
        let z: Vec<f64> = fields[4].split(',').map(|v| v.parse().unwrap()).collect();
        let parameters: Vec<_> = fields[5]
            .split('@')
            .find(|s| s.starts_with("rr:"))
            .unwrap()
            .split(':')
            .collect();
        let k: Vec<f64> = parameters[1]
            .split(',')
            .map(|v| v.parse().unwrap())
            .collect();
        let mut state = RachfordRiceIteration::new(z.len());
        state.split.vapor_fraction = parameters[2].parse().unwrap();
        model
            .rr_successive_substitution_iteration(
                fields[2].parse().unwrap(),
                fields[3].parse().unwrap(),
                &z,
                &k,
                &mut state,
                parameters[3] == "1",
            )
            .unwrap_or_else(|e| panic!("{}: {e}", columns[0]));
        counts[usize::from(state.converged) * 2 + usize::from(state.has_solution)] += 1;
        let mut actual = vec![f64::from(state.converged), f64::from(state.has_solution)];
        if state.has_solution {
            actual.extend([
                state.split.vapor_fraction,
                state.split.liquid_fraction,
                state.gibbs,
            ]);
            actual.extend(
                state
                    .split
                    .liquid
                    .iter()
                    .chain(&state.split.vapor)
                    .chain(&state.equilibrium_ratios)
                    .chain(&state.liquid_ln_fugacity)
                    .chain(&state.vapor_ln_fugacity),
            );
        }
        assert_eq!(actual.len(), columns.len() - 1, "{}", columns[0]);
        for (i, (a, b)) in actual.iter().zip(&columns[1..]).enumerate() {
            let b: f64 = b.parse().unwrap();
            let tolerance = if i < 2 {
                0.
            } else {
                1e-9 + 2e-8 * a.abs().max(b.abs())
            };
            assert!(
                (a - b).abs() <= tolerance,
                "{} field {i}: {a}, {b}",
                columns[0]
            );
        }
    }
    assert_eq!(counts, [1312, 792, 2, 4]);
}

#[test]
fn repeated_updates_close_material_and_chemical_potentials() {
    let z = [0.5, 0.3, 0.2];
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        let model = Cubic::new(kind, &["C1", "C2", "C3"]).unwrap();
        let mut k = model.wilson_k(180., 1e6).unwrap();
        let mut state = RachfordRiceIteration::new(3);
        for _ in 0..100 {
            model
                .rr_successive_substitution_iteration(180., 1e6, &z, &k, &mut state, false)
                .unwrap();
            assert!(state.has_solution);
            k.copy_from_slice(&state.equilibrium_ratios);
            if state.converged {
                break;
            }
        }
        assert!(state.converged, "{kind:?}");
        assert!((0. ..=1.).contains(&state.split.vapor_fraction));
        let l = model.prepare(180., &state.split.liquid).unwrap();
        let v = model.prepare(180., &state.split.vapor).unwrap();
        let mut pl = [0.; 3];
        let mut pv = [0.; 3];
        l.ln_fugacity_tp_into(1e6, l.volume(1e6, Phase::Liquid).unwrap(), &mut pl)
            .unwrap();
        v.ln_fugacity_tp_into(1e6, v.volume(1e6, Phase::Vapor).unwrap(), &mut pv)
            .unwrap();
        for (i, &zi) in z.iter().enumerate() {
            let xi = state.split.liquid[i];
            let yi = state.split.vapor[i];
            assert!(
                (zi - state.split.liquid_fraction * xi - state.split.vapor_fraction * yi).abs()
                    < 1e-12
            );
            assert!((xi.ln() + pl[i] - yi.ln() - pv[i]).abs() < 1e-9);
        }
        let old_phi = state.liquid_ln_fugacity.clone();
        model
            .rr_successive_substitution_iteration(180., 1e6, &z, &[1.; 3], &mut state, false)
            .unwrap();
        assert!(!state.has_solution && !state.converged);
        assert_eq!(state.liquid_ln_fugacity, old_phi);
        state.equilibrium_ratios.pop();
        assert!(
            model
                .rr_successive_substitution_iteration(180., 1e6, &z, &k, &mut state, false)
                .is_err()
        );
        assert!(!state.has_solution && !state.converged);
    }
}

#[test]
fn resizing_public_outputs_does_not_resize_the_owned_iteration_workspace() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2", "C3"]).unwrap();
    let mut state = RachfordRiceIteration::new(2);
    for buffer in [
        &mut state.split.liquid,
        &mut state.split.vapor,
        &mut state.equilibrium_ratios,
        &mut state.liquid_ln_fugacity,
        &mut state.vapor_ln_fugacity,
    ] {
        buffer.resize(3, 0.);
    }
    assert!(
        model
            .rr_successive_substitution_iteration(
                180.,
                1e6,
                &[0.5, 0.3, 0.2],
                &[4., 0.2, 0.2],
                &mut state,
                false
            )
            .is_err()
    );
}
