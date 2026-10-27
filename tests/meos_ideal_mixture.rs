use thermopack_rs::{MultiparameterIdealMixture, Phase};

fn close(actual: f64, expected: f64, tolerance: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale.max(expected.abs()).max(actual.abs()),
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance}, scale={scale}"
    );
}

fn check_derivatives(components: &[&str], x: &[f64], t: f64, p: f64, phase: Phase) {
    let model = MultiparameterIdealMixture::new(components).unwrap();
    let state = model.residual_tp(t, p, x, phase).unwrap();
    let ht = t * 1e-4;
    let hp = p * 1e-4;
    let tp = model.residual_tp(t + ht, p, x, phase).unwrap();
    let tm = model.residual_tp(t - ht, p, x, phase).unwrap();
    let pp = model.residual_tp(t, p + hp, x, phase).unwrap();
    let pm = model.residual_tp(t, p - hp, x, phase).unwrap();

    let fd = |a: f64, b: f64, h: f64| (a - b) / (2. * h);
    close(
        state.compressibility_factor_dt,
        fd(tp.compressibility_factor, tm.compressibility_factor, ht),
        2e-5,
        1. / t,
    );
    close(
        state.compressibility_factor_dp,
        fd(pp.compressibility_factor, pm.compressibility_factor, hp),
        2e-5,
        1. / p,
    );
    close(
        state.residual_enthalpy_dt,
        fd(tp.residual_enthalpy, tm.residual_enthalpy, ht),
        2e-5,
        1.,
    );
    close(
        state.residual_enthalpy_dp,
        fd(pp.residual_enthalpy, pm.residual_enthalpy, hp),
        2e-5,
        1e-8,
    );
    close(
        state.residual_entropy_dt,
        fd(tp.residual_entropy, tm.residual_entropy, ht),
        2e-5,
        1e-3,
    );
    close(
        state.residual_entropy_dp,
        fd(pp.residual_entropy, pm.residual_entropy, hp),
        2e-5,
        1e-8,
    );
    close(
        state.residual_gibbs_energy_dt,
        fd(tp.residual_gibbs_energy, tm.residual_gibbs_energy, ht),
        2e-5,
        1.,
    );
    close(
        state.residual_gibbs_energy_dp,
        fd(pp.residual_gibbs_energy, pm.residual_gibbs_energy, hp),
        2e-5,
        1e-8,
    );
    close(
        state.residual_gibbs_energy_dt,
        -state.residual_entropy,
        2e-12,
        1e-10,
    );
    close(
        state.residual_enthalpy,
        state.residual_gibbs_energy + t * state.residual_entropy,
        2e-12,
        1e-8,
    );
    let v_mix: f64 = x
        .iter()
        .zip(&state.component_volumes)
        .map(|(xi, vi)| xi * vi)
        .sum();
    close(
        state.residual_gibbs_energy_dp,
        v_mix - thermopack_rs::R * t / p,
        2e-12,
        1e-10,
    );
    for (i, lnphi) in state.ln_fugacity_coefficients.iter().enumerate() {
        close(
            state.ln_fugacity_coefficients_dt[i],
            -state.residual_enthalpy_composition_coefficients[i] / (thermopack_rs::R * t * t),
            2e-12,
            1e-12,
        );
        assert!(lnphi.is_finite());
        assert!(
            state.ln_fugacity_coefficients_composition_derivatives[i]
                .iter()
                .all(|v| *v == 0.)
        );
    }
}

#[test]
fn tp_derivatives_and_residual_identities_for_binary_and_ternary() {
    for phase in [Phase::Liquid, Phase::Vapor] {
        check_derivatives(&["C3", "R134A"], &[0.35, 0.65], 450., 1e6, phase);
        check_derivatives(&["N-H2", "O-H2", "P-H2"], &[0.2, 0.3, 0.5], 45., 1e6, phase);
    }
}

#[test]
fn composition_changes_follow_linear_coefficients() {
    let model = MultiparameterIdealMixture::new(&["C3", "R134A", "N-H2"]).unwrap();
    let x = [0.2, 0.35, 0.45];
    let step = 1e-5;
    for (i, j) in [(0, 1), (0, 2), (1, 2)] {
        let base = model.residual_tp(450., 1e6, &x, Phase::Vapor).unwrap();
        let mut plus = x;
        let mut minus = x;
        plus[i] += step;
        plus[j] -= step;
        minus[i] -= step;
        minus[j] += step;
        let a = model.residual_tp(450., 1e6, &plus, Phase::Vapor).unwrap();
        let b = model.residual_tp(450., 1e6, &minus, Phase::Vapor).unwrap();
        let check = |fd: f64, coefficients: &[f64]| {
            close(
                fd,
                coefficients[i] - coefficients[j],
                2e-5,
                coefficients[i].abs().max(coefficients[j].abs()).max(1e-8),
            );
        };
        check(
            (a.compressibility_factor - b.compressibility_factor) / (2. * step),
            &base.compressibility_factor_composition_coefficients,
        );
        check(
            (a.residual_enthalpy - b.residual_enthalpy) / (2. * step),
            &base.residual_enthalpy_composition_coefficients,
        );
        check(
            (a.residual_entropy - b.residual_entropy) / (2. * step),
            &base.residual_entropy_composition_coefficients,
        );
        check(
            (a.residual_gibbs_energy - b.residual_gibbs_energy) / (2. * step),
            &base.residual_gibbs_composition_coefficients,
        );
    }
}

#[test]
fn zero_weight_pure_limit_repeated_components_and_clone() {
    let mix = MultiparameterIdealMixture::new(&["C3", "R134A"]).unwrap();
    let pure_mix = mix.residual_tp(450., 1e5, &[1., 0.], Phase::Vapor).unwrap();
    let pure = MultiparameterIdealMixture::new(&["C3"])
        .unwrap()
        .residual_tp(450., 1e5, &[1.], Phase::Vapor)
        .unwrap();
    close(
        pure_mix.compressibility_factor,
        pure.compressibility_factor,
        1e-13,
        1.,
    );
    close(
        pure_mix.residual_enthalpy,
        pure.residual_enthalpy,
        1e-13,
        1.,
    );
    close(pure_mix.residual_entropy, pure.residual_entropy, 1e-13, 1.);
    close(
        pure_mix.residual_gibbs_energy,
        pure.residual_gibbs_energy,
        1e-13,
        1.,
    );
    assert!(pure_mix.component_volumes[1].is_finite());

    let repeated = MultiparameterIdealMixture::new(&["C3", "C3"]).unwrap();
    let repeated_state = repeated
        .residual_tp(450., 1e5, &[0.4, 0.6], Phase::Vapor)
        .unwrap();
    let clone = repeated.clone();
    let clone_state = clone
        .residual_tp(450., 1e5, &[0.4, 0.6], Phase::Vapor)
        .unwrap();
    close(
        repeated_state.residual_gibbs_energy,
        clone_state.residual_gibbs_energy,
        0.,
        1.,
    );
    close(
        repeated_state.compressibility_factor,
        repeated_state.component_compressibility_factors[0],
        1e-13,
        1.,
    );
}

#[test]
fn invalid_construction_and_state_inputs_are_rejected() {
    assert!(MultiparameterIdealMixture::new(&[]).is_err());
    assert!(MultiparameterIdealMixture::new(&["not-a-fluid"]).is_err());
    let model = MultiparameterIdealMixture::new(&["C3", "R134A"]).unwrap();
    for t in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(
            model
                .residual_tp(t, 1e5, &[0.5, 0.5], Phase::Vapor)
                .is_err()
        );
    }
    for p in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(
            model
                .residual_tp(450., p, &[0.5, 0.5], Phase::Vapor)
                .is_err()
        );
    }
    for x in [
        vec![-0.1, 1.1],
        vec![f64::NAN, 1.],
        vec![f64::INFINITY, 0.],
        vec![0.2, 0.2],
        vec![0.5],
    ] {
        assert!(model.residual_tp(450., 1e5, &x, Phase::Vapor).is_err());
    }
    assert!(
        model
            .residual_tp(450., 1e5, &[0.5, 0.5], Phase::MinimumGibbs)
            .is_err()
    );
}

#[test]
fn baseline_gibbs_pressure_derivative_retains_source_bug() {
    let state = MultiparameterIdealMixture::new(&["C3", "R134A"])
        .unwrap()
        .residual_tp(450., 1e6, &[0.35, 0.65], Phase::Vapor)
        .unwrap();
    assert_eq!(state.baseline_gibbs_dp, state.residual_gibbs_energy);
    assert!(
        (state.baseline_gibbs_dp - state.residual_gibbs_energy_dp).abs()
            > 1e-8 * state.residual_gibbs_energy.abs().max(1.)
    );
}
