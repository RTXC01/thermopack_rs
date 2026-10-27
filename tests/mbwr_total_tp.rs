use thermopack_rs::{
    EntropyReference, IdealGas, Mbwr, MbwrFluid, MbwrKind, Phase, R, component, mbwr_records,
};

const T: f64 = 150.0;
const P: f64 = 1.0e6;
const N: f64 = 0.7;
const RHOMAX: f64 = 30_000.0;

fn close(a: f64, b: f64, atol: f64, rtol: f64, what: &str) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{what}: {a} != {b}"
    );
}

fn derivatives(q: thermopack_rs::MbwrTpQuantity) -> thermopack_rs::MbwrTpDerivatives {
    q.derivatives.expect("regular physical TP response")
}

fn quantities(p: thermopack_rs::MbwrTpProperties) -> [thermopack_rs::MbwrTpQuantity; 6] {
    [
        p.internal_energy,
        p.enthalpy,
        p.entropy,
        p.helmholtz_energy,
        p.gibbs_energy,
        p.chemical_potential,
    ]
}

// Solve P(T,V,n)=Ptarget locally from a known root. Reusing the nearby root
// avoids a global phase selector changing branches during finite differences.
fn close_local_root(
    model: &MbwrFluid,
    t: f64,
    p: f64,
    n: f64,
    seed_volume: f64,
) -> thermopack_rs::MbwrTpProperties {
    let mut v = seed_volume;
    for _ in 0..20 {
        let q = model.caloric_tv(t, v, n).unwrap();
        let error = q.pressure.value - p;
        if error.abs() <= 1e-10_f64.max(2e-12 * p.abs()) {
            return model.caloric_tvp(t, v, n).unwrap();
        }
        assert!(q.pressure.dv.is_finite() && q.pressure.dv.abs() > 1e-14);
        let step = error / q.pressure.dv;
        // Local changes in T/P/n should need only small corrections; this
        // safeguard prevents a distant Newton jump if a fit is ill-conditioned.
        v = (v - step.clamp(-0.2 * v, 0.2 * v)).max(1e-16);
    }
    panic!("local pressure root did not converge for T={t}, P={p}, n={n}");
}

#[test]
fn physical_tp_response_matches_closed_pressure_finite_differences() {
    let model = MbwrFluid::from_record("mbwrcx1").unwrap();
    let state = model.caloric_tp(T, P, N, Phase::Vapor, RHOMAX).unwrap();
    assert!(state.volume.value > 0.0);
    close(state.pressure, P, 1e-5, 2e-10, "closed pressure");
    close(
        state.gibbs_energy.value,
        state.enthalpy.value - T * state.entropy.value,
        2e-7,
        2e-11,
        "G=H-TS",
    );
    close(
        state.chemical_potential.value,
        state.gibbs_energy.value / N,
        2e-7,
        2e-11,
        "mu=G/n",
    );

    let prepared = model.prepare_tp(T, RHOMAX).unwrap();
    let from_prepared = prepared.caloric_tp(P, N, Phase::Vapor).unwrap();
    close(
        from_prepared.volume.value,
        state.volume.value,
        0.0,
        0.0,
        "prepared volume",
    );
    close(
        from_prepared.enthalpy.value,
        state.enthalpy.value,
        0.0,
        0.0,
        "prepared H",
    );
    let h = 2e-5;
    let tp = model
        .caloric_tp(T * (1. + h), P, N, Phase::Vapor, RHOMAX)
        .unwrap();
    let tm = model
        .caloric_tp(T * (1. - h), P, N, Phase::Vapor, RHOMAX)
        .unwrap();
    let pp = model
        .caloric_tp(T, P * (1. + h), N, Phase::Vapor, RHOMAX)
        .unwrap();
    let pm = model
        .caloric_tp(T, P * (1. - h), N, Phase::Vapor, RHOMAX)
        .unwrap();
    let np = model
        .caloric_tp(T, P, N * (1. + h), Phase::Vapor, RHOMAX)
        .unwrap();
    let nm = model
        .caloric_tp(T, P, N * (1. - h), Phase::Vapor, RHOMAX)
        .unwrap();
    for (i, name) in ["U", "H", "S", "A", "G", "mu"].into_iter().enumerate() {
        let q = quantities(state)[i];
        let d = derivatives(q);
        let fd_t = (quantities(tp)[i].value - quantities(tm)[i].value) / (2. * T * h);
        let fd_p = (quantities(pp)[i].value - quantities(pm)[i].value) / (2. * P * h);
        let fd_n = (quantities(np)[i].value - quantities(nm)[i].value) / (2. * N * h);
        close(d.dt, fd_t, 3e-5, 5e-5, &format!("{name}_T|P"));
        close(d.dp, fd_p, 3e-8, 8e-5, &format!("{name}_P|T"));
        close(d.dn, fd_n, 3e-5, 8e-5, &format!("{name}_n|T,P"));
    }
    let d = derivatives(state.gibbs_energy);
    close(d.dt, -state.entropy.value, 3e-7, 3e-9, "G_T=-S");
    close(d.dp, state.volume.value, 3e-12, 3e-9, "G_P=V");
    close(d.dn, state.chemical_potential.value, 3e-7, 3e-9, "G_n=mu");
}

#[test]
fn all_default_fits_have_local_tp_enthalpy_and_entropy_responses() {
    let mut supported = 0;
    for record in mbwr_records() {
        let Ok(model) = MbwrFluid::from_record(record.source_name) else {
            continue;
        };
        let t = 300.;
        let p = 1.0e5;
        let n = 0.8;
        let q = model
            .caloric_tp(t, p, n, Phase::Vapor, RHOMAX)
            .unwrap_or_else(|e| panic!("{} TP root: {e}", record.source_name));
        let h = 1.0e-5;
        let tp = close_local_root(&model, t * (1. + h), p, n, q.volume.value);
        let tm = close_local_root(&model, t * (1. - h), p, n, q.volume.value);
        let pp = close_local_root(&model, t, p * (1. + h), n, q.volume.value);
        let pm = close_local_root(&model, t, p * (1. - h), n, q.volume.value);
        for (name, property, plus_t, minus_t, plus_p, minus_p) in [
            (
                "H",
                q.enthalpy,
                tp.enthalpy,
                tm.enthalpy,
                pp.enthalpy,
                pm.enthalpy,
            ),
            (
                "S", q.entropy, tp.entropy, tm.entropy, pp.entropy, pm.entropy,
            ),
        ] {
            let d = derivatives(property);
            close(
                d.dt,
                (plus_t.value - minus_t.value) / (2. * t * h),
                5e-5,
                1e-4,
                &format!("{} {name}_T", record.source_name),
            );
            close(
                d.dp,
                (plus_p.value - minus_p.value) / (2. * p * h),
                3e-7,
                2e-4,
                &format!("{} {name}_P", record.source_name),
            );
        }
        supported += 1;
    }
    assert_eq!(supported, 48);
}

#[test]
fn physical_tp_extensivity_phase_selection_and_validation() {
    let model = MbwrFluid::from_record("mbwrcx1").unwrap();
    let liquid = model.caloric_tp(T, P, N, Phase::Liquid, RHOMAX).unwrap();
    let vapor = model.caloric_tp(T, P, N, Phase::Vapor, RHOMAX).unwrap();
    let min = model
        .caloric_tp(T, P, N, Phase::MinimumGibbs, RHOMAX)
        .unwrap();
    close(
        min.gibbs_energy.value,
        liquid.gibbs_energy.value.min(vapor.gibbs_energy.value),
        1e-7,
        1e-11,
        "minimum-G root",
    );
    assert!(liquid.volume.value <= vapor.volume.value);
    let scaled = model
        .caloric_tp(T, P, 2.5 * N, Phase::Vapor, RHOMAX)
        .unwrap();
    for (i, (x, y)) in quantities(vapor)
        .into_iter()
        .zip(quantities(scaled))
        .enumerate()
    {
        close(
            y.value,
            if i == 5 { x.value } else { 2.5 * x.value },
            5e-7,
            3e-10,
            "extensive scaling",
        );
    }
    close(
        scaled.volume.value,
        2.5 * vapor.volume.value,
        1e-15,
        2e-10,
        "volume scaling",
    );
    close(
        scaled.compressibility.value,
        vapor.compressibility.value,
        1e-12,
        2e-10,
        "Z invariant",
    );
    assert!(model.caloric_tp(T, 0., N, Phase::Vapor, RHOMAX).is_err());
    assert!(model.caloric_tp(T, P, 0., Phase::Vapor, RHOMAX).is_err());
    assert!(model.prepare_tp(0., RHOMAX).is_err());
    assert!(model.prepare_tp(T, 0.).is_err());
}

#[test]
fn classical_source_tp_keeps_imposed_pressure_semantics_and_raw_iteration() {
    let model = MbwrFluid::from_record("mbwrcx1").unwrap();
    let v = N * R * T / P;
    let source = model.baseline_caloric_tp_at_volume(T, P, v, N).unwrap();
    // Compare source H/S to the existing qualified source residual API, then
    // add the same ideal mixture contributions used by the model.
    let residual = model.mbwr().baseline_residual_tp(T, P, v, N).unwrap();
    // The source total routine caches ideal H/S at 1 Pa, then applies the
    // explicit -R ln(P) pressure term itself.
    let ideal = model.ideal_gas().properties(T, 1.0).unwrap();
    close(
        source.enthalpy.value,
        residual.enthalpy.value + N * ideal.enthalpy,
        2e-7,
        2e-11,
        "source H",
    );
    close(
        source.entropy.value,
        residual.entropy.value + N * (ideal.entropy - R * P.ln()),
        2e-7,
        2e-11,
        "source S",
    );
    // Imposed P differs from the EOS pressure at the selected volume. This
    // API intentionally retains classical source H/S behavior at that input.
    let physical = model.caloric_tvp(T, v, N).unwrap();
    assert!((physical.pressure - P).abs() > 1.0);
    assert!((source.enthalpy.value - physical.enthalpy.value).abs() > 1e-6);

    let params = model.mbwr().source_density_parameters().unwrap();
    let prepared = model.prepare_baseline_tp(T, params).unwrap();
    let a = model.baseline_caloric_tp(T, P, N, Phase::Vapor).unwrap();
    let b = prepared.caloric_tp(P, N, Phase::Vapor).unwrap();
    assert_eq!(a.iteration.density, b.iteration.density);
    assert_eq!(a.iteration.path, b.iteration.path);
    let props = a.properties.unwrap();
    assert!(a.iteration.density > 0.0);
    assert!(props.enthalpy.value.is_finite() && props.entropy.value.is_finite());
}

#[test]
fn raw_reference_changes_and_type_ten_have_distinct_source_cp_response() {
    let mut model = MbwrFluid::from_record("mbwrcx1").unwrap();
    let t = 300.;
    let v = 7e-4;
    let n = 0.7;
    let original = model.caloric_tvp(t, v, n).unwrap();
    model.ideal_gas_mut().set_formation_enthalpy(1234.).unwrap();
    model
        .ideal_gas_mut()
        .set_standard_entropy(150., EntropyReference::Bar)
        .unwrap();
    let changed = model.caloric_tvp(t, v, n).unwrap();
    close(
        changed.enthalpy.value - original.enthalpy.value,
        n * (1234. - 10025.726002805888),
        2e-7,
        2e-10,
        "reference H shift",
    );
    let hydrogen = component("N-H2").unwrap();
    let ideal10 = IdealGas::with_integration_constants(hydrogen, 0.0, 0.0).unwrap();
    let mut model10 = MbwrFluid::with_ideal_gas(
        Mbwr::from_record("mbwrcx1").unwrap(),
        ideal10,
        hydrogen.molar_mass * 1.0e-3,
    )
    .unwrap();
    model10
        .ideal_gas_mut()
        .set_heat_capacity(10, &[2.0])
        .unwrap();
    let p = model10.caloric_tvp(t, v, n).unwrap();
    assert!(p.heat_capacity_p.is_some());
    let base = model10.baseline_caloric_tvp(t, v, n).unwrap();
    assert_ne!(p.heat_capacity_v, base.heat_capacity_v);
}

#[test]
fn critical_spinodal_retains_values_and_marks_tp_derivatives_unavailable() {
    // P=RT*rho-a*rho² has an exact spinodal at rho=1000 mol/m³.
    let t = 300.0;
    let rho = 1000.0;
    let mut coefficients = [0.0; 19];
    coefficients[1] = R * t * 1000.0 / (2.0 * rho);
    let eos = Mbwr::with_parameters(MbwrKind::Bender19, &coefficients, 1.0).unwrap();
    let ideal = IdealGas::with_integration_constants(component("AR").unwrap(), 0.0, 0.0).unwrap();
    let model = MbwrFluid::with_ideal_gas(eos, ideal, 0.04).unwrap();
    let q = model.caloric_tvp(t, 1.0 / rho, 1.0).unwrap();
    assert!(q.pressure.is_finite() && q.pressure > 0.0);
    for property in quantities(q).into_iter().take(4) {
        assert!(property.value.is_finite());
        assert!(
            property.derivatives.is_none(),
            "spinodal derivatives should be unavailable"
        );
    }
    assert!(q.volume.value.is_finite() && q.volume.derivatives.is_none());
    assert!(q.compressibility.value.is_finite() && q.compressibility.derivatives.is_none());
    assert!(q.ln_fugacity_coefficient.value.is_finite());
    assert!(q.gibbs_energy.derivatives.is_some());
    assert!(q.chemical_potential.derivatives.is_some());
    assert!(q.ln_fugacity_coefficient.derivatives.is_some());
}

#[test]
fn dilute_fugacity_tp_derivatives_retain_the_virial_response() {
    let model = MbwrFluid::from_record("mbwrcx1").unwrap();
    let t = 300.0;
    let rho = 1.0e-6;
    let n = 1.0;
    let q = model.caloric_tvp(t, n / rho, n).unwrap();
    let phi = q.ln_fugacity_coefficient;
    let d = derivatives(phi);
    let b = model.mbwr().virial_coefficients(t).unwrap().0;
    let dt = 0.02;
    let bp = model.mbwr().virial_coefficients(t + dt).unwrap().0;
    let bm = model.mbwr().virial_coefficients(t - dt).unwrap().0;
    let expected_t = rho * ((bp - bm) / (2.0 * dt) - b / t);
    let expected_p = b / (R * t);
    close(phi.value, b * rho, 1e-18, 2e-5, "dilute lnphi");
    close(d.dt, expected_t, 1e-18, 2e-3, "dilute lnphi_T|P");
    close(d.dp, expected_p, 1e-20, 2e-5, "dilute lnphi_P|T");
    assert!(d.dt.is_finite() && d.dp.is_finite());
    assert_ne!(d.dp, 0.0);
}
