use thermopack_rs::{
    EntropyReference, IdealGas, Mbwr, MbwrFluid, MbwrKind, R, component, mbwr_records,
};

const T: f64 = 300.0;
const V: f64 = 7.0e-4;
const N: f64 = 0.7;

#[test]
fn dilute_fugacity_retains_the_first_virial_term_below_pressure_roundoff() {
    let fluid = MbwrFluid::from_record("mbwrcx1").unwrap();
    let b = fluid.mbwr().virial_coefficients(300.).unwrap().0;
    assert!(b.abs() > 1e-8);
    for rho in [1e-12, 1e-20] {
        let q = fluid.caloric_tv(300., 1. / rho, 1.).unwrap();
        let expected = b * rho;
        let actual = q.ln_fugacity_coefficient.unwrap();
        assert!(
            (actual - expected).abs() <= 1e-10 * expected.abs(),
            "rho={rho}: {actual} != {expected}"
        );
    }
}

fn close(actual: f64, expected: f64, atol: f64, rtol: f64, context: &str) {
    assert!(
        actual.is_finite()
            && expected.is_finite()
            && (actual - expected).abs() <= atol + rtol * actual.abs().max(expected.abs()),
        "{context}: {actual} != {expected} (atol={atol}, rtol={rtol})"
    );
}

fn quantities(
    p: thermopack_rs::mbwr_caloric::MbwrTvProperties,
) -> [thermopack_rs::mbwr_caloric::MbwrTvQuantity; 6] {
    [
        p.internal_energy,
        p.enthalpy,
        p.entropy,
        p.helmholtz_energy,
        p.gibbs_energy,
        p.chemical_potential,
    ]
}

fn check_fd(actual: f64, fd: f64, context: &str) {
    // Central differences of extensive energies involve subtracting close
    // values; use a small unit-aware absolute floor plus a relative bound.
    close(actual, fd, 2.0e-5, 3.0e-5, context);
}

#[test]
fn database_backed_caloric_tv_derivatives_and_thermodynamic_identities() {
    let mut supported = 0;
    for record in mbwr_records() {
        let Ok(model) = MbwrFluid::from_record(record.source_name) else {
            continue;
        };
        assert_eq!(
            model.mbwr().record().unwrap().source_name,
            record.source_name
        );
        supported += 1;

        let p = model.caloric_tv(T, V, N).unwrap();
        let h = 2.0e-5;
        let tp = model.caloric_tv(T + T * h, V, N).unwrap();
        let tm = model.caloric_tv(T - T * h, V, N).unwrap();
        let vp = model.caloric_tv(T, V * (1.0 + h), N).unwrap();
        let vm = model.caloric_tv(T, V * (1.0 - h), N).unwrap();
        let np = model.caloric_tv(T, V, N * (1.0 + h)).unwrap();
        let nm = model.caloric_tv(T, V, N * (1.0 - h)).unwrap();
        let q = quantities(p);
        let qt_p = quantities(tp);
        let qt_m = quantities(tm);
        let qv_p = quantities(vp);
        let qv_m = quantities(vm);
        let qn_p = quantities(np);
        let qn_m = quantities(nm);
        let names = ["U", "H", "S", "A", "G", "mu"];
        for i in 0..6 {
            check_fd(
                q[i].dt,
                (qt_p[i].value - qt_m[i].value) / (2.0 * T * h),
                &format!("{} {}_T", record.source_name, names[i]),
            );
            check_fd(
                q[i].dv,
                (qv_p[i].value - qv_m[i].value) / (2.0 * V * h),
                &format!("{} {}_V", record.source_name, names[i]),
            );
            check_fd(
                q[i].dn,
                (qn_p[i].value - qn_m[i].value) / (2.0 * N * h),
                &format!("{} {}_n", record.source_name, names[i]),
            );
        }
        close(
            p.pressure.dt,
            (tp.pressure.value - tm.pressure.value) / (2.0 * T * h),
            1.0,
            3.0e-5,
            "P_T",
        );
        close(
            p.pressure.dv,
            (vp.pressure.value - vm.pressure.value) / (2.0 * V * h),
            5.0e4,
            3.0e-5,
            "P_V",
        );
        close(
            p.pressure.dn,
            (np.pressure.value - nm.pressure.value) / (2.0 * N * h),
            1.0,
            3.0e-5,
            "P_n",
        );

        close(
            p.internal_energy.value,
            p.helmholtz_energy.value + T * p.entropy.value,
            2.0e-7,
            2.0e-12,
            "U=A+TS",
        );
        close(
            p.enthalpy.value,
            p.internal_energy.value + p.pressure.value * V,
            2.0e-7,
            2.0e-12,
            "H=U+PV",
        );
        close(
            p.gibbs_energy.value,
            p.helmholtz_energy.value + p.pressure.value * V,
            2.0e-7,
            2.0e-12,
            "G=A+PV",
        );
        close(
            p.gibbs_energy.value,
            N * p.chemical_potential.value,
            2.0e-7,
            2.0e-11,
            "G=n mu",
        );
        close(
            p.helmholtz_energy.dt,
            -p.entropy.value,
            2.0e-7,
            2.0e-11,
            "A_T=-S",
        );
        close(
            p.helmholtz_energy.dv,
            -p.pressure.value,
            2.0e-7,
            2.0e-11,
            "A_V=-P",
        );
        close(
            p.helmholtz_energy.dn,
            p.chemical_potential.value,
            2.0e-7,
            2.0e-11,
            "A_n=mu",
        );
        close(
            p.chemical_potential.dt,
            -p.entropy.dn,
            2.0e-7,
            2.0e-10,
            "mu_T=-S_n",
        );
        close(
            p.chemical_potential.dv,
            -p.pressure.dn,
            2.0e-7,
            2.0e-10,
            "mu_V=-P_n",
        );
        close(
            p.heat_capacity_v,
            p.internal_energy.dt,
            1.0e-12,
            1.0e-12,
            "Cv=U_T",
        );
        if let Some(cp) = p.heat_capacity_p {
            close(
                cp,
                p.heat_capacity_v - T * p.pressure.dt.powi(2) / p.pressure.dv,
                1.0e-10,
                1.0e-12,
                "Cp response identity",
            );
        }
    }
    assert_eq!(
        supported, 48,
        "number of database records with default ideal metadata"
    );
}

#[test]
fn total_reduced_helmholtz_first_and_second_derivatives_match_local_differences() {
    let model = MbwrFluid::from_record("mbwrcx1").unwrap();
    let f = model.reduced_helmholtz(T, V, N).unwrap();
    let h = 2.0e-5;
    let tlo = model.reduced_helmholtz(T * (1.0 - h), V, N).unwrap();
    let thi = model.reduced_helmholtz(T * (1.0 + h), V, N).unwrap();
    let vlo = model.reduced_helmholtz(T, V * (1.0 - h), N).unwrap();
    let vhi = model.reduced_helmholtz(T, V * (1.0 + h), N).unwrap();
    let nlo = model.reduced_helmholtz(T, V, N * (1.0 - h)).unwrap();
    let nhi = model.reduced_helmholtz(T, V, N * (1.0 + h)).unwrap();
    close(
        f.dt,
        (thi.value - tlo.value) / (2.0 * T * h),
        1.0e-10,
        3.0e-5,
        "F_T",
    );
    close(
        f.dv,
        (vhi.value - vlo.value) / (2.0 * V * h),
        1.0e-5,
        3.0e-5,
        "F_V",
    );
    close(
        f.dn,
        (nhi.value - nlo.value) / (2.0 * N * h),
        1.0e-10,
        3.0e-5,
        "F_n",
    );
    for (actual, plus, minus, step, floor, name) in [
        (f.dtt, thi.dt, tlo.dt, T * h, 1e-12, "F_TT"),
        (f.dtv, vhi.dt, vlo.dt, V * h, 1e-7, "F_TV"),
        (f.dtn, nhi.dt, nlo.dt, N * h, 1e-11, "F_Tn"),
        (f.dvv, vhi.dv, vlo.dv, V * h, 1e-3, "F_VV"),
        (f.dvn, vhi.dn, vlo.dn, V * h, 1e-6, "F_Vn"),
        (f.dnn, nhi.dn, nlo.dn, N * h, 1e-10, "F_nn"),
    ] {
        close(actual, (plus - minus) / (2.0 * step), floor, 5.0e-5, name);
    }
}

#[test]
fn caloric_values_are_extensive_and_intensive_properties_are_invariant() {
    let model = MbwrFluid::from_record("mbwrcx1").unwrap();
    let scale = 3.5;
    let a = model.caloric_tv(T, V, N).unwrap();
    let b = model.caloric_tv(T, V * scale, N * scale).unwrap();
    for ((x, y), name) in quantities(a)
        .into_iter()
        .zip(quantities(b))
        .zip(["U", "H", "S", "A", "G", "mu"])
    {
        let extensive_a = name != "mu";
        close(
            y.value,
            if extensive_a {
                scale * x.value
            } else {
                x.value
            },
            2.0e-7,
            2.0e-11,
            name,
        );
        close(
            y.dt,
            if extensive_a { scale * x.dt } else { x.dt },
            2.0e-7,
            2.0e-10,
            "scaled d/dT",
        );
        close(
            y.dv,
            if extensive_a { x.dv } else { x.dv / scale },
            2.0e-7,
            2.0e-10,
            "scaled d/dV",
        );
        close(
            y.dn,
            if extensive_a { x.dn } else { x.dn / scale },
            2.0e-7,
            2.0e-10,
            "scaled d/dn",
        );
    }
    close(
        b.pressure.value,
        a.pressure.value,
        1.0e-7,
        2.0e-12,
        "P scaling",
    );
    close(b.pressure.dt, a.pressure.dt, 1.0e-7, 2.0e-12, "P_T scaling");
    close(
        b.pressure.dv,
        a.pressure.dv / scale,
        1.0e-7,
        2.0e-12,
        "P_V scaling",
    );
    close(
        b.pressure.dn,
        a.pressure.dn / scale,
        1.0e-7,
        2.0e-12,
        "P_n scaling",
    );
    close(
        b.heat_capacity_v,
        scale * a.heat_capacity_v,
        1.0e-7,
        2.0e-12,
        "Cv scaling",
    );
    if let (Some(x), Some(y)) = (a.heat_capacity_p, b.heat_capacity_p) {
        close(y, scale * x, 1.0e-7, 2.0e-12, "Cp scaling");
    }
    assert_eq!(a.speed_of_sound.is_some(), b.speed_of_sound.is_some());
    if let (Some(x), Some(y)) = (a.speed_of_sound, b.speed_of_sound) {
        close(y, x, 1.0e-10, 2.0e-12, "sound speed scaling");
    }
}

#[test]
fn raw_ideal_references_clone_ownership_and_setter_semantics() {
    let mut model = MbwrFluid::from_record("mbwrcx1").unwrap();
    let original = model.clone();
    let ideal = model.ideal_gas().properties(298.15, 100000.0).unwrap();
    close(
        ideal.enthalpy,
        10025.726002805888,
        2.0e-8,
        2.0e-12,
        "raw H reference",
    );
    close(
        ideal.entropy,
        112.14282031044416,
        2.0e-10,
        2.0e-12,
        "raw S reference",
    );

    let before = model.caloric_tv(T, V, N).unwrap();
    let cp = *model.ideal_gas().heat_capacity_correlation().1;
    model
        .ideal_gas_mut()
        .set_formation_enthalpy(1234.0)
        .unwrap();
    model
        .ideal_gas_mut()
        .set_standard_entropy(150.0, EntropyReference::Bar)
        .unwrap();
    let referenced = model.ideal_gas().properties(298.15, 100000.0).unwrap();
    close(referenced.enthalpy, 1234.0, 1.0e-9, 0.0, "set H reference");
    close(referenced.entropy, 150.0, 1.0e-9, 0.0, "set S reference");

    model.ideal_gas_mut().set_heat_capacity(8, &[2.5]).unwrap();
    assert_eq!(model.ideal_gas().heat_capacity_correlation().1[0], 2.5);
    // Cp edits retain the integration constants, so values move by the
    // difference between the new and old primitive correlations.
    let retained = model.ideal_gas().properties(298.15, 100000.0).unwrap();
    close(
        retained.enthalpy,
        1234.0 - 10025.726002805888 + 2.5 * R * 298.15,
        1.0e-8,
        2.0e-12,
        "Cp edit retains H offset",
    );
    close(
        retained.entropy,
        150.0 - 112.14282031044416 + 2.5 * R * 298.15_f64.ln() - R * 100000.0_f64.ln(),
        1.0e-9,
        2.0e-12,
        "Cp edit retains S offset",
    );
    assert_eq!(original.ideal_gas().heat_capacity_correlation().1, &cp);
    assert_ne!(
        model.caloric_tv(T, V, N).unwrap().helmholtz_energy.value,
        before.helmholtz_energy.value
    );
}

#[test]
fn explicit_ideal_metadata_can_evaluate_every_source_fit() {
    let ideal = IdealGas::new(component("C1").unwrap()).unwrap();
    let mut count = 0;
    for record in mbwr_records() {
        let eos = Mbwr::from_record(record.source_name).unwrap();
        let model = MbwrFluid::with_ideal_gas(eos, ideal.clone(), 0.016).unwrap();
        let state = model.caloric_tv(T, N, N).unwrap();
        assert!(state.pressure.value.is_finite(), "{}", record.source_name);
        assert!(
            state.chemical_potential.value.is_finite(),
            "{}",
            record.source_name
        );
        count += 1;
    }
    assert_eq!(count, 65);
}

#[test]
fn explicit_model_validation_and_type_ten_ideal_derivative_convention() {
    let gas = IdealGas::with_integration_constants(component("N-H2").unwrap(), 0.0, 0.0).unwrap();
    let custom = Mbwr::with_parameters(MbwrKind::Bender19, &[0.0; 19], 1.0).unwrap();
    let model = MbwrFluid::with_ideal_gas(custom, gas, 0.014).unwrap();
    let t = 300.0;
    let v = 1.0e-3;
    let n = 1.0;
    let state = model.caloric_tv(t, v, n).unwrap();
    let baseline = model.baseline_caloric_tv(t, v, n).unwrap();
    assert_eq!(model.ideal_gas().heat_capacity_correlation().0, 10);
    assert_eq!(
        state.helmholtz_energy.value,
        baseline.helmholtz_energy.value
    );
    assert_eq!(state.internal_energy.value, baseline.internal_energy.value);
    assert_ne!(state.heat_capacity_v, baseline.heat_capacity_v);
    let dt = t * 1.0e-5;
    let u_plus = model
        .caloric_tv(t + dt, v, n)
        .unwrap()
        .internal_energy
        .value;
    let u_minus = model
        .caloric_tv(t - dt, v, n)
        .unwrap()
        .internal_energy
        .value;
    close(
        state.heat_capacity_v,
        (u_plus - u_minus) / (2.0 * dt),
        1.0e-5,
        2.0e-6,
        "type-10 physical Cv",
    );
    assert!(
        IdealGas::with_integration_constants(component("N-H2").unwrap(), f64::NAN, 0.0).is_err()
    );
    assert!(
        IdealGas::with_integration_constants(component("N-H2").unwrap(), 0.0, f64::INFINITY)
            .is_err()
    );
    let valid_ideal = IdealGas::new(component("N-H2").unwrap()).unwrap();
    let custom = Mbwr::with_parameters(MbwrKind::Bender19, &[0.0; 19], 1.0).unwrap();
    assert!(MbwrFluid::with_ideal_gas(custom, valid_ideal, 0.0).is_err());
    let valid_ideal = IdealGas::new(component("N-H2").unwrap()).unwrap();
    let custom = Mbwr::with_parameters(MbwrKind::Bender19, &[0.0; 19], 1.0).unwrap();
    assert!(MbwrFluid::with_ideal_gas(custom, valid_ideal, f64::NAN).is_err());
    assert!(MbwrFluid::new(MbwrKind::Mbwr32, "no-such-fluid").is_err());
    assert!(MbwrFluid::from_record("no-such-record").is_err());
    assert!(model.prepare(0.0).is_err());
    assert!(model.prepare(f64::NAN).is_err());
    assert!(model.caloric_tv(t, 0.0, n).is_err());
    assert!(model.caloric_tv(t, v, -1.0).is_err());
    assert!(model.caloric_tv(t, f64::INFINITY, n).is_err());
    assert!(model.caloric_tv(t, v, f64::NAN).is_err());
}
