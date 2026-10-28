use thermopack_rs::{EntropyReference, Mbwr, MbwrFluid, MbwrKind, R};

fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a} != {b} (atol={atol}, rtol={rtol})"
    );
}

fn check_coexistence(model: &Mbwr, fluid: &MbwrFluid, t: f64, max_density: f64) {
    let sat = model.saturation_pressure(t, max_density).unwrap();
    assert!(sat.pressure > 0. && sat.dpdt > 0.);
    assert!(sat.liquid_volume > 0. && sat.vapor_volume > sat.liquid_volume);

    let liquid = model.pressure_tv(t, sat.liquid_volume, 1.).unwrap();
    let vapor = model.pressure_tv(t, sat.vapor_volume, 1.).unwrap();
    close(liquid.value, sat.pressure, 0.02, 2e-9);
    close(vapor.value, sat.pressure, 0.02, 2e-9);
    // P_V<0 is the local mechanical stability condition on each outer branch.
    assert!(liquid.dv < 0. && vapor.dv < 0.);

    let branches = model.prepare_equilibrium(t, max_density).unwrap();
    let spinodals = branches.spinodal_densities();
    assert!(spinodals.len() >= 2);
    let vapor_density = sat.vapor_volume.recip();
    let liquid_density = sat.liquid_volume.recip();
    assert!(spinodals.windows(2).any(|pair| {
        let barrier_density = 0.5 * (pair[0] + pair[1]);
        barrier_density > vapor_density
            && barrier_density < liquid_density
            && model.pressure_density(t, barrier_density).unwrap().drho < 0.
    }));

    let l = fluid.caloric_tv(t, sat.liquid_volume, 1.).unwrap();
    let v = fluid.caloric_tv(t, sat.vapor_volume, 1.).unwrap();
    close(l.pressure.value, sat.pressure, 0.02, 2e-9);
    close(v.pressure.value, sat.pressure, 0.02, 2e-9);
    close(
        l.chemical_potential.value,
        v.chemical_potential.value,
        2e-5,
        2e-9,
    );
    close(l.gibbs_energy.value, v.gibbs_energy.value, 2e-5, 2e-9);
    let stable_roots = branches.densities(sat.pressure).unwrap();
    assert!(
        stable_roots
            .iter()
            .any(|rho| (rho - vapor_density).abs() < 1e-5 * vapor_density)
    );
    assert!(
        stable_roots
            .iter()
            .any(|rho| (rho - liquid_density).abs() < 1e-5 * liquid_density)
    );
    for rho in stable_roots {
        let state = fluid.caloric_tv(t, 1. / rho, 1.).unwrap();
        assert!(state.gibbs_energy.value >= l.gibbs_energy.value - 2e-5);
    }

    let latent_h = v.enthalpy.value - l.enthalpy.value;
    let latent_s = v.entropy.value - l.entropy.value;
    assert!(latent_h > 0. && latent_s > 0.);
    close(latent_h, t * latent_s, 2e-5, 2e-8);
    close(
        sat.dpdt,
        latent_h / (t * (sat.vapor_volume - sat.liquid_volume)),
        2e-5,
        2e-8,
    );
}

#[test]
fn ordinary_states_close_pressure_chemical_potential_and_clapeyron() {
    for (id, t) in [("C1", 150.), ("C3", 290.), ("R134a", 300.)] {
        let model = Mbwr::new(MbwrKind::Mbwr32, id).unwrap();
        let fluid = MbwrFluid::new(MbwrKind::Mbwr32, id).unwrap();
        let critical = model.critical_point().unwrap();
        let max_density = 4. * critical.volume.recip();
        check_coexistence(&model, &fluid, t, max_density);
    }
}

#[test]
fn clapeyron_tangent_matches_independent_curve_difference() {
    let model = Mbwr::new(MbwrKind::Mbwr32, "C1").unwrap();
    let tc = model.critical_point().unwrap().temperature;
    let max_density = 4. * model.critical_point().unwrap().volume.recip();
    for tr in [0.7, 0.9] {
        let t = tr * tc;
        let dt = 2e-5 * t;
        let s = model.saturation_pressure(t, max_density).unwrap();
        let hi = model.saturation_pressure(t + dt, max_density).unwrap();
        let lo = model.saturation_pressure(t - dt, max_density).unwrap();
        close(s.dpdt, (hi.pressure - lo.pressure) / (2. * dt), 2e-3, 2e-5);
    }
}

#[test]
fn prepared_and_convenience_equilibrium_solvers_agree() {
    let model = Mbwr::new(MbwrKind::Mbwr32, "C1").unwrap();
    let critical = model.critical_point().unwrap();
    let max_density = 4. * critical.volume.recip();
    let prepared = model
        .prepare_equilibrium(150., max_density)
        .unwrap()
        .saturation_pressure()
        .unwrap();
    let convenience = model.saturation_pressure(150., max_density).unwrap();
    close(prepared.temperature, convenience.temperature, 0., 0.);
    close(prepared.pressure, convenience.pressure, 0., 0.);
    close(prepared.liquid_volume, convenience.liquid_volume, 0., 0.);
    close(prepared.vapor_volume, convenience.vapor_volume, 0., 0.);
    close(prepared.dpdt, convenience.dpdt, 0., 0.);
}

#[test]
fn inverse_temperature_is_confined_to_explicit_bracket_and_roundtrips() {
    for (kind, id) in [
        (MbwrKind::Mbwr32, "C1"),
        (MbwrKind::Mbwr32, "C3"),
        (MbwrKind::Mbwr32, "R134a"),
    ] {
        let model = Mbwr::new(kind, id).unwrap();
        let critical = model.critical_point().unwrap();
        let max_density = 4. * critical.volume.recip();
        let t = 0.82 * critical.temperature;
        let forward = model.saturation_pressure(t, max_density).unwrap();
        let inverse = model
            .saturation_temperature(
                forward.pressure,
                0.75 * critical.temperature,
                0.9 * critical.temperature,
                max_density,
            )
            .unwrap();
        close(inverse.temperature, t, 2e-7, 2e-9);
        close(inverse.pressure, forward.pressure, 1e-4, 2e-10);
        assert!(
            model
                .saturation_temperature(
                    forward.pressure,
                    0.5 * critical.temperature,
                    0.7 * critical.temperature,
                    max_density
                )
                .is_err()
        );
    }
}

#[test]
fn near_critical_states_remain_distinct_and_supercritical_states_fail() {
    for (kind, id) in [
        (MbwrKind::Mbwr32, "C1"),
        (MbwrKind::Mbwr32, "C3"),
        (MbwrKind::Mbwr32, "R134a"),
    ] {
        let model = Mbwr::new(kind, id).unwrap();
        let critical = model.critical_point().unwrap();
        let max_density = 2. * critical.volume.recip();
        for delta in [1e-3, 1e-4, 1e-5, 1e-6] {
            let t = critical.temperature * (1. - delta);
            let sat = model
                .saturation_pressure(t, max_density)
                .unwrap_or_else(|e| {
                    panic!("{id}, 1-Tr={delta}: expected two physical outer branches: {e}")
                });
            assert!(sat.vapor_volume > sat.liquid_volume);
            assert!((sat.vapor_volume - sat.liquid_volume) / sat.vapor_volume > 1e-7);
            let p_l = model.pressure_tv(t, sat.liquid_volume, 1.).unwrap();
            let p_v = model.pressure_tv(t, sat.vapor_volume, 1.).unwrap();
            close(p_l.value, sat.pressure, 0.05, 1e-8);
            close(p_v.value, sat.pressure, 0.05, 1e-8);
            assert!(p_l.dv < 0. && p_v.dv < 0.);
        }
        assert!(
            model
                .saturation_pressure(critical.temperature * (1. + 1e-5), max_density)
                .is_err()
        );
    }
}

#[test]
fn total_reference_shifts_cancel_between_coexisting_phases() {
    let model = Mbwr::new(MbwrKind::Mbwr32, "C1").unwrap();
    let mut fluid = MbwrFluid::new(MbwrKind::Mbwr32, "C1").unwrap();
    let sat = model.saturation_pressure(150., 35000.).unwrap();
    let before_l = fluid.caloric_tv(150., sat.liquid_volume, 1.).unwrap();
    let before_v = fluid.caloric_tv(150., sat.vapor_volume, 1.).unwrap();
    fluid
        .ideal_gas_mut()
        .set_formation_enthalpy(12345.)
        .unwrap();
    fluid
        .ideal_gas_mut()
        .set_standard_entropy(85., EntropyReference::Bar)
        .unwrap();
    let after_l = fluid.caloric_tv(150., sat.liquid_volume, 1.).unwrap();
    let after_v = fluid.caloric_tv(150., sat.vapor_volume, 1.).unwrap();

    // Reference changes shift both phases equally while leaving coexistence
    // residuals, pressure, and the EOS-only saturation answer untouched.
    close(
        after_l.chemical_potential.value - before_l.chemical_potential.value,
        after_v.chemical_potential.value - before_v.chemical_potential.value,
        2e-7,
        2e-9,
    );
    close(
        after_l.gibbs_energy.value - before_l.gibbs_energy.value,
        after_v.gibbs_energy.value - before_v.gibbs_energy.value,
        2e-7,
        2e-9,
    );
    close(after_l.pressure.value, before_l.pressure.value, 1e-8, 1e-12);
    close(after_v.pressure.value, before_v.pressure.value, 1e-8, 1e-12);
    let sat_after = fluid.saturation_pressure(150., 35000.).unwrap();
    close(sat_after.pressure, sat.pressure, 1e-8, 1e-12);
    close(sat_after.dpdt, sat.dpdt, 1e-8, 1e-12);
}

#[test]
fn explicit_domains_and_pressure_windows_are_enforced() {
    let model = Mbwr::new(MbwrKind::Mbwr32, "C1").unwrap();
    let good = model.saturation_pressure(150., 35000.).unwrap();
    assert!(model.saturation_pressure(150., 1000.).is_err());
    assert!(
        model
            .saturation_pressure_bounded(150., 1.1 * good.pressure, 2. * good.pressure, 35000.)
            .is_err()
    );
    assert!(
        model
            .saturation_pressure_bounded(150., 0.5 * good.pressure, 0.9 * good.pressure, 35000.)
            .is_err()
    );
    for (lo, hi) in [
        (0., good.pressure),
        (good.pressure, good.pressure),
        (2. * good.pressure, good.pressure),
        (f64::NAN, 2. * good.pressure),
    ] {
        assert!(
            model
                .saturation_pressure_bounded(150., lo, hi, 35000.)
                .is_err()
        );
    }
    for t in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(model.saturation_pressure(t, 35000.).is_err());
    }
    for p in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(model.saturation_temperature(p, 100., 180., 35000.).is_err());
    }
    assert!(
        model
            .saturation_temperature(good.pressure, 180., 100., 35000.)
            .is_err()
    );
    assert!(
        model
            .saturation_temperature(good.pressure, 0., 180., 35000.)
            .is_err()
    );
    assert!(
        model
            .saturation_temperature(good.pressure, 100., 180., 1000.)
            .is_err()
    );
}

#[test]
fn nitrogen_cold_coexistence_survives_pressure_roundoff_near_density_bound() {
    let model = Mbwr::from_record("MBWRcx9").unwrap();
    let t = 82.02848069145;
    let max_density = 44844.05784;
    let sat = model.saturation_pressure(t, max_density).unwrap();
    assert!(sat.pressure > 0.);
    assert!(sat.vapor_volume > sat.liquid_volume);
    let liquid = model.pressure_tv(t, sat.liquid_volume, 1.).unwrap();
    let vapor = model.pressure_tv(t, sat.vapor_volume, 1.).unwrap();
    close(liquid.value, sat.pressure, 0.02, 2e-8);
    close(vapor.value, sat.pressure, 0.02, 2e-8);
    assert!(liquid.dv < 0. && vapor.dv < 0.);
    assert!((sat.vapor_volume - sat.liquid_volume) / sat.vapor_volume > 1e-3);
    assert!(sat.liquid_volume.recip() < max_density);
}

#[test]
fn independently_tabulated_critical_fits_have_two_resolved_near_critical_spinodals() {
    // The 49 fits in this fixture have independent source critical roots.
    // Test the local pressure-derivative structure instead of assuming every
    // extrapolated fit has a physically meaningful full coexistence curve.
    let mut checked = 0;
    for row in include_str!("fixtures/mbwr_critical.csv").lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let model = Mbwr::from_record(fields[1]).unwrap();
        let t_critical: f64 = columns[1].parse().unwrap();
        let rho_critical: f64 = columns[3].parse().unwrap();
        let t = t_critical * (1. - 1e-6);
        let eq = model.prepare_equilibrium(t, 2. * rho_critical).unwrap();
        let spinodals = eq.spinodal_densities();
        let near: Vec<_> = spinodals
            .iter()
            .copied()
            .filter(|rho| (rho / rho_critical - 1.).abs() < 0.1)
            .collect();
        assert_eq!(near.len(), 2, "{} spinodals {spinodals:?}", fields[1]);
        for rho in near {
            let p = model.pressure_density(t, rho).unwrap();
            close(p.drho, 0., 2e-4 * R * t, 0.);
        }
        checked += 1;
    }
    assert_eq!(checked, 49);
}
