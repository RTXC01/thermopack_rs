use thermopack_rs::{Mbwr, MbwrKind, Phase, R};

fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a} != {b}"
    );
}

#[test]
fn stored_critical_and_triple_points_are_exact() {
    for (kind, id, expected_tc, expected_pc) in [
        (MbwrKind::Mbwr32, "C1", 190.5554, 4598817.681727696),
        (MbwrKind::Mbwr32, "C2", 305.3407, 4871417.856003743),
        (MbwrKind::Mbwr32, "C3", 369.8508, 4247648.246647852),
        (MbwrKind::Mbwr32, "O2", 154.5792, 5042983.469714031),
        (MbwrKind::Mbwr32, "N2", 126.1987, 3400008.9046928287),
        (MbwrKind::Mbwr32, "R134a", 374.183, 4055998.3116944283),
        (MbwrKind::Bender19, "C1", 190.7699, 4626986.189794588),
        (MbwrKind::Bender19, "H2O", 647.8465, 22258653.904915676),
    ] {
        let model = Mbwr::new(kind, id).unwrap();
        let cp = model.source_critical_point().unwrap();
        close(cp.temperature, expected_tc, 0.01, 1e-4);
        close(cp.pressure, expected_pc, 1.0, 1e-4);
        assert!(cp.volume > 0. && cp.volume < 0.001);
        assert!(cp.compressibility > 0.2 && cp.compressibility < 0.35);

        let tp = model.source_triple_point().unwrap();
        assert!(tp.0 > 0. && tp.0 < cp.temperature);
        assert!(tp.1 > 0. && tp.1 < cp.pressure);
    }
}

#[test]
fn analytic_critical_points_satisfy_inflection_and_positive_pressure() {
    for (kind, id) in [
        (MbwrKind::Mbwr32, "C1"),
        (MbwrKind::Mbwr32, "C2"),
        (MbwrKind::Mbwr32, "C3"),
        (MbwrKind::Mbwr32, "O2"),
        (MbwrKind::Mbwr32, "N2"),
        (MbwrKind::Mbwr32, "R134a"),
        (MbwrKind::Bender19, "C1"),
        (MbwrKind::Bender19, "C2"),
        (MbwrKind::Bender19, "C3"),
        (MbwrKind::Bender19, "O2"),
        (MbwrKind::Bender19, "N2"),
    ] {
        let model = Mbwr::new(kind, id).unwrap();
        let cp = model.critical_point().unwrap();
        assert!(cp.temperature > 50. && cp.temperature < 800.);
        assert!(cp.pressure > 1e5 && cp.pressure < 3e7);
        assert!(cp.volume > 0. && cp.volume < 0.001);
        assert!(cp.compressibility > 0.15 && cp.compressibility < 0.4);

        if let Some(source_cp) = model.source_critical_point() {
            close(cp.temperature, source_cp.temperature, 0.5, 0.01);
            close(cp.pressure, source_cp.pressure, 5e4, 0.02);
            close(cp.volume, source_cp.volume, 1e-5, 0.05);
        }

        let rho = 1. / cp.volume;
        let p = model.pressure_density(cp.temperature, rho).unwrap();
        close(p.value, cp.pressure, 0.01, 1e-7);
        assert!(
            p.drho.abs() < 1e-6 * (R * cp.temperature),
            "{:?} {id} dP/drho = {}",
            model.kind(),
            p.drho
        );
        assert!(
            p.drhorho.abs() < 1e-6 * (R * cp.temperature / rho),
            "{:?} {id} d2P/drho2 = {}",
            model.kind(),
            p.drhorho
        );
    }
}

#[test]
fn branch_isolation_and_volume_roots_close_pressure() {
    let model = Mbwr::new(MbwrKind::Mbwr32, "C1").unwrap();
    let subcritical_t = 150.0;
    let max_rho = 35000.0;
    let branches = model.prepare_equilibrium(subcritical_t, max_rho).unwrap();
    assert_eq!(branches.maximum_density(), max_rho);
    assert_eq!(branches.spinodal_densities().len(), 2);

    let rho_sp1 = branches.spinodal_densities()[0];
    let rho_sp2 = branches.spinodal_densities()[1];
    assert!(rho_sp1 < rho_sp2);
    let p_sp1 = model.pressure_density(subcritical_t, rho_sp1).unwrap();
    let p_sp2 = model.pressure_density(subcritical_t, rho_sp2).unwrap();
    assert!(p_sp1.drho.abs() < 1e-6 * (R * subcritical_t));
    assert!(p_sp2.drho.abs() < 1e-6 * (R * subcritical_t));
    assert!(p_sp1.value > p_sp2.value);

    for p in [1e5, 5e5, 1e6] {
        let roots = branches.densities(p).unwrap();
        assert_eq!(roots.len(), 2);
        assert!(roots[0] < rho_sp1);
        assert!(roots[1] > rho_sp2);

        let v_vap = branches.volume(p, Phase::Vapor).unwrap();
        let v_liq = branches.volume(p, Phase::Liquid).unwrap();
        assert!(v_vap > v_liq);

        let p_vap = model.pressure_tv(subcritical_t, v_vap, 1.0).unwrap();
        let p_liq = model.pressure_tv(subcritical_t, v_liq, 1.0).unwrap();
        close(p_vap.value, p, 0.01, 1e-7);
        close(p_liq.value, p, 0.01, 1e-7);

        let v_min = branches.volume(p, Phase::MinimumGibbs).unwrap();
        assert!(v_min == v_vap || v_min == v_liq);
    }

    let supercritical_t = 250.0;
    let super_branches = model.prepare_equilibrium(supercritical_t, max_rho).unwrap();
    assert!(super_branches.spinodal_densities().is_empty());
    for p in [1e5, 1e6, 1e7] {
        let roots = super_branches.densities(p).unwrap();
        assert_eq!(roots.len(), 1);
        let v = super_branches.volume(p, Phase::Vapor).unwrap();
        let p_calc = model.pressure_tv(supercritical_t, v, 1.0).unwrap();
        close(p_calc.value, p, 0.01, 1e-7);
    }
}

#[test]
fn invalid_inputs_and_unphysical_roots_are_rejected() {
    let model = Mbwr::new(MbwrKind::Mbwr32, "C1").unwrap();
    assert!(model.prepare_equilibrium(0.0, 30000.0).is_err());
    assert!(model.prepare_equilibrium(-100.0, 30000.0).is_err());
    assert!(model.prepare_equilibrium(150.0, 0.0).is_err());
    assert!(model.prepare_equilibrium(150.0, -100.0).is_err());

    let branches = model.prepare_equilibrium(150.0, 30000.0).unwrap();
    assert!(branches.densities(0.0).is_err());
    assert!(branches.densities(-1e5).is_err());
    assert!(branches.densities(f64::NAN).is_err());

    assert!(model.critical_point_from(0.0, 10000.0).is_err());
    assert!(model.critical_point_from(-10.0, 10000.0).is_err());
    assert!(model.critical_point_from(190.0, 0.0).is_err());
    assert!(model.critical_point_from(190.0, -100.0).is_err());
    assert!(model.critical_point_from(f64::NAN, 1000.0).is_err());
}
