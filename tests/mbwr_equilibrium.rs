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
    assert!(model.prepare_equilibrium(150., f64::MAX).is_err());
}

#[test]
fn custom_ideal_pressure_and_explicit_density_domain_are_respected() {
    let model = Mbwr::with_parameters(MbwrKind::Bender19, &[0.; 19], 0.1).unwrap();
    assert!(model.critical_point().is_err());
    assert!(model.critical_point_from(300., 1000.).is_err());
    assert!(model.source_critical_point().is_none());
    assert!(model.source_triple_point().is_none());
    let branches = model.prepare_equilibrium(300., 1234.).unwrap();
    assert!(branches.spinodal_densities().is_empty());
    for rho in [1e-8, 1., 1234.] {
        let p = rho * R * 300.;
        let roots = branches.densities(p).unwrap();
        assert_eq!(roots.len(), 1);
        close(roots[0], rho, 1e-14, 1e-13);
        for phase in [Phase::Vapor, Phase::Liquid, Phase::MinimumGibbs] {
            close(branches.volume(p, phase).unwrap(), 1. / rho, 1e-14, 1e-13);
        }
    }
    assert!(
        branches
            .densities(1234. * R * 300. * 1.0001)
            .unwrap()
            .is_empty()
    );
    assert!(
        branches
            .volume(1234. * R * 300. * 1.0001, Phase::Vapor)
            .is_err()
    );
}

#[test]
fn source_metadata_and_independent_original_pressure_critical_roots() {
    let rows = include_str!("fixtures/mbwr_metadata.csv");
    let mut count = 0;
    for row in rows.lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let model = Mbwr::from_record(fields[1]).unwrap();
        let expected: Vec<f64> = columns[1..].iter().map(|s| s.parse().unwrap()).collect();
        let mut actual = [0.; 8];
        if let Some(c) = model.source_critical_point() {
            actual[..5].copy_from_slice(&[
                1.,
                c.temperature,
                c.pressure,
                1. / c.volume,
                c.compressibility,
            ]);
        }
        if let Some((t, p)) = model.source_triple_point() {
            actual[5..].copy_from_slice(&[1., t, p]);
        }
        for (&a, b) in actual.iter().zip(expected) {
            close(a, b, 1e-12, 2e-14);
        }
        count += 1;
    }
    assert_eq!(count, 65);
    let mut count = 0;
    for row in include_str!("fixtures/mbwr_critical.csv").lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let model = Mbwr::from_record(fields[1]).unwrap();
        let c = model
            .critical_point_from(fields[2].parse().unwrap(), fields[3].parse().unwrap())
            .unwrap();
        let actual = [c.temperature, c.pressure, 1. / c.volume, c.compressibility];
        let atol = [2e-5, 0.1, 0.01, 1e-7];
        for (i, a) in actual.into_iter().enumerate() {
            close(a, columns[i + 1].parse().unwrap(), atol[i], 2e-7);
        }
        let rho = 1. / c.volume;
        let default = model.critical_point().unwrap();
        close(default.temperature, c.temperature, 2e-5, 2e-7);
        close(default.volume, c.volume, 1e-10, 2e-7);
        let q = model.pressure_density(c.temperature, rho).unwrap();
        assert!(q.drho.abs() < 1e-8 * R * c.temperature);
        assert!((rho * q.drhorho).abs() < 1e-8 * R * c.temperature);
        // Resolve the physical spinodal pair arbitrarily closer than a fixed
        // discovery grid can distinguish; verify both original derivatives.
        let t = c.temperature * (1. - 1e-7);
        let branches = model.prepare_equilibrium(t, 2. * rho).unwrap();
        let near: Vec<_> = branches
            .spinodal_densities()
            .iter()
            .filter(|&&r| r > 0.9 * rho && r < 1.1 * rho)
            .collect();
        assert_eq!(
            near.len(),
            2,
            "{} critical {rho}, spinodals {:?}",
            fields[1],
            branches.spinodal_densities()
        );
        for &&r in &near {
            assert!(model.pressure_density(t, r).unwrap().drho.abs() < 1e-7 * R * t);
        }
        count += 1;
    }
    assert_eq!(count, 49);
}

#[test]
fn original_finite_closed_stable_density_iterates_are_in_bounded_root_sets() {
    // Consecutive source liquid/vapor requests share an isotherm/domain. Cache
    // its prepared native branches; checking roots at several P reuses them.
    let mut previous = String::new();
    let mut model = None;
    let mut branches = None;
    let mut count = 0;
    for row in include_str!("fixtures/mbwr_roots.csv").lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        let f: Vec<_> = columns[0].split('|').collect();
        let key = format!("{}|{}|{}", f[1], f[2], f[4]);
        let t: f64 = f[2].parse().unwrap();
        if key != previous {
            drop(branches.take());
            model = Some(Mbwr::from_record(f[1]).unwrap());
            branches = Some(
                model
                    .as_ref()
                    .unwrap()
                    .prepare_equilibrium(t, f[4].parse().unwrap())
                    .unwrap(),
            );
            previous = key;
        }
        let p: f64 = f[3].parse().unwrap();
        let expected: f64 = columns[1].parse().unwrap();
        let roots = branches.as_ref().unwrap().densities(p).unwrap();
        assert!(
            roots
                .iter()
                .any(|&r| (r - expected).abs() <= 1e-3 + 2e-6 * r.abs().max(expected.abs())),
            "{} expected {expected}, roots {roots:?}",
            columns[0]
        );
        for rho in roots {
            let q = model.as_ref().unwrap().pressure_density(t, rho).unwrap();
            close(q.value, p, 0.01, 2e-9);
            assert!(q.drho >= -1e-8 * R * t);
        }
        count += 1;
    }
    assert_eq!(count, 3406);
}
