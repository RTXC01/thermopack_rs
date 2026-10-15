use thermopack_rs::{Cubic, CubicKind, FlashPhase, Phase, numerics};

#[test]
fn flash_material_and_fugacity_closure() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        let model = Cubic::new(kind, &["C1", "C2", "C3"]).unwrap();
        let z = [0.5, 0.3, 0.2];
        let mut two_phase = 0;
        for t in [150.0, 200.0, 250.0, 300.0, 400.0] {
            for p in [1e5, 1e6, 5e6, 1e7] {
                let f = model
                    .tp_flash(t, p, &z)
                    .unwrap_or_else(|e| panic!("{kind:?} T={t} P={p}: {e}"));
                if f.phase == FlashPhase::TwoPhase {
                    two_phase += 1;
                    assert!(f.vapor_fraction > 0.0 && f.vapor_fraction < 1.0);
                    for (i, zi) in z.iter().enumerate() {
                        assert!(
                            (zi - f.liquid_fraction * f.liquid[i] - f.vapor_fraction * f.vapor[i])
                                .abs()
                                < 1e-10
                        );
                    }
                    assert!((f.liquid.iter().sum::<f64>() - 1.0).abs() < 1e-12);
                    assert!((f.vapor.iter().sum::<f64>() - 1.0).abs() < 1e-12);
                    let l = model.prepare(t, &f.liquid).unwrap();
                    let v = model.prepare(t, &f.vapor).unwrap();
                    let mut pl = [0.0; 3];
                    let mut pv = [0.0; 3];
                    l.ln_fugacity_into(l.volume(p, Phase::Liquid).unwrap(), &mut pl)
                        .unwrap();
                    v.ln_fugacity_into(v.volume(p, Phase::Vapor).unwrap(), &mut pv)
                        .unwrap();
                    for i in 0..3 {
                        let residual = f.liquid[i].ln() + pl[i] - f.vapor[i].ln() - pv[i];
                        assert!(residual.abs() < 2e-9, "{kind:?} {t} {p}: {residual}");
                    }
                }
            }
        }
        assert!(two_phase > 0);
    }
}

#[test]
fn pivoted_linear_solve() {
    let mut a = [0.0, 2.0, 1.0, 1.0, -2.0, -3.0, 2.0, 3.0, 1.0];
    let mut b = [3.0, 0.0, 7.0];
    numerics::solve_linear_in_place(&mut a, &mut b).unwrap();
    for (actual, expected) in b.iter().zip([1.0, 2.0, -1.0]) {
        assert!((actual - expected).abs() < 1e-12, "{b:?}");
    }
}

#[test]
fn zero_and_trace_feed_components() {
    let model = Cubic::new(CubicKind::PengRobinson, &["N2", "C1", "C2", "C3"]).unwrap();
    let z = [0.0, 0.7, 0.2999999999, 1e-10];
    for t in [180.0, 240.0, 400.0] {
        let result = model.tp_flash(t, 1e6, &z).unwrap();
        assert_eq!(result.liquid[0], 0.0);
        assert_eq!(result.vapor[0], 0.0);
        if result.phase == FlashPhase::TwoPhase {
            for (i, zi) in z.iter().enumerate() {
                assert!(
                    (zi - result.liquid_fraction * result.liquid[i]
                        - result.vapor_fraction * result.vapor[i])
                        .abs()
                        < 1e-12
                );
            }
        }
    }
}

#[test]
#[allow(clippy::excessive_precision)] // Keep the upstream regression input verbatim.
fn near_critical_upstream_regression() {
    let model = Cubic::new(
        CubicKind::PengRobinson,
        &["CO2", "N2", "O2", "AR", "CO", "H2S", "H2", "C1"],
    )
    .unwrap();
    let z = [
        0.98199899999065965,
        6.0000000025901721e-3,
        1.0000000003932113e-6,
        5.0000000019329808e-4,
        4.0000000015823992e-4,
        9.9999999995230694e-5,
        1.0000000006101501e-2,
        1.0000000003014362e-3,
    ];
    let result = model
        .tp_flash(302.46137724669268, 7705980.6837845352, &z)
        .unwrap();
    if result.phase == FlashPhase::TwoPhase {
        for (i, zi) in z.iter().enumerate() {
            assert!(
                (zi - result.liquid_fraction * result.liquid[i]
                    - result.vapor_fraction * result.vapor[i])
                    .abs()
                    < 1e-10
            );
        }
        assert!(result.fugacity_residual < 1e-9);
    }
}
