use thermopack_rs::{Cubic, CubicKind, CubicOptions, MixingRule};

fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a} vs {b}"
    );
}

#[test]
fn unrecorded_hv_pair_matches_exact_source() {
    // Exact-source PR/HV1 methane/water, no interGE record. cbSelect requests
    // "Classic" interactions; the fork only stores "vdW", so kij is zero.
    let eos = Cubic::with_options(
        CubicKind::PengRobinson,
        &["C1", "H2O"],
        CubicOptions {
            mixing_rule: MixingRule::HuronVidal1,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(eos.kij(0, 1).unwrap(), 0.);
    for (t, v, p, s) in [
        (250., 0.0002, 134917.599618027, 143.27005927420456),
        (250., 0.002, 913422.8257895749, 167.22880433634913),
        (400., 0.0002, 9665699.502878487, 157.28031574053063),
        (400., 0.002, 1573815.2122684529, 180.07198593949337),
    ] {
        let c = eos.prepare(t, &[0.5, 0.5]).unwrap().caloric_tv(v).unwrap();
        close(c.pressure, p, 1e-5, 1e-12);
        close(c.entropy, s, 1e-10, 1e-12);
    }
}

#[test]
fn explicit_fallback_interaction_and_temperature_derivatives_match_classical_rule() {
    for mixing_rule in [MixingRule::HuronVidal1, MixingRule::HuronVidal2] {
        let mut hv = Cubic::with_options(
            CubicKind::PengRobinson,
            &["C1", "H2O"],
            CubicOptions {
                mixing_rule,
                ..Default::default()
            },
        )
        .unwrap();
        let mut classic = Cubic::new(CubicKind::PengRobinson, &["C1", "H2O"]).unwrap();
        for kij in [0., 0.2, -0.05] {
            hv.set_kij(0, 1, kij).unwrap();
            classic.set_kij(0, 1, kij).unwrap();
            for t in [250., 400.] {
                let a = hv.prepare(t, &[0.3, 0.7]).unwrap();
                let b = classic.prepare(t, &[0.3, 0.7]).unwrap();
                let ca = a.caloric_tv(2e-4).unwrap();
                let cb = b.caloric_tv(2e-4).unwrap();
                let pa = a.pressure(2e-4).unwrap();
                let pb = b.pressure(2e-4).unwrap();
                for (x, y) in [
                    (ca.pressure, cb.pressure),
                    (ca.entropy, cb.entropy),
                    (ca.heat_capacity_v, cb.heat_capacity_v),
                    (pa.dt, pb.dt),
                    (pa.dv, pb.dv),
                ] {
                    close(x, y, 1e-7, 2e-12);
                }
            }
        }
    }
}
