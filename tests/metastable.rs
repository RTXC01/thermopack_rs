use thermopack_rs::metastable::{IsentropeEnd, IsentropeOptions};
use thermopack_rs::{Cubic, CubicKind, CubicOptions, MixingRule, Phase};

fn model(label: &str) -> (Cubic, Vec<f64>) {
    let c: Vec<_> = label.split('|').collect();
    let names: Vec<_> = c[2].split(',').collect();
    let kind = match c[1] {
        "SRK" => CubicKind::SoaveRedlichKwong,
        "RK" => CubicKind::RedlichKwong,
        "VDW" => CubicKind::VanDerWaals,
        _ => CubicKind::PengRobinson,
    };
    let eos = if c[3] == "quantum" {
        Cubic::quantum(&names).unwrap()
    } else {
        Cubic::with_options(
            kind,
            &names,
            CubicOptions {
                volume_translation: c[3] == "translated",
                mixing_rule: match c[3] {
                    "hv1" => MixingRule::HuronVidal1,
                    "unifac" => MixingRule::Unifac,
                    _ => MixingRule::VdW,
                },
                ..Default::default()
            },
        )
        .unwrap()
    };
    (eos, c[5].split(',').map(|s| s.parse().unwrap()).collect())
}

fn close(a: f64, b: f64, atol: f64, rtol: f64, label: &str) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{label}: {a} vs {b}"
    );
}

#[test]
fn local_ps_matches_every_qualified_source_return() {
    let mut count = 0;
    for row in include_str!("fixtures/homogeneous_ps.csv").lines().skip(1) {
        let f: Vec<_> = row.split(';').collect();
        let (eos, z) = model(f[0]);
        let c: Vec<_> = f[0].split('|').collect();
        let q: Vec<f64> = c[4]
            .split(':')
            .skip(1)
            .map(|s| s.parse().unwrap())
            .collect();
        let expected: Vec<f64> = f[1..].iter().map(|s| s.parse().unwrap()).collect();
        let state = eos.homogeneous_ps(q[0], q[1], &z, q[2], q[3]).unwrap();
        close(state.temperature, expected[0], 2e-5, 2e-6, f[0]);
        close(state.volume, expected[1], 1e-10, 2e-6, f[0]);
        close(
            state.pressure_residual,
            0.,
            1.1e-10 * q[0].abs().max(1e5),
            0.,
            f[0],
        );
        close(
            state.entropy_residual,
            0.,
            1.1e-10 * q[1].abs().max(1.),
            0.,
            f[0],
        );
        count += 1;
    }
    assert_eq!(count, 298);
}

#[test]
fn baseline_maps_match_complete_source_curves() {
    let mut label = "";
    let mut points = Vec::new();
    let mut states = 0;
    let mut curves = 0;
    let mut warning_coordinates = 0;
    for row in include_str!("fixtures/metastable_isentrope.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = row.split(';').collect();
        if label != f[0] {
            label = f[0];
            let (eos, z) = model(label);
            let c: Vec<_> = label.split('|').collect();
            let q: Vec<f64> = c[4]
                .split(':')
                .skip(1)
                .map(|s| s.parse().unwrap())
                .collect();
            points = eos
                .baseline_metastable_isentrope(&z, q[0], q[1], q[2], q[3] as usize)
                .unwrap()
                .points;
            assert_eq!(points.len(), 12);
            curves += 1;
        }
        let i: usize = f[1].parse().unwrap();
        let e: Vec<f64> = f[2..].iter().map(|s| s.parse().unwrap()).collect();
        let p = points[i];
        for (j, a) in [p.temperature, p.volume, p.pressure]
            .into_iter()
            .enumerate()
        {
            close(a, e[j], [2e-5, 1e-10, 2.][j], 2e-6, label);
        }
        if i > 0 && i + 1 < points.len() {
            close(
                p.pressure_residual,
                0.,
                1.1e-10 * p.pressure.abs().max(1e5),
                0.,
                label,
            );
            close(
                p.entropy_residual,
                0.,
                1.1e-10 * e[5].abs().max(1.),
                0.,
                label,
            );
        } else if e[3].abs() > 1e-7 * e[2].abs().max(1e5) || e[4].abs() > 1e-7 * e[5].abs().max(1.)
        {
            // Source first and terminal coordinates are reported without a PS
            // correction. Keep the warning visible rather than forcing closure.
            warning_coordinates += 1;
        }
        states += 1;
    }
    assert_eq!((curves, states, warning_coordinates), (97, 1164, 25));
}

#[test]
fn ps_jacobian_identities_hold_for_fallback_mixing_and_translation() {
    for config in ["classic", "translated", "hv1", "unifac"] {
        let names = if config == "hv1" { "C1,H2O" } else { "C1,C3" };
        let label = format!("META|PR|{names}|{config}|PS|0.4,0.6");
        let (eos, z) = model(&label);
        for (t, v) in [(250., 2e-4), (400., 2e-3)] {
            let state = eos.prepare(t, &z).unwrap();
            let p = state.pressure(v).unwrap();
            let caloric = state.caloric_tv(v).unwrap();
            let dt = t * 1e-5;
            let dv = v * 1e-5;
            let dsdt = (eos
                .prepare(t + dt, &z)
                .unwrap()
                .caloric_tv(v)
                .unwrap()
                .entropy
                - eos
                    .prepare(t - dt, &z)
                    .unwrap()
                    .caloric_tv(v)
                    .unwrap()
                    .entropy)
                / (2. * dt);
            let dsdv = (state.caloric_tv(v + dv).unwrap().entropy
                - state.caloric_tv(v - dv).unwrap().entropy)
                / (2. * dv);
            close(dsdt, caloric.heat_capacity_v / t, 1e-8, 2e-7, &label);
            close(dsdv, p.dt, 1e-5, 2e-7, &label);
        }
    }
}

#[test]
fn fallback_hv_pair_matches_source_pressure_and_entropy() {
    let (eos, z) = model("META|PR|C1,H2O|hv1|PS|0.5,0.5");
    let state = eos.prepare(250., &z).unwrap();
    close(
        state.pressure(2e-4).unwrap().value,
        134917.599618027,
        1e-5,
        1e-12,
        "fallback P",
    );
    close(
        state.caloric_tv(2e-4).unwrap().entropy,
        143.27005927420456,
        1e-10,
        1e-12,
        "fallback S",
    );
}

#[test]
fn physical_isentrope_closes_ps_and_refines_stability_endpoint() {
    let eos = Cubic::new(CubicKind::PengRobinson, &["C1"]).unwrap();
    let z = [1.];
    let state = eos.prepare(180., &z).unwrap();
    let v = state.volume(6e6, Phase::Liquid).unwrap();
    let s = state.caloric_tv(v).unwrap().entropy;
    let seed = eos.homogeneous_ps(6e6, s, &z, 180., v).unwrap();
    for count in [2, 50] {
        let curve = eos
            .homogeneous_isentrope(
                &z,
                s,
                seed,
                IsentropeOptions {
                    minimum_pressure: -1e7,
                    points: count,
                    stop_at_spinodal: true,
                },
            )
            .unwrap();
        assert_eq!(curve.end, IsentropeEnd::Spinodal);
        for (i, p) in curve.points.iter().enumerate() {
            close(
                p.pressure_residual,
                0.,
                1.1e-10 * p.pressure.abs().max(1e5),
                0.,
                "physical P",
            );
            close(
                p.entropy_residual,
                0.,
                1.1e-10 * s.abs().max(1.),
                0.,
                "physical S",
            );
            let lambda = eos
                .prepare(p.temperature, &z)
                .unwrap()
                .mixture_stability(p.volume)
                .unwrap()
                .eigenvalue;
            if i + 1 == curve.points.len() {
                assert!(lambda.abs() < 1e-9);
            } else {
                assert!(lambda > 0.);
            }
            if i > 0 {
                assert!(p.pressure < curve.points[i - 1].pressure);
            }
        }
    }
    let curve = eos
        .homogeneous_isentrope(
            &z,
            s,
            seed,
            IsentropeOptions {
                minimum_pressure: 1e5,
                points: 2,
                stop_at_spinodal: true,
            },
        )
        .unwrap();
    assert_eq!(curve.end, IsentropeEnd::MinimumPressure);
    close(
        curve.points.last().unwrap().pressure,
        1e5,
        0.,
        0.,
        "minimum pressure",
    );
}

#[test]
fn local_unstable_and_zero_feed_roots_do_not_trigger_equilibrium_flash() {
    let eos = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
    let z = [1., 0.];
    let state = eos.prepare(150., &z).unwrap();
    let mut negative_pressure = false;
    let mut unstable = false;
    for v in [5e-5, 1.5e-4, 2e-3] {
        let c = state.caloric_tv(v).unwrap();
        negative_pressure |= c.pressure < 0.;
        unstable |= state.mixture_stability(v).unwrap().eigenvalue < 0.;
        let p = eos
            .homogeneous_ps(c.pressure, c.entropy, &z, 150. * 1.002, v * 0.998)
            .unwrap();
        close(p.temperature, 150., 1e-7, 1e-9, "local root T");
        close(p.volume, v, 1e-13, 1e-9, "local root V");
    }
    assert!(negative_pressure && unstable);
}

#[test]
fn mixture_isentrope_stops_at_local_composition_instability() {
    let eos = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
    let z = [0.5, 0.5];
    let state = eos.prepare(250., &z).unwrap();
    let v = state.volume(1e7, Phase::Liquid).unwrap();
    let s = state.caloric_tv(v).unwrap().entropy;
    let seed = eos.homogeneous_ps(1e7, s, &z, 250., v).unwrap();
    let curve = eos
        .homogeneous_isentrope(
            &z,
            s,
            seed,
            IsentropeOptions {
                minimum_pressure: -2e7,
                points: 8,
                stop_at_spinodal: true,
            },
        )
        .unwrap();
    assert_eq!(curve.end, IsentropeEnd::Spinodal);
    let end = curve.points.last().unwrap();
    let end_state = eos.prepare(end.temperature, &z).unwrap();
    assert!(
        end_state
            .mixture_stability(end.volume)
            .unwrap()
            .eigenvalue
            .abs()
            < 1e-9
    );
    // A mixture loses composition stability before its mechanical pressure
    // derivative vanishes; testing only dP/dV would miss this boundary.
    assert!(end_state.pressure(end.volume).unwrap().dv < 0.);
    close(
        end.pressure_residual,
        0.,
        1.1e-10 * end.pressure.abs().max(1e5),
        0.,
        "mixture P",
    );
    close(
        end.entropy_residual,
        0.,
        1.1e-10 * s.abs().max(1.),
        0.,
        "mixture S",
    );
}

#[test]
fn invalid_specs_and_point_counts_are_rejected() {
    let eos = Cubic::new(CubicKind::PengRobinson, &["C1"]).unwrap();
    for (p, s, t, v, z) in [
        (f64::NAN, 100., 200., 1e-3, vec![1.]),
        (1e5, f64::INFINITY, 200., 1e-3, vec![1.]),
        (1e5, 100., 0., 1e-3, vec![1.]),
        (1e5, 100., 200., -1., vec![1.]),
        (1e5, 100., 200., 1e-3, vec![0.]),
        (1e5, 100., 200., 1e-3, vec![1., 1.]),
    ] {
        assert!(eos.homogeneous_ps(p, s, &z, t, v).is_err());
    }
    assert!(
        eos.baseline_metastable_isentrope(&[1.], 6e6, 108.47763364585803, 1e5, 0)
            .is_err()
    );
    assert!(
        eos.baseline_metastable_isentrope(&[1.], 6e6, 108.47763364585803, 7e6, 12)
            .is_err()
    );
    let one = eos
        .baseline_metastable_isentrope(&[1.], 6e6, 108.47763364585803, 1e5, 1)
        .unwrap();
    assert_eq!((one.points.len(), one.end), (1, IsentropeEnd::PointLimit));
    let seed = one.points[0];
    for options in [
        IsentropeOptions {
            points: 1,
            ..Default::default()
        },
        IsentropeOptions {
            minimum_pressure: 6e6,
            ..Default::default()
        },
        IsentropeOptions {
            minimum_pressure: f64::NAN,
            ..Default::default()
        },
    ] {
        assert!(
            eos.homogeneous_isentrope(&[1.], seed.entropy, seed, options)
                .is_err()
        );
    }
}
