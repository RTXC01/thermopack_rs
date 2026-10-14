use thermopack_rs::{Cubic, CubicKind, Phase, R, component, components};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}

#[test]
fn database_and_input_validation() {
    assert_eq!(components().len(), 109);
    assert_eq!(component("Methane").unwrap().ident, "C1");
    assert_eq!(component(" ch4 ").unwrap().ident, "C1");
    assert!(component("unobtainium").is_err());
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "CO2"]).unwrap();
    assert!(model.prepare(0.0, &[0.5, 0.5]).is_err());
    assert!(model.prepare(f64::NAN, &[0.5, 0.5]).is_err());
    assert!(model.prepare(300.0, &[0.1, 0.1]).is_err());
    assert!(model.prepare(300.0, &[1.0, -0.0, 0.0]).is_err());
    close(model.kij(0, 1).unwrap(), 0.092, 1e-14);
    let state = model.prepare(300.0, &[0.5, 0.5]).unwrap();
    assert!(state.pressure(state.covolume()).is_err());
    assert!(state.volume(-1.0, Phase::Vapor).is_err());
}

#[test]
fn pressure_derivatives_and_root_closure() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        let model = Cubic::new(kind, &["C1", "CO2", "N2"]).unwrap();
        for t in [180.0, 300.0, 600.0] {
            let state = model.prepare(t, &[0.6, 0.3, 0.1]).unwrap();
            for v in [state.covolume() * 1.5, 0.001, 0.1] {
                let p = state.pressure(v).unwrap();
                let dt = t * 1e-5;
                let dv = v * 1e-5;
                let pt = (model
                    .prepare(t + dt, &[0.6, 0.3, 0.1])
                    .unwrap()
                    .pressure(v)
                    .unwrap()
                    .value
                    - model
                        .prepare(t - dt, &[0.6, 0.3, 0.1])
                        .unwrap()
                        .pressure(v)
                        .unwrap()
                        .value)
                    / (2.0 * dt);
                let pv = (state.pressure(v + dv).unwrap().value
                    - state.pressure(v - dv).unwrap().value)
                    / (2.0 * dv);
                close(p.dt, pt, 1e-7);
                close(p.dv, pv, 1e-7);
            }
            for p in [1.0, 1e3, 1e5, 1e6, 1e7, 1e8] {
                for phase in [Phase::Liquid, Phase::Vapor, Phase::MinimumGibbs] {
                    let v = state.volume(p, phase).unwrap();
                    close(state.pressure(v).unwrap().value, p, 2e-8);
                }
            }
        }
    }
}

#[test]
fn caloric_derivatives_and_gibbs_identity() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
    for t in [200.0, 300.0, 500.0] {
        let x = [0.7, 0.3];
        let state = model.prepare(t, &x).unwrap();
        let p = 1e6;
        let v = state.volume(p, Phase::Vapor).unwrap();
        let (cp, cv) = state.heat_capacities(v).unwrap();
        let dt = 0.001;
        let plus = model.prepare(t + dt, &x).unwrap();
        let minus = model.prepare(t - dt, &x).unwrap();
        let vp = plus.volume(p, Phase::Vapor).unwrap();
        let vm = minus.volume(p, Phase::Vapor).unwrap();
        close(
            cp,
            (plus.enthalpy(vp).unwrap() - minus.enthalpy(vm).unwrap()) / (2.0 * dt),
            2e-7,
        );
        close(
            cp / t,
            (plus.entropy(vp).unwrap() - minus.entropy(vm).unwrap()) / (2.0 * dt),
            2e-7,
        );
        close(
            cv,
            ((plus.enthalpy(v).unwrap() - plus.pressure(v).unwrap().value * v)
                - (minus.enthalpy(v).unwrap() - minus.pressure(v).unwrap().value * v))
                / (2.0 * dt),
            2e-7,
        );
        let residual = state.residual(v).unwrap();
        let mut phi = [0.0; 2];
        state.ln_fugacity_into(v, &mut phi).unwrap();
        close(
            residual.enthalpy - t * residual.entropy,
            R * t * (x[0] * phi[0] + x[1] * phi[1]),
            1e-10,
        );
        assert!(state.speed_of_sound(v).unwrap() > 0.0);
    }
}

#[test]
fn independent_models_can_run_on_threads() {
    let model = Cubic::new(CubicKind::PengRobinson, &["CO2", "C1"]).unwrap();
    let expected = model
        .prepare(300.0, &[0.8, 0.2])
        .unwrap()
        .volume(5e6, Phase::Vapor)
        .unwrap();
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..8)
            .map(|_| {
                scope.spawn(|| {
                    model
                        .prepare(300.0, &[0.8, 0.2])
                        .unwrap()
                        .volume(5e6, Phase::Vapor)
                        .unwrap()
                })
            })
            .collect();
        for handle in handles {
            assert_eq!(handle.join().unwrap(), expected);
        }
    });
}

#[test]
fn rachford_rice_closure_and_polynomial_roots() {
    use thermopack_rs::numerics::{cubic_roots, rachford_rice};
    let (roots, n) = cubic_roots(-6.0, 11.0, -6.0);
    assert_eq!(n, 3);
    for (a, b) in roots.into_iter().zip([1.0, 2.0, 3.0]) {
        close(a, b, 1e-13);
    }
    let small = 1e-10;
    let (roots, n) = cubic_roots(
        -(1.0 + 3.0 * small),
        3.0 * small + 2.0 * small * small,
        -2.0 * small * small,
    );
    assert_eq!(n, 3);
    for (a, b) in roots.into_iter().zip([small, 2.0 * small, 1.0]) {
        assert!((a - b).abs() < 1e-8 * b.abs());
    }
    let z = [0.4, 0.6];
    let k = [3.0, 0.2];
    let beta = rachford_rice(&z, &k).unwrap();
    close(beta, 0.2, 1e-13);
    let liquid: f64 = z
        .iter()
        .zip(k)
        .map(|(zi, ki)| zi / (1.0 + beta * (ki - 1.0)))
        .sum();
    let vapor: f64 = z
        .iter()
        .zip(k)
        .map(|(zi, ki)| zi * ki / (1.0 + beta * (ki - 1.0)))
        .sum();
    close(liquid, 1.0, 1e-13);
    close(vapor, 1.0, 1e-13);
}

#[test]
fn pure_saturation_fugacity_and_clapeyron() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        for name in ["C1", "CO2", "C3"] {
            let model = Cubic::new(kind, &[name]).unwrap();
            let tc = model.components()[0].critical_temperature;
            for tr in [0.65, 0.8, 0.95, 0.999] {
                let t = tc * tr;
                let sat = model
                    .saturation_pressure(t)
                    .unwrap_or_else(|e| panic!("{kind:?} {name} Tr={tr}: {e}"));
                assert!(sat.liquid_volume < sat.vapor_volume);
                let state = model.prepare(t, &[1.0]).unwrap();
                let (mut l, mut v) = ([0.0], [0.0]);
                state.ln_fugacity_into(sat.liquid_volume, &mut l).unwrap();
                state.ln_fugacity_into(sat.vapor_volume, &mut v).unwrap();
                close(l[0], v[0], 1e-9);
                let dt = tc * 1e-6;
                let slope = (model.saturation_pressure(t + dt).unwrap().pressure
                    - model.saturation_pressure(t - dt).unwrap().pressure)
                    / (2.0 * dt);
                close(sat.dpdt, slope, 1e-5);
            }
            assert!(model.saturation_pressure(tc * 1.01).is_err());
        }
    }
}
