use thermopack_rs::SolidWater;

fn close(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a} != {b}"
    );
}

#[test]
fn all_ice_source_functions_and_reference_sequences_match() {
    let mut count = 0;
    for row in include_str!("fixtures/solid_water.csv").lines().skip(1) {
        let c: Vec<_> = row.split(';').collect();
        let f: Vec<_> = c[0].split('|').collect();
        let mut model = SolidWater::new();
        if f[1] != "raw" {
            for pair in f[1].split(':') {
                let v: Vec<f64> = pair.split(',').map(|s| s.parse().unwrap()).collect();
                model.initialize(v[0], v[1]).unwrap();
            }
        }
        let t = f[2].parse().unwrap();
        let p = f[3].parse().unwrap();
        let s = model.state(t, p).unwrap();
        let g = s.gibbs;
        let actual = [
            g.value,
            g.dp,
            g.dpp,
            g.dt,
            g.dtt,
            g.dtp,
            s.volume,
            s.entropy,
            s.enthalpy,
            s.internal_energy,
            s.helmholtz_energy,
            s.baseline_heat_capacity_p,
        ];
        let scalar = [
            model.gibbs(t, p).unwrap(),
            model.dg_dp(t, p).unwrap(),
            model.d2g_dp2(t, p).unwrap(),
            model.dg_dt(t, p).unwrap(),
            model.d2g_dt2(t, p).unwrap(),
            model.d2g_dtdp(t, p).unwrap(),
            model.specific_volume(t, p).unwrap(),
            model.entropy(t, p).unwrap(),
            model.enthalpy(t, p).unwrap(),
            model.internal_energy(t, p).unwrap(),
            model.helmholtz(t, p).unwrap(),
            model.baseline_heat_capacity(t, p).unwrap(),
        ];
        let atol = [
            1e-6, 1e-16, 1e-25, 1e-9, 1e-12, 1e-17, 1e-16, 1e-9, 1e-6, 1e-6, 1e-6, 1e-10,
        ];
        for (i, a) in actual.into_iter().enumerate() {
            close(a, c[i + 1].parse().unwrap(), atol[i], 1e-11);
            close(scalar[i], c[i + 1].parse().unwrap(), atol[i], 1e-11);
        }
        count += 1;
    }
    assert_eq!(count, 780);
}

#[test]
fn ice_gibbs_derivatives_close_independent_differences_and_caloric_identities() {
    let mut model = SolidWater::new();
    model.initialize(100., -10000.).unwrap();
    for t in [100., 200., 250., 273.15] {
        for p in [1e5, 1e8] {
            let q = model.state(t, p).unwrap();
            let dt = 0.001;
            let dp = 1000.;
            let hi = model.state(t + dt, p).unwrap();
            let lo = model.state(t - dt, p).unwrap();
            close(
                q.gibbs.dt,
                (hi.gibbs.value - lo.gibbs.value) / (2. * dt),
                1e-7,
                1e-7,
            );
            close(
                q.gibbs.dtt,
                (hi.gibbs.dt - lo.gibbs.dt) / (2. * dt),
                1e-10,
                1e-7,
            );
            close(
                q.gibbs.dtp,
                (hi.volume - lo.volume) / (2. * dt),
                1e-17,
                1e-6,
            );
            close(
                q.heat_capacity_p,
                (hi.enthalpy - lo.enthalpy) / (2. * dt),
                1e-6,
                1e-7,
            );
            let hp = model.state(t, p + dp).unwrap();
            let lp = model.state(t, p - dp).unwrap();
            close(
                q.gibbs.dp,
                (hp.gibbs.value - lp.gibbs.value) / (2. * dp),
                1e-13,
                1e-7,
            );
            close(
                q.gibbs.dpp,
                (hp.volume - lp.volume) / (2. * dp),
                1e-24,
                1e-6,
            );
            close(
                q.gibbs.dtp,
                (-hp.entropy + lp.entropy) / (2. * dp),
                1e-17,
                1e-6,
            );
            close(q.gibbs.value, q.enthalpy - t * q.entropy, 1e-10, 1e-13);
            close(q.internal_energy, q.enthalpy - p * q.volume, 1e-10, 1e-13);
            close(
                q.helmholtz_energy,
                q.internal_energy - t * q.entropy,
                1e-10,
                1e-13,
            );
            close(
                q.baseline_heat_capacity_p,
                q.heat_capacity_p * 18.01528e-3,
                1e-15,
                1e-14,
            );
            assert!(q.heat_capacity_p > q.baseline_heat_capacity_p * 50.);
        }
    }
}

#[test]
fn ice_reference_alignment_is_owned_and_invalid_shifts_are_transactional() {
    let original = SolidWater::new();
    let mut shifted = original;
    shifted.initialize(185., -16000.).unwrap();
    let (t, p) = (
        SolidWater::REFERENCE_TEMPERATURE,
        SolidWater::REFERENCE_PRESSURE,
    );
    close(
        shifted.entropy(t, p).unwrap(),
        185. - SolidWater::MELTING_ENTROPY,
        1e-10,
        1e-12,
    );
    close(shifted.gibbs(t, p).unwrap(), -16000., 1e-9, 1e-12);
    let first = shifted.state(250., 1e5).unwrap();
    shifted.initialize(185., -16000.).unwrap();
    close(
        shifted.gibbs(250., 1e5).unwrap(),
        first.gibbs.value,
        1e-9,
        1e-12,
    );
    close(
        shifted.heat_capacity(250., 1e5).unwrap(),
        original.heat_capacity(250., 1e5).unwrap(),
        1e-12,
        1e-12,
    );
    let before = shifted;
    assert!(shifted.initialize(f64::NAN, 0.).is_err());
    assert_eq!(shifted, before);
    assert!(shifted.initialize(0., f64::INFINITY).is_err());
    assert_eq!(shifted, before);
    assert_ne!(original, shifted);
    for t in [-1., f64::NAN, f64::INFINITY] {
        assert!(shifted.prepare(t).is_err());
    }
    assert!(shifted.state(250., f64::NAN).is_err());
    assert!(shifted.state(250., f64::MAX).is_err());
    let zero = shifted.state(0., p).unwrap();
    assert_eq!(zero.heat_capacity_p, 0.);
    let prepared = shifted.prepare(250.).unwrap();
    for p in [0., 1e5, 1e8] {
        close(
            prepared.state(p).unwrap().gibbs.value,
            shifted.gibbs(250., p).unwrap(),
            0.,
            0.,
        );
    }
}
