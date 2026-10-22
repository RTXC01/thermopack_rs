use thermopack_rs::spinodal_properties::{
    BaselineSpinodalExtrema, SpinodalProperty, baseline_pure_spinodal_extrema,
};
use thermopack_rs::{Cubic, CubicKind, CubicOptions};

fn property(code: &str) -> SpinodalProperty {
    match code {
        "1" => SpinodalProperty::Entropy,
        "2" => SpinodalProperty::Volume,
        "3" => SpinodalProperty::Enthalpy,
        "4" => SpinodalProperty::Temperature,
        "5" => SpinodalProperty::Pressure,
        _ => SpinodalProperty::InternalEnergy,
    }
}
fn model(c: &[&str]) -> Cubic {
    Cubic::with_options(
        match c[1] {
            "SRK" => CubicKind::SoaveRedlichKwong,
            "RK" => CubicKind::RedlichKwong,
            "VDW" => CubicKind::VanDerWaals,
            _ => CubicKind::PengRobinson,
        },
        &[c[2]],
        CubicOptions {
            volume_translation: c[3] == "translated",
            ..Default::default()
        },
    )
    .unwrap()
}
fn values(e: BaselineSpinodalExtrema) -> [f64; 6] {
    [
        e.minimum.temperature,
        e.minimum.volume,
        e.minimum.pressure,
        e.maximum.temperature,
        e.maximum.volume,
        e.maximum.pressure,
    ]
}
fn matches(a: &[f64], b: &[f64]) -> bool {
    a.iter()
        .zip(b)
        .enumerate()
        .all(|(j, (a, b))| (a - b).abs() <= [2e-5, 1e-11, 2.][j % 3] + 2e-6 * a.abs().max(b.abs()))
}

#[test]
fn all_finite_source_candidates_match_fresh_and_prepared_queries() {
    let mut group = String::new();
    let mut eos = None;
    let mut search = None;
    let mut counts = [0; 2];
    // Keep owned models outside the loop so the prepared workspace can borrow
    // one immutably across every property query for that grid.
    let mut models = Vec::new();
    let mut keys = Vec::new();
    for row in include_str!("fixtures/spinodal_extrema.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = row.split(';').collect();
        let c: Vec<_> = f[0].split('|').collect();
        let key = format!("{}|{}|{}|{}", c[1], c[2], c[3], c[5]);
        if keys.last() != Some(&key) {
            models.push(model(&c));
            keys.push(key);
        }
    }
    let mut index = 0;
    for row in include_str!("fixtures/spinodal_extrema.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = row.split(';').collect();
        let c: Vec<_> = f[0].split('|').collect();
        let key = format!("{}|{}|{}|{}", c[1], c[2], c[3], c[5]);
        if key != group {
            group = key;
            eos = Some(&models[index]);
            index += 1;
            search = Some(
                eos.unwrap()
                    .prepare_pure_spinodal_search(&[1.], c[5].parse().unwrap(), None)
                    .unwrap(),
            );
        }
        let model = eos.unwrap();
        let property = property(c[4].split(':').nth(1).unwrap());
        let expected: Vec<f64> = f[1..].iter().map(|s| s.parse().unwrap()).collect();
        let ambiguous = expected[12] != 0.;
        let alternative = [
            expected[3],
            expected[4],
            expected[13],
            expected[0],
            expected[1],
            expected[14],
        ];
        for actual in [
            baseline_pure_spinodal_extrema(model, &[1.], property, c[5].parse().unwrap()).unwrap(),
            search.as_ref().unwrap().extrema(property).unwrap(),
        ] {
            assert_eq!(actual.ordering_ambiguous, ambiguous, "{}", f[0]);
            let native = values(actual);
            assert!(
                matches(&native, &expected[..6]) || ambiguous && matches(&native, &alternative),
                "{} {:?} {:?}",
                f[0],
                native,
                expected
            );
            // Preserve and independently verify the source's mixed-coordinate
            // minimum pressure rather than treating it as physical closure.
            let mixed = model
                .prepare(actual.maximum.temperature, &[1.])
                .unwrap()
                .pressure(actual.minimum.volume)
                .unwrap()
                .value;
            assert!((mixed - actual.minimum.pressure).abs() < 1e-7);
        }
        counts[usize::from(ambiguous)] += 1;
    }
    assert_eq!(counts, [379, 83]);
}

#[test]
fn monotone_branches_return_first_points_and_report_pressure_coordinate_mismatch() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1"]).unwrap();
    let e =
        baseline_pure_spinodal_extrema(&model, &[1.], SpinodalProperty::Temperature, 1e6).unwrap();
    assert!(e.minimum.temperature < e.maximum.temperature);
    assert!(e.maximum.temperature < 0.99 * model.components()[0].critical_temperature);
    let actual = model
        .prepare(e.minimum.temperature, &[1.])
        .unwrap()
        .pressure(e.minimum.volume)
        .unwrap()
        .value;
    assert!((actual - e.minimum.pressure).abs() > 1e5);
    assert!(!e.ordering_ambiguous);
    let pressure =
        baseline_pure_spinodal_extrema(&model, &[1.], SpinodalProperty::Pressure, 1e6).unwrap();
    assert!(pressure.ordering_ambiguous);
}

#[test]
fn invalid_grids_and_nonpure_feeds_are_rejected() {
    let pure = Cubic::new(CubicKind::PengRobinson, &["C1"]).unwrap();
    for p in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(
            baseline_pure_spinodal_extrema(&pure, &[1.], SpinodalProperty::Entropy, p).is_err()
        );
    }
    for z in [&[0.][..], &[-1.], &[1., 1.]] {
        assert!(baseline_pure_spinodal_extrema(&pure, z, SpinodalProperty::Entropy, 1e6).is_err());
    }
    let mixed = Cubic::new(CubicKind::PengRobinson, &["C1", "C2"]).unwrap();
    assert!(
        baseline_pure_spinodal_extrema(&mixed, &[0.5, 0.5], SpinodalProperty::Entropy, 1e6)
            .is_err()
    );
}

#[test]
fn refined_caloric_candidates_are_stationary_at_the_corrected_spinodal_volume() {
    use thermopack_rs::mixture_spinodal::spinodal_volume_at_temperature;
    let model = Cubic::new(CubicKind::PengRobinson, &["C1"]).unwrap();
    for property in [
        SpinodalProperty::Entropy,
        SpinodalProperty::Enthalpy,
        SpinodalProperty::InternalEnergy,
    ] {
        let point = baseline_pure_spinodal_extrema(&model, &[1.], property, 1e6)
            .unwrap()
            .maximum;
        let dt = point.temperature * 1e-5;
        let low =
            spinodal_volume_at_temperature(&model, point.temperature - dt, &[1.], point.volume)
                .unwrap();
        let high =
            spinodal_volume_at_temperature(&model, point.temperature + dt, &[1.], point.volume)
                .unwrap();
        let a = property.evaluate(&model, &[1.], low).unwrap();
        let b = property.evaluate(&model, &[1.], high).unwrap();
        let gradient = (b - a) / (2. * dt);
        assert!(
            gradient.abs() < 1e-5 * a.abs().max(b.abs()).max(1.) / point.temperature,
            "{property:?}: {gradient}"
        );
        // Reported source volumes need not themselves lie on the spinodal.
        assert!(
            model
                .prepare(point.temperature, &[1.])
                .unwrap()
                .mixture_stability(point.volume)
                .unwrap()
                .eigenvalue
                .abs()
                > 1e-6
        );
    }
}
