use thermopack_rs::{EntropyReference, IdealGas, Mbwr, MbwrFluid, R, component};
#[path = "../rust/support/mbwr_caloric_coordinates.rs"]
mod protocol;

fn floors(t: f64, v: f64, n: f64) -> Vec<f64> {
    let rt = R * t;
    let e = 2e-10;
    let mut s = vec![
        n,
        n / t,
        n / v,
        1.,
        n / t.powi(2),
        n / t / v,
        1. / t,
        n / v.powi(2),
        1. / v,
        1. / n,
    ];
    s.extend([n * rt / v, n * R / v, n * rt / v.powi(2), rt / v]);
    let energy = [n * rt, n * R, n * rt / v, rt];
    s.extend(energy);
    s.extend(energy);
    s.extend([n * R, n * R / t, n * R / v, R]);
    s.extend(energy);
    s.extend(energy);
    s.extend([rt, R, rt / v, rt / n]);
    s.push(n * R);
    let mut result: Vec<_> = s.into_iter().map(|s| s * e).collect();
    result.extend([0., e * n * R, 0., e, 0., 2e-8]);
    result
}

#[test]
fn total_tv_source_and_physical_views_replay_for_every_fit() {
    let mut count = 0;
    for line in include_str!("fixtures/mbwr_caloric.csv").lines().skip(1) {
        let row: Vec<_> = line.split(';').collect();
        let case: Vec<_> = row[0].split('|').collect();
        let t = case[1].parse().unwrap();
        let v = case[2].parse().unwrap();
        let n = case[3].parse().unwrap();
        let mut model = if case[4] == "DEFAULT" {
            MbwrFluid::from_record(case[0]).unwrap()
        } else {
            let c = component(case[4]).unwrap();
            MbwrFluid::with_ideal_gas(
                Mbwr::from_record(case[0]).unwrap(),
                IdealGas::with_integration_constants(c, 0., 0.).unwrap(),
                c.molar_mass * 1e-3,
            )
            .unwrap()
        };
        let profile: usize = case[5].parse().unwrap();
        if profile == 1 || profile == 3 {
            model
                .ideal_gas_mut()
                .set_formation_enthalpy(-12345.678)
                .unwrap();
            model
                .ideal_gas_mut()
                .set_standard_entropy(210.25, EntropyReference::Atmosphere)
                .unwrap();
        }
        if profile == 2 || profile == 3 {
            model.ideal_gas_mut().set_heat_capacity(8, &[3.5]).unwrap();
        }
        let margin = floors(t, v, n);
        for (column, result) in [
            model.baseline_caloric_tv(t, v, n),
            model.caloric_tv(t, v, n),
        ]
        .into_iter()
        .enumerate()
        {
            let q = result.unwrap_or_else(|e| panic!("{}: {e}", row[0]));
            let actual = protocol::coordinates(&q);
            let expected: Vec<f64> = row[column + 1]
                .split(',')
                .map(|s| s.parse().unwrap())
                .collect();
            assert_eq!(actual.len(), 45);
            assert_eq!(expected.len(), 45);
            assert_eq!(margin.len(), 45);
            for (i, ((a, b), floor)) in actual.iter().zip(expected).zip(&margin).enumerate() {
                assert!(
                    a.is_finite() && (a - b).abs() <= floor + 2e-7 * a.abs().max(b.abs()),
                    "{} view {column} field {i}: {a} != {b}",
                    row[0]
                );
            }
        }
        count += 1;
    }
    assert_eq!(count, 5994);
}
