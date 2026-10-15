use thermopack_rs::{Cubic, CubicKind, FlashPhase, Phase};

fn kind(name: &str) -> CubicKind {
    match name {
        "PR" => CubicKind::PengRobinson,
        "SRK" => CubicKind::SoaveRedlichKwong,
        "RK" => CubicKind::RedlichKwong,
        "VDW" => CubicKind::VanDerWaals,
        _ => panic!("unknown fixture model"),
    }
}

#[test]
fn fugacity_derivatives_match_reference() {
    let mut count = 0;
    for line in include_str!("fixtures/fugacity_derivatives.csv")
        .lines()
        .skip(1)
    {
        let columns: Vec<_> = line.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let names: Vec<_> = fields[1].split(',').collect();
        let model = Cubic::new(kind(fields[0]), &names).unwrap();
        let x: Vec<f64> = fields[4].split(',').map(|s| s.parse().unwrap()).collect();
        let n = x.len();
        let state = model.prepare(fields[2].parse().unwrap(), &x).unwrap();
        let v = state
            .volume(
                fields[3].parse().unwrap(),
                if fields[5].ends_with('L') {
                    Phase::Liquid
                } else {
                    Phase::Vapor
                },
            )
            .unwrap();
        let (mut dt, mut dp, mut dn) = (vec![0.0; n], vec![0.0; n], vec![0.0; n * n]);
        state
            .fugacity_derivatives_into(v, &mut dt, &mut dp, &mut dn)
            .unwrap();
        let values: Vec<_> = dt.into_iter().chain(dp).chain(dn).collect();
        assert_eq!(values.len(), columns.len() - 1);
        for (i, (a, b)) in values.iter().zip(&columns[1..]).enumerate() {
            let b: f64 = b.parse().unwrap();
            let abs = if i < n {
                1e-9
            } else if i < 2 * n {
                1e-11
            } else {
                1e-8
            };
            assert!(
                (a - b).abs() <= abs + 2e-7 * a.abs().max(b.abs()),
                "{} quantity {i}: {a} {b}",
                columns[0]
            );
        }
        count += 1;
    }
    assert_eq!(count, 1320);
}

#[test]
fn tp_flash_matches_reference_physical_splits() {
    let mut count = 0;
    let mut permutations = 0;
    for line in include_str!("fixtures/tp_flash.csv").lines().skip(1) {
        let columns: Vec<_> = line.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let names: Vec<_> = fields[1].split(',').collect();
        let model = Cubic::new(kind(fields[0]), &names).unwrap();
        let t: f64 = fields[2].parse().unwrap();
        let p: f64 = fields[3].parse().unwrap();
        let z: Vec<f64> = fields[4].split(',').map(|s| s.parse().unwrap()).collect();
        let result = model
            .tp_flash(t, p, &z)
            .unwrap_or_else(|e| panic!("{}: {e}", columns[0]));
        let reference: Vec<f64> = columns[2..].iter().map(|s| s.parse().unwrap()).collect();
        let mut actual = vec![
            result.phase as i32 as f64,
            result.vapor_fraction,
            result.liquid_fraction,
        ];
        actual.extend(result.liquid.iter().chain(&result.vapor));
        let matches = |values: &[f64]| {
            values
                .iter()
                .zip(&reference)
                .enumerate()
                .all(|(i, (a, b))| {
                    (a - b).abs()
                        <= if i == 0 {
                            0.0
                        } else {
                            2e-7 + 2e-7 * a.abs().max(b.abs())
                        }
                })
        };
        if columns[1] == "1" && result.phase == FlashPhase::TwoPhase && !matches(&actual) {
            // Permit exchange only with independent evidence that both phase
            // compositions have a unique mechanically stable volume root.
            for composition in [&result.liquid, &result.vapor] {
                assert_eq!(
                    model.prepare(t, composition).unwrap().volumes(p).unwrap().1,
                    1
                );
            }
            actual = vec![
                result.phase as i32 as f64,
                result.liquid_fraction,
                result.vapor_fraction,
            ];
            actual.extend(result.vapor.iter().chain(&result.liquid));
            permutations += 1;
        }
        assert_eq!(actual.len(), reference.len());
        assert!(
            matches(&actual),
            "{}\nRust {actual:?}\nreference {reference:?}",
            columns[0]
        );
        if result.phase == FlashPhase::TwoPhase {
            for (i, zi) in z.iter().enumerate() {
                assert!(
                    (zi - result.liquid_fraction * result.liquid[i]
                        - result.vapor_fraction * result.vapor[i])
                        .abs()
                        < 1e-10
                );
            }
        }
        count += 1;
    }
    assert_eq!(count, 1043);
    assert_eq!(permutations, 1);
}
