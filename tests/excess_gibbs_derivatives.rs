use thermopack_rs::{Cubic, CubicKind, CubicOptions, MixingRule, Phase};

#[test]
fn every_qualified_activity_model_tp_derivative_matches_the_source() {
    let mut count = 0;
    for row in include_str!("fixtures/ge_fugacity_derivatives.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = row.split(';').collect();
        let c: Vec<_> = f[0].split('|').collect();
        let s: Vec<_> = c[5].split(':').collect();
        let names: Vec<_> = c[2].split(',').collect();
        let z: Vec<f64> = c[4].split(',').map(|s| s.parse().unwrap()).collect();
        let n = z.len();
        let mixing_rule = match s[3] {
            "HV1" => MixingRule::HuronVidal1,
            "HV2" => MixingRule::HuronVidal2,
            "NRTL" => MixingRule::Nrtl,
            _ => MixingRule::Unifac,
        };
        let model = Cubic::with_options(
            if c[1] == "SRK" {
                CubicKind::SoaveRedlichKwong
            } else {
                CubicKind::PengRobinson
            },
            &names,
            CubicOptions {
                mixing_rule,
                volume_translation: c[3] == "translated",
                parameter_reference: s[4],
                ..Default::default()
            },
        )
        .unwrap();
        let prepared = model.prepare(s[0].parse().unwrap(), &z).unwrap();
        let v = prepared
            .volume(
                s[1].parse().unwrap(),
                if s[2] == "L" {
                    Phase::Liquid
                } else {
                    Phase::Vapor
                },
            )
            .unwrap();
        let mut dt = vec![0.; n];
        let mut dp = vec![0.; n];
        let mut dn = vec![0.; n * n];
        prepared
            .fugacity_derivatives_into(v, &mut dt, &mut dp, &mut dn)
            .unwrap();
        for (i, (a, b)) in dt
            .into_iter()
            .chain(dp)
            .chain(dn)
            .zip(f[1..].iter().map(|s| s.parse::<f64>().unwrap()))
            .enumerate()
        {
            let floor = if i < n {
                1e-9
            } else if i < 2 * n {
                1e-11
            } else {
                1e-8
            };
            assert!(
                (a - b).abs() <= floor + 2e-7 * a.abs().max(b.abs()),
                "{} {i}: {a} vs {b}",
                f[0]
            );
        }
        count += 1;
    }
    assert_eq!(count, 252);
}

#[test]
fn activity_mixing_fugacity_mole_derivatives_match_independent_differences() {
    for (kind, names, mixing_rule, reference) in [
        (
            CubicKind::PengRobinson,
            vec!["C1", "C3"],
            MixingRule::Unifac,
            "DEFAULT",
        ),
        (
            CubicKind::SoaveRedlichKwong,
            vec!["C1", "C3"],
            MixingRule::Unifac,
            "DEFAULT",
        ),
        (
            CubicKind::PengRobinson,
            vec!["C1", "C2", "C3"],
            MixingRule::Unifac,
            "DEFAULT",
        ),
        (
            CubicKind::SoaveRedlichKwong,
            vec!["C1", "H2O"],
            MixingRule::HuronVidal1,
            "DEFAULT",
        ),
        (
            CubicKind::SoaveRedlichKwong,
            vec!["C1", "H2O"],
            MixingRule::HuronVidal2,
            "DEFAULT",
        ),
        (
            CubicKind::PengRobinson,
            vec!["H2O", "MEG"],
            MixingRule::HuronVidal1,
            "Hemmingsen2011",
        ),
        (
            CubicKind::PengRobinson,
            vec!["H2S", "C3"],
            MixingRule::Nrtl,
            "Dicko2012",
        ),
    ] {
        let model = Cubic::with_options(
            kind,
            &names,
            CubicOptions {
                mixing_rule,
                parameter_reference: reference,
                ..Default::default()
            },
        )
        .unwrap();
        let n = names.len();
        let z = vec![1. / n as f64; n];
        let (t, p) = (300., 1e5);
        for phase in [Phase::Liquid, Phase::Vapor] {
            let state = model.prepare(t, &z).unwrap();
            let v = state.volume(p, phase).unwrap();
            let mut dt = vec![0.; n];
            let mut dp = vec![0.; n];
            let mut dn = vec![0.; n * n];
            state
                .fugacity_derivatives_into(v, &mut dt, &mut dp, &mut dn)
                .unwrap();
            for j in 0..n {
                let h = 1e-5;
                let mut plus = z.clone();
                let mut minus = z.clone();
                plus[j] += h;
                minus[j] -= h;
                let plus_sum = plus.iter().sum::<f64>();
                let minus_sum = minus.iter().sum::<f64>();
                for v in &mut plus {
                    *v /= plus_sum;
                }
                for v in &mut minus {
                    *v /= minus_sum;
                }
                let mut phi_plus = vec![0.; n];
                let mut phi_minus = vec![0.; n];
                let a = model.prepare(t, &plus).unwrap();
                let b = model.prepare(t, &minus).unwrap();
                a.ln_fugacity_tp_into(p, a.volume(p, phase).unwrap(), &mut phi_plus)
                    .unwrap();
                b.ln_fugacity_tp_into(p, b.volume(p, phase).unwrap(), &mut phi_minus)
                    .unwrap();
                for i in 0..n {
                    let fd = (phi_plus[i] - phi_minus[i]) / (2. * h);
                    assert!(
                        (dn[i * n + j] - fd).abs() < 2e-6 + 2e-6 * fd.abs(),
                        "{kind:?} {names:?} {mixing_rule:?} {phase:?} {i},{j}: {} vs {fd}",
                        dn[i * n + j]
                    );
                    assert!((dn[i * n + j] - dn[j * n + i]).abs() < 1e-9);
                }
            }
            for i in 0..n {
                assert!(
                    z.iter()
                        .enumerate()
                        .map(|(j, z)| z * dn[i * n + j])
                        .sum::<f64>()
                        .abs()
                        < 1e-9
                );
            }
        }
    }
}
