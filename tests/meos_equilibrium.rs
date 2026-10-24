use std::collections::BTreeMap;
use thermopack_rs::{LennardJonesKind, Multiparameter, Phase};

fn model(config: &str) -> Multiparameter {
    let (e, c) = config.split_once(':').unwrap();
    match e {
        "NIST_MEOS" => Multiparameter::new(c).unwrap(),
        "LJ_MEOS" => Multiparameter::lennard_jones(c, LennardJonesKind::Full).unwrap(),
        "LJTS_MEOS" => {
            Multiparameter::lennard_jones(c, LennardJonesKind::TruncatedShifted).unwrap()
        }
        "LJTS_REDUCED" => Multiparameter::reduced_lennard_jones_truncated_shifted(),
        _ => panic!("unknown model {config}"),
    }
}
fn close(actual: f64, expected: f64, tolerance: f64, context: &str) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{context}: {actual} != {expected}, tolerance {tolerance}"
    );
}
fn groups(fixture: &str) -> BTreeMap<String, Vec<Vec<&str>>> {
    let mut groups = BTreeMap::<String, Vec<Vec<&str>>>::new();
    for line in fixture.lines().skip(1) {
        let row: Vec<_> = line.split(';').collect();
        let case: Vec<_> = row[0].split('|').collect();
        groups.entry(case[1..4].join("|")).or_default().push(row);
    }
    groups
}
#[test]
fn analytic_pressure_density_jets_match_independent_source_differences() {
    let mut count = 0;
    for line in include_str!("fixtures/meos_pressure_jet.csv")
        .lines()
        .skip(1)
    {
        let row: Vec<_> = line.split(';').collect();
        let c: Vec<_> = row[0].split('|').collect();
        let m = model(c[1]);
        let q = m
            .prepare(c[2].parse().unwrap())
            .unwrap()
            .pressure_density_jet(c[3].parse().unwrap())
            .unwrap();
        let actual = [
            q.value,
            q.drho,
            q.drhorho,
            q.drhorhorho,
            q.dt,
            q.dtrho,
            q.dtrhorho,
        ];
        for (i, value) in actual.into_iter().enumerate() {
            close(
                value,
                row[1 + i].parse().unwrap(),
                row[8 + i].parse().unwrap(),
                row[0],
            );
        }
        count += 1;
    }
    assert_eq!(count, 672);
}
#[test]
fn local_and_default_critical_points_match_source_pressure_equations() {
    let mut count = 0;
    for line in include_str!("fixtures/meos_critical.csv").lines().skip(1) {
        let row: Vec<_> = line.split(';').collect();
        let c: Vec<_> = row[0].split('|').collect();
        let m = model(c[1]);
        let q = if c[5] == "CRITICAL_AUTO" {
            m.critical_point()
        } else {
            m.critical_point_from(c[2].parse().unwrap(), c[3].parse().unwrap())
        }
        .unwrap();
        for (i, value) in [q.temperature, q.pressure, 1. / q.volume, q.compressibility]
            .into_iter()
            .enumerate()
        {
            close(
                value,
                row[1 + i].parse().unwrap(),
                row[5 + i].parse().unwrap(),
                row[0],
            );
        }
        let j = m
            .prepare(q.temperature)
            .unwrap()
            .pressure_density_jet(1. / q.volume)
            .unwrap();
        assert!(j.drhorhorho > 0. && j.value > 0.);
        count += 1;
    }
    assert_eq!(count, 52);
}
#[test]
fn full_stable_root_sets_and_phase_choices_match_source_kernels() {
    let mut count = 0;
    for (group, rows) in groups(include_str!("fixtures/meos_pressure_branches.csv")) {
        let f: Vec<_> = group.split('|').collect();
        let m = model(f[0]);
        let t = f[1].parse().unwrap();
        let maximum = f[2].parse().unwrap();
        let prepared = m.prepare(t).unwrap();
        let eq = prepared.equilibrium(maximum).unwrap();
        let mut buffer = vec![0.; eq.maximum_root_count()];
        for row in rows {
            let p = row[0]
                .split('|')
                .nth(5)
                .unwrap()
                .strip_prefix("ROOTS:")
                .unwrap()
                .parse()
                .unwrap();
            let expected: Vec<f64> = row[1].split(',').map(|x| x.parse().unwrap()).collect();
            let n = eq
                .densities_into(p, &mut buffer)
                .unwrap_or_else(|e| panic!("{}: {e}", row[0]));
            assert_eq!(n, expected.len(), "{}", row[0]);
            for (&a, &b) in buffer[..n].iter().zip(&expected) {
                close(a, b, 2e-8 * m.rc() + 2e-7 * a.abs().max(b.abs()), row[0]);
                let (actual_p, stiffness) = prepared.pressure_density(a).unwrap();
                assert!(stiffness > 0. && a <= maximum);
                close(
                    actual_p,
                    p,
                    2e-8 * p + 256. * f64::EPSILON * m.gas_constant() * t * a,
                    row[0],
                );
            }
            for (i, phase) in [Phase::Vapor, Phase::Liquid, Phase::MinimumGibbs]
                .into_iter()
                .enumerate()
            {
                let v = eq.volume(p, phase).unwrap();
                let expected: f64 = row[i + 2].parse().unwrap();
                close(v, expected, 1e-13 + 3e-7 * expected.abs(), row[0]);
            }
            count += 1;
        }
    }
    assert_eq!(count, 804);
}
#[test]
fn all_independent_source_spinodals_have_native_enclosures() {
    let mut count = 0;
    for (group, rows) in groups(include_str!("fixtures/meos_spinodals.csv")) {
        let f: Vec<_> = group.split('|').collect();
        let m = model(f[0]);
        let eq = m
            .prepare_equilibrium(f[1].parse().unwrap(), f[2].parse().unwrap())
            .unwrap();
        for row in rows {
            let rho: f64 = row[1].parse().unwrap();
            let tolerance = 2e-7 * m.rc();
            assert!(
                eq.spinodal_brackets()
                    .iter()
                    .any(|b| b.minimum_density - tolerance <= rho
                        && rho <= b.maximum_density + tolerance),
                "{}: {rho}",
                row[0]
            );
            count += 1;
        }
    }
    assert_eq!(count, 212);
}
#[test]
fn unrepresentable_extrapolated_root_sets_remain_explicit_errors() {
    let mut count = 0;
    for line in include_str!("fixtures/meos_pressure_unresolved.csv")
        .lines()
        .skip(1)
    {
        let row: Vec<_> = line.split(';').collect();
        let c: Vec<_> = row[0].split('|').collect();
        let m = model(c[1]);
        let eq = m
            .prepare_equilibrium(c[2].parse().unwrap(), c[3].parse().unwrap())
            .unwrap();
        let p = c[5].strip_prefix("ROOTS:").unwrap().parse().unwrap();
        assert!(eq.densities(p).is_err(), "{}", row[0]);
        for (i, phase) in [Phase::Vapor, Phase::Liquid, Phase::MinimumGibbs]
            .into_iter()
            .enumerate()
        {
            let volume = eq.volume(p, phase);
            if row[4 + i].is_empty() {
                assert!(volume.is_err(), "{}", row[0]);
            } else {
                let expected: f64 = row[4 + i].parse().unwrap();
                close(
                    volume.unwrap(),
                    expected,
                    1e-13 + 3e-7 * expected.abs(),
                    row[0],
                );
            }
        }
        count += 1;
    }
    assert_eq!(count, 8);
}
#[test]
fn bounded_domains_buffers_spinodal_windows_and_invalid_inputs_are_explicit() {
    let m = Multiparameter::new("C3").unwrap();
    for value in [0., -1., f64::INFINITY, f64::NAN] {
        assert!(m.prepare_equilibrium(280., value).is_err());
        assert!(m.prepare_equilibrium(value, 8. * m.rc()).is_err());
        assert!(m.critical_point_from(value, m.rc()).is_err());
        assert!(m.critical_point_from(m.tc(), value).is_err());
    }
    let eq = m.prepare_equilibrium(280., 8. * m.rc()).unwrap();
    assert!(eq.densities_into(1e6, &mut []).is_err());
    let mut checked = 0;
    for b in eq.spinodal_brackets() {
        let p = (b.minimum_pressure + b.maximum_pressure) * 0.5;
        if p > 0. {
            assert!(eq.densities(p).is_err());
            checked += 1;
        }
    }
    assert!(checked > 0);
    let small = m.prepare_equilibrium(280., 0.01 * m.rc()).unwrap();
    assert!(small.densities(1e9).unwrap().is_empty());
    assert!(small.volume(1e9, Phase::MinimumGibbs).is_err());
    let t = 1.5 * m.tc();
    let maximum = 0.2 * m.rc();
    let p = m.prepare(t).unwrap().pressure_density(maximum).unwrap().0;
    let bounded = m.prepare_equilibrium(t, maximum).unwrap();
    let roots = bounded.densities(p).unwrap();
    assert_eq!(roots.len(), 1);
    close(roots[0], maximum, 1e-10 * m.rc(), "root at search boundary");
}
