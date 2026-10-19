//! Raw source rounding is compared in the reference oracle's process-main
//! context; physical/domain tests use the normal threaded harness.
use thermopack_rs::{Mbwr, MbwrDensityParameters, Phase};
fn main() {
    let mut previous = String::new();
    let mut prepared = None;
    let mut count = 0;
    let mut paths = [0; 5];
    let mut nonclosing = 0;
    for row in include_str!("fixtures/mbwr_density.csv").lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        let f: Vec<_> = columns[0].split('|').collect();
        let key = format!("{}|{}|{}", f[1], f[2], f[5]);
        if key != previous {
            let a: Vec<f64> = f[5]
                .split('@')
                .nth(1)
                .unwrap()
                .split(',')
                .map(|v| v.parse().unwrap())
                .collect();
            let data = MbwrDensityParameters::new(a[0], a[1], a[2], a[3]).unwrap();
            prepared = Some(
                Mbwr::from_record(f[1])
                    .unwrap()
                    .prepare_baseline_density(f[2].parse().unwrap(), data)
                    .unwrap(),
            );
            previous = key;
        }
        let p: f64 = f[3].parse().unwrap();
        let result = prepared
            .as_ref()
            .unwrap()
            .density(
                p,
                if f[4] == "L" {
                    Phase::Liquid
                } else {
                    Phase::Vapor
                },
            )
            .unwrap();
        let actual = [
            result.density,
            if result.phase == Phase::Liquid {
                1.
            } else {
                2.
            },
            result.pressure,
            result.pressure_density_derivative,
        ];
        let atol = [1e-12, 0., 1e-5, 1e-8];
        for (i, a) in actual.into_iter().enumerate() {
            let b: f64 = columns[i + 1].parse().unwrap();
            assert!(
                (a - b).abs() <= atol[i] + (if i == 1 { 0. } else { 2e-9 }) * a.abs().max(b.abs()),
                "{} field {i}: {a} != {b}",
                columns[0]
            );
        }
        assert_eq!(
            result.pressure_converged,
            result.relative_pressure_error.abs() <= 1e-6
        );
        assert_eq!(
            result.mechanically_stable,
            result.density > 0. && result.pressure_density_derivative >= 0.
        );
        if !result.pressure_converged {
            nonclosing += 1;
        }
        paths[result.path as usize] += 1;
        count += 1;
    }
    assert_eq!(count, 8839);
    assert!(paths.iter().all(|&n| n > 0));
    assert!(nonclosing > 400);
    // The one non-finite original call is retained in fixture metadata. This
    // process-main comparison asserts a native error rather than a false state.
    let model = Mbwr::from_record("mbwrcx58").unwrap();
    let data = MbwrDensityParameters::new(
        190.76990000000714,
        4626986.189794588,
        10097.864021341442,
        0.01131,
    )
    .unwrap();
    assert!(
        model
            .baseline_density(133.538930000005, 0.01, Phase::Liquid, data)
            .is_err()
    );
    println!(
        "{count} raw MBWR density states verified, paths {paths:?}, {nonclosing} raw unclosed returns retained"
    );
}
