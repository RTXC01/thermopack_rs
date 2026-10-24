use std::{hint::black_box, time::Instant};
use thermopack_rs::{LennardJonesKind, Multiparameter, MultiparameterEquilibrium, Phase};

fn fields(eq: &MultiparameterEquilibrium<'_>, p: f64, buffer: &mut [f64]) -> [f64; 5] {
    let n = eq.densities_into(black_box(p), buffer).unwrap();
    assert!(n > 0);
    [
        n as f64,
        buffer[..n].iter().sum(),
        buffer[0],
        buffer[n - 1],
        1. / eq.volume(p, Phase::MinimumGibbs).unwrap(),
    ]
}
fn encoded(q: [f64; 5]) -> String {
    q.map(|v| format!("{v:.17e}")).join(",")
}
fn main() {
    let check = std::env::args().any(|s| s == "--check");
    for (profile, (config, tr)) in [
        ("NIST_MEOS:C3", 0.95),
        ("NIST_MEOS:R134A", 1.05),
        ("NIST_MEOS:N-H2", 0.95),
        ("LJ_MEOS:AR", 0.95),
        ("LJTS_MEOS:KR", 0.95),
    ]
    .into_iter()
    .enumerate()
    {
        let (e, c) = config.split_once(':').unwrap();
        let m = match e {
            "LJ_MEOS" => Multiparameter::lennard_jones(c, LennardJonesKind::Full),
            "LJTS_MEOS" => Multiparameter::lennard_jones(c, LennardJonesKind::TruncatedShifted),
            _ => Multiparameter::new(c),
        }
        .unwrap();
        let t = tr * m.tc();
        let maximum = 8. * m.rc();
        let inputs: Vec<f64> = include_str!("../../tests/fixtures/meos_pressure_branches.csv")
            .lines()
            .skip(1)
            .filter_map(|line| {
                let case: Vec<_> = line.split(';').next().unwrap().split('|').collect();
                if case[1] != config
                    || (case[2].parse::<f64>().unwrap() / m.tc() - tr).abs() > 1e-12
                {
                    return None;
                }
                Some(case[5].strip_prefix("ROOTS:").unwrap().parse().unwrap())
            })
            .collect();
        assert!(!inputs.is_empty());
        let eq = m.prepare_equilibrium(t, maximum).unwrap();
        let mut buffer = vec![0.; eq.maximum_root_count()];
        if check {
            for (i, &p) in inputs.iter().enumerate() {
                println!(
                    "check_{}_{i},{}",
                    profile + 1,
                    encoded(fields(&eq, p, &mut buffer))
                );
            }
            continue;
        }
        let setup_iterations = if cfg!(debug_assertions) { 1 } else { 10 };
        let start = Instant::now();
        for _ in 0..setup_iterations {
            black_box(
                black_box(&m)
                    .prepare_equilibrium(black_box(t), black_box(maximum))
                    .unwrap(),
            );
        }
        println!(
            "setup_{},{:.6}",
            profile + 1,
            start.elapsed().as_nanos() as f64 / f64::from(setup_iterations)
        );
        for fresh in [true, false] {
            let iterations = if cfg!(debug_assertions) {
                2
            } else if fresh {
                100
            } else {
                1000
            };
            let start = Instant::now();
            let mut sums = [0.; 5];
            for i in 0..iterations {
                let p = black_box(inputs[i % inputs.len()]);
                let q = if fresh {
                    let eq = black_box(&m)
                        .prepare_equilibrium(black_box(t), black_box(maximum))
                        .unwrap();
                    fields(&eq, p, &mut buffer)
                } else {
                    fields(black_box(&eq), p, &mut buffer)
                };
                for (sum, value) in sums.iter_mut().zip(black_box(q)) {
                    *sum += value;
                }
            }
            println!(
                "{}_{},{:.6},{}",
                if fresh { "fresh" } else { "prepared" },
                profile + 1,
                start.elapsed().as_nanos() as f64 / iterations as f64,
                encoded(sums)
            );
        }
    }
}
