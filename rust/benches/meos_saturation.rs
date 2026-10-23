use std::{hint::black_box, time::Instant};
use thermopack_rs::{LennardJonesKind, Multiparameter};

fn main() {
    let check = std::env::args().any(|a| a == "--check");
    let iterations = if cfg!(debug_assertions) { 3 } else { 1000 };
    for (profile, config) in [
        "NIST_MEOS:C3",
        "NIST_MEOS:R134A",
        "LJ_MEOS:KR",
        "LJ_MEOS:AR",
        "LJTS_MEOS:KR",
    ]
    .into_iter()
    .enumerate()
    {
        let (eos, component) = config.split_once(':').unwrap();
        let model = match eos {
            "LJ_MEOS" => Multiparameter::lennard_jones(component, LennardJonesKind::Full),
            "LJTS_MEOS" => {
                Multiparameter::lennard_jones(component, LennardJonesKind::TruncatedShifted)
            }
            _ => Multiparameter::new(component),
        }
        .unwrap();
        for mode in ["SAT_P", "SAT_T"] {
            let inputs: Vec<f64> = include_str!("../../tests/fixtures/meos_saturation.csv")
                .lines()
                .skip(1)
                .filter_map(|row| {
                    let case: Vec<_> = row.split(';').next().unwrap().split('|').collect();
                    (case[1] == config && case[5] == mode).then(|| case[2].parse().unwrap())
                })
                .collect();
            assert!(!inputs.is_empty());
            let solve = |i: usize| {
                let value = black_box(inputs[i % inputs.len()]);
                let q = if mode == "SAT_P" {
                    black_box(&model).saturation_pressure(value)
                } else {
                    black_box(&model).saturation_temperature(value)
                }
                .unwrap();
                [
                    q.temperature,
                    q.pressure,
                    q.liquid_volume,
                    q.vapor_volume,
                    q.dpdt,
                ]
            };
            if check {
                for i in 0..inputs.len() {
                    println!("check_{}_{mode}_{i},{}", profile + 1, encoded(solve(i)));
                }
            } else {
                let start = Instant::now();
                let mut sums = [0.; 5];
                for i in 0..iterations {
                    for (sum, value) in sums.iter_mut().zip(black_box(solve(i))) {
                        *sum += value;
                    }
                }
                println!(
                    "native_{}_{mode},{:.6},{}",
                    profile + 1,
                    start.elapsed().as_nanos() as f64 / f64::from(iterations as u32),
                    encoded(sums)
                );
            }
        }
    }
}

fn encoded(values: [f64; 5]) -> String {
    values.map(|v| format!("{v:.17e}")).join(",")
}
