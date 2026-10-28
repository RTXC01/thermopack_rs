use std::{env, hint::black_box, time::Instant};
use thermopack_rs::{PcSaft, PcSaftKind};

fn main() {
    let a: Vec<String> = env::args().collect();
    // variant, components, temperature, volume, moles, mode, warmup, iterations
    assert_eq!(a.len(), 9);
    let kind = if a[1] == "simplified" {
        PcSaftKind::SPcSaft
    } else {
        PcSaftKind::PcSaft
    };
    let components: Vec<&str> = a[2].split(',').collect();
    let t: f64 = a[3].parse().unwrap();
    let v: f64 = a[4].parse().unwrap();
    let n: Vec<f64> = a[5].split(',').map(|x| x.parse().unwrap()).collect();
    let sum: f64 = n.iter().sum();
    let x: Vec<f64> = n.iter().map(|v| v / sum).collect();
    let rho = sum / v;
    let warmup: usize = a[7].parse().unwrap();
    let iterations: usize = a[8].parse().unwrap();
    assert!(components.len() > 1 && iterations > 0);
    let model = PcSaft::new_with_options(kind, &components, "Default").unwrap();
    let prepared = model.prepare(t, &x).unwrap();
    let evaluate = || {
        if a[6] == "fresh" {
            model
                .prepare(black_box(t), black_box(&x))
                .unwrap()
                .fugacity_coefficients(black_box(rho), 1.0)
        } else {
            prepared.fugacity_coefficients(black_box(rho), 1.0)
        }
    };
    let vector = evaluate();
    for _ in 0..warmup {
        black_box(evaluate());
    }
    let mut checksum = 0.0;
    let start = Instant::now();
    for _ in 0..iterations {
        for (i, value) in evaluate().iter().enumerate() {
            checksum += black_box(*value) * (i + 1) as f64;
        }
    }
    let seconds = start.elapsed().as_secs_f64();
    println!(
        "RESULT={{\"seconds\":{seconds:.17e},\"checksum\":{checksum:.17e},\"vector\":[{}]}}",
        vector
            .iter()
            .map(|v| format!("{v:.17e}"))
            .collect::<Vec<_>>()
            .join(",")
    );
}
