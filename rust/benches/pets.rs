use std::{hint::black_box, time::Instant};
use thermopack_rs::{Pets, PureHelmholtz};
fn values(f: PureHelmholtz) -> [f64; 10] {
    [
        f.value, f.dt, f.dv, f.dn, f.dtt, f.dtv, f.dtn, f.dvv, f.dvn, f.dnn,
    ]
}
fn main() {
    let model = Pets::new().unwrap();
    let prepared = model.prepare(300.).unwrap();
    if std::env::args().any(|s| s == "--check") {
        for i in 0..1000 {
            let v = 5e-5 + i as f64 * 1e-6;
            let f = prepared.residual_helmholtz(v, 1.).unwrap();
            println!(
                "check_{i},{}",
                values(f)
                    .iter()
                    .map(|v| format!("{v:.17e}"))
                    .collect::<Vec<_>>()
                    .join(",")
            );
        }
        return;
    }
    let iterations = if cfg!(debug_assertions) {
        2000
    } else {
        1_000_000
    };
    for fresh in [false, true] {
        let mut checksum = [0.; 10];
        let start = Instant::now();
        for i in 0..iterations {
            let v = black_box(5e-5 + (i % 1000) as f64 * 1e-6);
            let f = if fresh {
                black_box(&model).residual_helmholtz(black_box(300.), v, black_box(1.))
            } else {
                black_box(&prepared).residual_helmholtz(v, black_box(1.))
            }
            .unwrap();
            for (sum, value) in checksum.iter_mut().zip(values(f)) {
                *sum += value;
            }
        }
        println!(
            "{},{:.6},{}",
            if fresh {
                "with_temperature_preparation"
            } else {
                "prepared"
            },
            start.elapsed().as_nanos() as f64 / iterations as f64,
            checksum
                .iter()
                .map(|v| format!("{v:.17e}"))
                .collect::<Vec<_>>()
                .join(",")
        );
    }
}
