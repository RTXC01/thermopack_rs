use std::hint::black_box;
use std::time::Instant;
use thermopack_rs::{Cubic, CubicKind, Phase};

fn main() {
    let model = Cubic::new(CubicKind::PengRobinson, &["CO2", "C1", "N2"]).unwrap();
    let state = model.prepare(300.0, &[0.7, 0.2, 0.1]).unwrap();
    let iterations = 1_000_000;
    for prepare in [false, true] {
        let mut checksum = [0.0; 3];
        let start = Instant::now();
        for i in 0..iterations {
            let prepared;
            let current = if prepare {
                prepared = model
                    .prepare(black_box(300.0), black_box(&[0.7, 0.2, 0.1]))
                    .unwrap();
                &prepared
            } else {
                &state
            };
            let p = current
                .pressure(black_box(0.0001 + (i % 1000) as f64 * 1e-6))
                .unwrap();
            checksum[0] += p.value;
            checksum[1] += p.dt;
            checksum[2] += p.dv;
        }
        println!(
            "pressure_{},{:.6},{:.17e},{:.17e},{:.17e}",
            if prepare { "with_mixing" } else { "prepared" },
            start.elapsed().as_nanos() as f64 / iterations as f64,
            checksum[0],
            checksum[1],
            checksum[2]
        );
    }
    for prepare in [false, true] {
        let mut checksum = 0.0;
        let start = Instant::now();
        for i in 0..iterations {
            let prepared;
            let current = if prepare {
                prepared = model
                    .prepare(black_box(300.0), black_box(&[0.7, 0.2, 0.1]))
                    .unwrap();
                &prepared
            } else {
                &state
            };
            checksum += current
                .volume(black_box(1e6 + (i % 1000) as f64 * 1000.0), Phase::Vapor)
                .unwrap();
        }
        println!(
            "volume_{},{:.6},{checksum:.17e}",
            if prepare { "with_mixing" } else { "prepared" },
            start.elapsed().as_nanos() as f64 / iterations as f64
        );
    }
}
