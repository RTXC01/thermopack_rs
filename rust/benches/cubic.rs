use std::hint::black_box;
use std::time::Instant;
use thermopack_rs::{Cubic, CubicKind, Phase};

fn main() {
    let model = Cubic::new(CubicKind::PengRobinson, &["CO2", "C1", "N2"]).unwrap();
    let state = model.prepare(300.0, &[0.7, 0.2, 0.1]).unwrap();
    let iterations = 1_000_000;
    let start = Instant::now();
    for i in 0..iterations {
        black_box(
            state
                .pressure(black_box(0.0001 + (i % 1000) as f64 * 1e-6))
                .unwrap(),
        );
    }
    println!(
        "prepared pressure: {:.1} ns/state",
        start.elapsed().as_nanos() as f64 / iterations as f64
    );
    let start = Instant::now();
    for i in 0..iterations {
        black_box(
            state
                .volume(black_box(1e6 + (i % 1000) as f64 * 1000.0), Phase::Vapor)
                .unwrap(),
        );
    }
    println!(
        "prepared volume roots: {:.1} ns/state",
        start.elapsed().as_nanos() as f64 / iterations as f64
    );
}
