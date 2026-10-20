//! Measures the cost of the transcendental functions of the native math
//! runtime. The reference build uses the MinGW universal C runtime while the
//! default Windows toolchain of this crate uses the MinGW GNU C runtime, and
//! the two differ in the speed of the elementary functions.
use std::hint::black_box;
use std::time::Instant;

fn main() {
    let iterations = 1_000_000;
    let mut checksum = 0.0;

    let start = Instant::now();
    for i in 0..iterations {
        let x = black_box(150.0 + (i % 1000) as f64 * 0.7);
        checksum += black_box((-2.0 * x / 300.0).exp());
    }
    println!(
        "exp,{:.6},{checksum:.17e}",
        start.elapsed().as_nanos() as f64 / iterations as f64
    );

    let start = Instant::now();
    for i in 0..iterations {
        let x = black_box(150.0 + (i % 1000) as f64 * 0.7);
        checksum += black_box(x.ln());
    }
    println!(
        "log,{:.6},{checksum:.17e}",
        start.elapsed().as_nanos() as f64 / iterations as f64
    );

    let start = Instant::now();
    for i in 0..iterations {
        let x = black_box(150.0 + (i % 1000) as f64 * 0.7);
        checksum += black_box((x / 300.0).sinh());
    }
    println!(
        "sinh,{:.6},{checksum:.17e}",
        start.elapsed().as_nanos() as f64 / iterations as f64
    );

    let start = Instant::now();
    for i in 0..iterations {
        let x = black_box(150.0 + (i % 1000) as f64 * 0.7);
        checksum += black_box((-2.0 * x / 300.0).exp().ln_1p());
    }
    println!(
        "log1p,{:.6},{checksum:.17e}",
        start.elapsed().as_nanos() as f64 / iterations as f64
    );
}
