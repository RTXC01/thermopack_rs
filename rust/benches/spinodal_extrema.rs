use std::{hint::black_box, time::Instant};
use thermopack_rs::spinodal_properties::{
    BaselineSpinodalExtrema, SpinodalProperty, baseline_pure_spinodal_extrema,
};
use thermopack_rs::{Cubic, CubicKind};
fn values(e: BaselineSpinodalExtrema) -> [f64; 6] {
    assert!(!e.ordering_ambiguous);
    [
        e.minimum.temperature,
        e.minimum.volume,
        e.minimum.pressure,
        e.maximum.temperature,
        e.maximum.volume,
        e.maximum.pressure,
    ]
}
fn main() {
    let check = std::env::args().any(|s| s == "--check");
    let properties = [
        SpinodalProperty::Entropy,
        SpinodalProperty::Enthalpy,
        SpinodalProperty::InternalEnergy,
    ];
    for (profile, (kind, name)) in [
        (CubicKind::PengRobinson, "C1"),
        (CubicKind::SoaveRedlichKwong, "N2"),
    ]
    .into_iter()
    .enumerate()
    {
        let model = Cubic::new(kind, &[name]).unwrap();
        let start = Instant::now();
        let cached = model
            .prepare_pure_spinodal_search(&[1.], 1e6, None)
            .unwrap();
        let setup = start.elapsed().as_nanos();
        for mode in ["fresh", "prepared"] {
            let solve = |i: usize| {
                values(if mode == "fresh" {
                    baseline_pure_spinodal_extrema(
                        black_box(&model),
                        &[1.],
                        black_box(properties[i % 3]),
                        1e6,
                    )
                    .unwrap()
                } else {
                    cached.extrema(black_box(properties[i % 3])).unwrap()
                })
            };
            if check {
                for i in 0..3 {
                    let q = solve(i);
                    println!(
                        "check_{}_{mode}_{i},{}",
                        profile + 1,
                        q.iter()
                            .map(|v| format!("{v:.17e}"))
                            .collect::<Vec<_>>()
                            .join(",")
                    );
                }
            } else {
                let iterations = if cfg!(debug_assertions) { 3 } else { 100 };
                let mut sums = [0.; 6];
                let start = Instant::now();
                for i in 0..iterations {
                    let q = black_box(solve(i));
                    for j in 0..6 {
                        sums[j] += q[j];
                    }
                }
                println!(
                    "native_{}_{mode},{:.6},{}",
                    profile + 1,
                    start.elapsed().as_nanos() as f64 / iterations as f64,
                    sums.iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                );
            }
        }
        if !check {
            println!("setup_{},{setup}", profile + 1);
        }
    }
}
