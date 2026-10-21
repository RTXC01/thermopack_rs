use std::{hint::black_box, time::Instant};
use thermopack_rs::{Cubic, CubicKind, MixtureSaturation};
fn fields(s: &MixtureSaturation) -> [f64; 8] {
    [
        s.temperature,
        s.pressure,
        s.liquid_volume,
        s.vapor_volume,
        s.liquid[0],
        s.liquid[1],
        s.vapor[0],
        s.vapor[1],
    ]
}
fn main() {
    let check = std::env::args().any(|s| s == "--check");
    let iterations = if cfg!(debug_assertions) { 10 } else { 1000 };
    for (profile, kind) in [CubicKind::PengRobinson, CubicKind::SoaveRedlichKwong]
        .into_iter()
        .enumerate()
    {
        let model = Cubic::new(kind, &["C1", "C2"]).unwrap();
        let z = [0.5, 0.5];
        for mode in ["BP", "DP", "BT", "DT"] {
            let solve = |i: usize| {
                let t = black_box(160. + 0.4 * (i % 100) as f64);
                let p = black_box(1e5 + 3e4 * (i % 100) as f64);
                let model = black_box(&model);
                match mode {
                    "BP" => model.bubble_pressure(t, black_box(&z)),
                    "DP" => model.dew_pressure(t, black_box(&z)),
                    "BT" => model.bubble_temperature(p, black_box(&z)),
                    _ => model.dew_temperature(p, black_box(&z)),
                }
                .unwrap()
            };
            if check {
                for i in 0..100 {
                    println!(
                        "check_{}_{mode}_{i},{}",
                        profile + 1,
                        fields(&solve(i))
                            .iter()
                            .map(|v| format!("{v:.17e}"))
                            .collect::<Vec<_>>()
                            .join(",")
                    );
                }
                continue;
            }
            let mut sums = [0.; 8];
            let start = Instant::now();
            for i in 0..iterations {
                let s = black_box(fields(&solve(i)));
                for (sum, v) in sums.iter_mut().zip(s) {
                    *sum += v;
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
}
