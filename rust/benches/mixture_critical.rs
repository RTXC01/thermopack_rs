use std::{hint::black_box, time::Instant};
use thermopack_rs::{Cubic, CubicKind, critical::critical_point};

fn main() {
    let check = std::env::args().any(|a| a == "--check");
    for (profile, (kind, names)) in [
        (CubicKind::PengRobinson, vec!["C1", "C2"]),
        (CubicKind::SoaveRedlichKwong, vec!["C1", "C2"]),
        (CubicKind::PengRobinson, vec!["C1", "C2", "C3"]),
    ]
    .into_iter()
    .enumerate()
    {
        let model = Cubic::new(kind, &names).unwrap();
        let solve = |i: usize| {
            let x = black_box((i % 25) as f64 / 24.0);
            let z = if profile == 2 {
                vec![0.1 + 0.6 * x, 0.2, 0.7 - 0.6 * x]
            } else {
                vec![0.1 + 0.8 * x, 0.9 - 0.8 * x]
            };
            let p = critical_point(black_box(&model), black_box(&z)).unwrap();
            [p.temperature, p.volume, p.pressure]
        };
        if check {
            for i in 0..25 {
                let p = solve(i);
                println!(
                    "check_{}_{i},{:.17e},{:.17e},{:.17e}",
                    profile + 1,
                    p[0],
                    p[1],
                    p[2]
                );
            }
            continue;
        }
        let iterations = if cfg!(debug_assertions) { 25 } else { 200 };
        let mut sums = [0.0; 3];
        let start = Instant::now();
        for i in 0..iterations {
            let p = black_box(solve(i));
            for j in 0..3 {
                sums[j] += p[j];
            }
        }
        println!(
            "native_{},{:.6},{:.17e},{:.17e},{:.17e}",
            profile + 1,
            start.elapsed().as_nanos() as f64 / iterations as f64,
            sums[0],
            sums[1],
            sums[2]
        );
    }
}
