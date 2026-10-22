use std::{hint::black_box, time::Instant};
use thermopack_rs::{Cubic, CubicKind, Phase};

fn main() {
    let check = std::env::args().any(|s| s == "--check");
    for profile in 1..=4 {
        let names: &[&str] = match profile {
            1 | 4 => &["C1", "C2"],
            2 => &["CO2", "C3"],
            _ => &["C1"],
        };
        let model = Cubic::new(
            if profile == 2 {
                CubicKind::SoaveRedlichKwong
            } else {
                CubicKind::PengRobinson
            },
            names,
        )
        .unwrap();
        let z = if names.len() == 1 {
            vec![1.]
        } else {
            vec![0.5, 0.5]
        };
        let queries = if profile <= 2 { 25 } else { 5 };
        let mut specs = Vec::new();
        for i in 0..queries {
            let t = match profile {
                1 => 250. + i as f64,
                2 => 350. + i as f64,
                3 => 180. + 0.5 * i as f64,
                _ => 250. + 0.5 * i as f64,
            };
            let state = model.prepare(t, &z).unwrap();
            let p0 = if profile == 3 { 6e6 } else { 1e7 };
            let v = if profile <= 2 {
                2e-3 + 2e-5 * i as f64
            } else {
                state.volume(p0, Phase::Liquid).unwrap()
            };
            let c = state.caloric_tv(v).unwrap();
            specs.push((t, v, if profile <= 2 { c.pressure } else { p0 }, c.entropy));
        }
        let solve = |i: usize| {
            let (t, v, p, s) = black_box(specs[i % queries]);
            if profile <= 2 {
                let q = model.homogeneous_ps(p, s, &z, 1.02 * t, 0.98 * v).unwrap();
                vec![[
                    q.temperature,
                    q.volume,
                    q.pressure + q.pressure_residual,
                    q.entropy,
                ]]
            } else {
                model
                    .baseline_metastable_isentrope(&z, p, s, 1e5, 12)
                    .unwrap()
                    .points
                    .into_iter()
                    .map(|q| [q.temperature, q.volume, q.pressure, q.entropy])
                    .collect()
            }
        };
        if check {
            for i in 0..queries {
                for (j, q) in solve(i).into_iter().enumerate() {
                    println!(
                        "check_{profile}_{i}_{j},{:.17e},{:.17e},{:.17e},{:.17e}",
                        q[0], q[1], q[2], q[3]
                    );
                }
            }
            continue;
        }
        let iterations = if cfg!(debug_assertions) {
            queries
        } else if profile <= 2 {
            10000
        } else {
            20
        };
        let mut sums = [0.; 4];
        let start = Instant::now();
        for i in 0..iterations {
            for q in black_box(solve(i)) {
                for j in 0..4 {
                    sums[j] += q[j];
                }
            }
        }
        println!(
            "native_{profile},{:.6},{:.17e},{:.17e},{:.17e},{:.17e}",
            start.elapsed().as_nanos() as f64 / iterations as f64,
            sums[0],
            sums[1],
            sums[2],
            sums[3]
        );
    }
}
