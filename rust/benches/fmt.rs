use std::{hint::black_box, time::Instant};
use thermopack_rs::fmt::{FmtDensity, FmtModel};

fn fields(s: FmtDensity) -> [f64; 43] {
    let mut v = [0.; 43];
    v[0] = s.value;
    v[1..7].copy_from_slice(&s.gradient);
    for i in 0..6 {
        v[7 + i * 6..13 + i * 6].copy_from_slice(&s.hessian[i]);
    }
    v
}
fn main() {
    let input: [[f64; 6]; 100] = std::array::from_fn(|i| {
        let s = 1. + 0.01 * i as f64;
        [
            0.1 * s,
            0.2 * s,
            0.3 * s,
            0.05 + 0.008 * i as f64,
            0.02 * s,
            -0.04 * s,
        ]
    });
    let check = std::env::args().any(|s| s == "--check");
    let batches = if cfg!(debug_assertions) { 1 } else { 1000 };
    for model in [
        FmtModel::Rosenfeld,
        FmtModel::WhiteBear,
        FmtModel::WhiteBearII,
    ] {
        let mut output = [model.energy_density(input[0]).unwrap(); 100];
        for mode in ["analytic", "baseline"] {
            let mut evaluate = || {
                let model = black_box(model);
                if mode == "analytic" {
                    model
                        .energy_density_grid(black_box(&input), black_box(&mut output))
                        .unwrap();
                } else {
                    model
                        .baseline_energy_density_grid(black_box(&input), black_box(&mut output))
                        .unwrap();
                }
            };
            if check {
                evaluate();
                for (i, s) in output.into_iter().enumerate() {
                    println!(
                        "check_{mode}_{}_{i},{}",
                        model as i32,
                        fields(s)
                            .iter()
                            .map(|v| format!("{v:.17e}"))
                            .collect::<Vec<_>>()
                            .join(",")
                    );
                }
                continue;
            }
            let start = Instant::now();
            let mut checksum = 0.;
            for _ in 0..batches {
                // Reborrow output per call; all 43 fields contribute to the checksum.
                let model = black_box(model);
                if mode == "analytic" {
                    model
                        .energy_density_grid(black_box(&input), black_box(&mut output))
                        .unwrap();
                } else {
                    model
                        .baseline_energy_density_grid(black_box(&input), black_box(&mut output))
                        .unwrap();
                }
                for s in black_box(&output) {
                    checksum += fields(*s).iter().sum::<f64>();
                }
            }
            println!(
                "{mode}_{},{:.6},{checksum:.17e}",
                model as i32,
                start.elapsed().as_nanos() as f64 / (100 * batches) as f64
            );
        }
    }
}
