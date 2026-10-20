use std::{hint::black_box, time::Instant};
use thermopack_rs::{HelmholtzJet, LennardJonesKind, Multiparameter, MultiparameterState};
fn jet(j: HelmholtzJet) -> [f64; 6] {
    [
        j.alpha,
        j.d_tau,
        j.d_tautau,
        j.d_delta,
        j.d_deltadelta,
        j.d_deltatau,
    ]
}
fn fields(s: MultiparameterState) -> [f64; 24] {
    let mut v = [0.; 24];
    v[..6].copy_from_slice(&jet(s.ideal));
    v[6..12].copy_from_slice(&jet(s.residual));
    v[12..].copy_from_slice(&[
        s.pressure.value,
        s.pressure.dt,
        s.pressure.dv,
        s.internal_energy,
        s.entropy,
        s.enthalpy,
        s.helmholtz_energy,
        s.gibbs_energy,
        s.heat_capacity_v,
        s.heat_capacity_p,
        s.ln_fugacity.unwrap(),
        s.baseline_speed_of_sound.unwrap(),
    ]);
    v
}
fn main() {
    let models = [
        Multiparameter::new("C3").unwrap(),
        Multiparameter::new("R134A").unwrap(),
        Multiparameter::new("N-H2").unwrap(),
        Multiparameter::lennard_jones("AR", LennardJonesKind::Full).unwrap(),
        Multiparameter::lennard_jones("AR", LennardJonesKind::TruncatedShifted).unwrap(),
        Multiparameter::reduced_lennard_jones_truncated_shifted(),
    ];
    let check = std::env::args().any(|s| s == "--check");
    let iterations = if cfg!(debug_assertions) { 100 } else { 100000 };
    for (profile, eos) in models.iter().enumerate() {
        let t = 1.4 * eos.tc();
        let prepared = eos.prepare(t).unwrap();
        for mode in ["fresh", "prepared"] {
            let evaluate = |v| {
                if mode == "fresh" {
                    black_box(eos)
                        .prepare(black_box(t))
                        .unwrap()
                        .state(v)
                        .unwrap()
                } else {
                    black_box(&prepared).state(v).unwrap()
                }
            };
            if check {
                for i in 0..100 {
                    let v = 1. / (eos.rc() * (0.05 + 0.015 * i as f64));
                    println!(
                        "check_{mode}_{}_{i},{}",
                        profile + 1,
                        fields(evaluate(v))
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
            for i in 0..iterations {
                let v = black_box(1. / (eos.rc() * (0.05 + 0.015 * (i % 100) as f64)));
                checksum += black_box(fields(evaluate(v))).iter().sum::<f64>();
            }
            println!(
                "{mode}_{},{:.6},{checksum:.17e}",
                profile + 1,
                start.elapsed().as_nanos() as f64 / iterations as f64
            );
        }
    }
}
