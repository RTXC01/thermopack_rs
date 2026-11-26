use std::{hint::black_box, time::Instant};
use thermopack_rs::{
    Phase,
    hardsphere_wca::{cavity_distribution_function, cavity_integral_lj_fres_tvn},
    ljs::{LjUfPotential, LjsBh, LjsUf, LjsUv},
    ljs_wca::{LjsWca, wca_term},
};

fn bench_ns<F: FnMut(usize) -> f64>(name: &str, rounds: usize, iters: usize, mut f: F) {
    // Discarded warm-up
    let mut acc = 0.0;
    for i in 0..(iters / 10).max(10) {
        acc += f(i);
    }
    black_box(acc);

    let mut samples = Vec::with_capacity(rounds);
    let mut checksum = 0.0;
    for _ in 0..rounds {
        let mut s = 0.0;
        let start = Instant::now();
        for i in 0..iters {
            s += f(i);
        }
        let ns = start.elapsed().as_nanos() as f64 / iters as f64;
        checksum = black_box(s);
        samples.push(ns);
    }
    samples.sort_by(f64::total_cmp);
    let median = samples[rounds / 2];
    println!(
        "{name}: median={median:.1} ns (min={:.1}, max={:.1}) checksum={checksum:.6e}",
        samples[0],
        samples[rounds - 1]
    );
}

fn main() {
    let bh = LjsBh::new();
    let wca = LjsWca::new();
    let uv = LjsUv::new();
    let uf = LjsUf::new(LjUfPotential::LennardJonesSplined);
    let rounds = 5;

    // --- Matched baseline workloads (identical to cb6efa1 baseline) ---
    bench_ns("bh_helmholtz_deriv", rounds, 20_000, |i| {
        let v = 1.0e-4 + (i % 100) as f64 * 1.0e-5;
        let (a, at) = black_box(&bh)
            .residual_helmholtz_energy_derivative(black_box(150.0), black_box(v), black_box(1.0))
            .unwrap();
        a + at
    });

    bench_ns("bh_pressure", rounds, 20_000, |i| {
        let v = 1.0e-4 + (i % 100) as f64 * 1.0e-5;
        let (p, z) = black_box(&bh)
            .pressure(black_box(150.0), black_box(v), black_box(1.0))
            .unwrap();
        p * 1.0e-6 + z
    });

    bench_ns("bh_state", rounds, 20_000, |i| {
        let v = 1.0e-4 + (i % 100) as f64 * 1.0e-5;
        let st = black_box(&bh)
            .state(black_box(150.0), black_box(v), black_box(1.0))
            .unwrap();
        st.pressure * 1.0e-6 + st.h_total * 1.0e-3 + st.s_total
    });

    bench_ns("bh_volume_vapor", rounds, 100, |i| {
        let p = 5.0e5 + (i % 20) as f64 * 1.0e5;
        black_box(&bh)
            .volume(black_box(200.0), black_box(p), black_box(Phase::Vapor))
            .unwrap()
    });

    bench_ns("wca_term_a2", rounds, 200_000, |i| {
        let rho = 0.1 + (i % 100) as f64 * 0.008;
        let t = wca_term(2, black_box(124.0), black_box(0.05), black_box(rho)).unwrap();
        t.0 + t.1
    });

    bench_ns("cavity_dist_fn", rounds, 200_000, |i| {
        let eta = 0.05 + (i % 100) as f64 * 0.005;
        let c = cavity_distribution_function(
            black_box(eta),
            black_box(1.05),
            black_box(1.0e-4),
            black_box(2.0e-7),
        );
        c.value + c.packing_fraction
    });

    // --- New entry points ---
    let bh_prep = bh.prepare(150.0).unwrap();
    bench_ns("bh_prepared_tvn", rounds, 200_000, |i| {
        let v = 1.0e-4 + (i % 100) as f64 * 1.0e-5;
        let f = black_box(&bh_prep)
            .residual_helmholtz_tvn(black_box(v), black_box(1.0))
            .unwrap();
        f.value + f.dtt + f.dvv
    });

    bench_ns("bh_volume_liquid", rounds, 2_000, |i| {
        let p = 1.0e6 + (i % 20) as f64 * 1.0e5;
        black_box(&bh)
            .volume(black_box(85.0), black_box(p), black_box(Phase::Liquid))
            .unwrap()
    });

    bench_ns("bh_baseline_volume", rounds, 2_000, |i| {
        let p = 1.0e6 + (i % 20) as f64 * 1.0e5;
        black_box(&bh)
            .baseline_volume_iteration(black_box(85.0), black_box(p), black_box(Phase::Liquid))
            .unwrap()
            .volume
    });

    bench_ns("cavity_integral_tvn", rounds, 5_000, |i| {
        let v = 8.0e-5 + (i % 100) as f64 * 1.0e-6;
        let f = cavity_integral_lj_fres_tvn(
            black_box(3.42e-10),
            black_box(124.0),
            black_box(120.0),
            black_box(v),
            black_box(1.2),
        )
        .unwrap();
        f.value + f.dtt + f.dvv
    });

    let wca_prep = wca.prepare(120.0).unwrap();
    bench_ns("wca_tvn_fresh", rounds, 5_000, |i| {
        let v = 8.0e-5 + (i % 100) as f64 * 1.0e-6;
        let f = black_box(&wca)
            .residual_helmholtz_tvn(black_box(120.0), black_box(v), black_box(1.2))
            .unwrap();
        f.value + f.dtt + f.dvv
    });

    bench_ns("wca_tvn_prepared", rounds, 10_000, |i| {
        let v = 8.0e-5 + (i % 100) as f64 * 1.0e-6;
        let f = black_box(&wca_prep)
            .residual_helmholtz_tvn(black_box(v), black_box(1.2))
            .unwrap();
        f.value + f.dtt + f.dvv
    });

    let wca_prep_liq = wca.prepare(85.0).unwrap();
    bench_ns("wca_volume_liquid_prep", rounds, 200, |i| {
        let p = 1.5e6 + (i % 20) as f64 * 1.0e5;
        black_box(&wca_prep_liq)
            .volume(black_box(p), black_box(Phase::Liquid))
            .unwrap()
    });

    let uv_prep = uv.prepare(120.0).unwrap();
    bench_ns("uv_tvn_prepared", rounds, 10_000, |i| {
        let v = 8.0e-5 + (i % 100) as f64 * 1.0e-6;
        let f = black_box(&uv_prep)
            .residual_helmholtz_tvn(black_box(v), black_box(1.2))
            .unwrap();
        f.value + f.dtt + f.dvv
    });

    let uv_prep_vap = uv.prepare(130.0).unwrap();
    bench_ns("uv_volume_vapor_prep", rounds, 200, |i| {
        let p = 5.0e5 + (i % 20) as f64 * 1.0e5;
        black_box(&uv_prep_vap)
            .volume(black_box(p), black_box(Phase::Vapor))
            .unwrap()
    });

    let uf_prep = uf.prepare(120.0).unwrap();
    bench_ns("uf_tvn_prepared", rounds, 10_000, |i| {
        let v = 8.0e-5 + (i % 100) as f64 * 1.0e-6;
        let f = black_box(&uf_prep)
            .residual_helmholtz_tvn(black_box(v), black_box(1.2))
            .unwrap();
        f.value + f.dtt + f.dvv
    });
}
