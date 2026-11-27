use std::{hint::black_box, time::Instant};
use thermopack_rs::{
    ApparentComposition, Cpa, CpaOptions, Cubic, CubicKind, CubicOptions, HvPairParams, MixingRule,
    Phase, init_umr, init_vtpr,
};

fn bench_ns<F: FnMut(usize) -> f64>(name: &str, rounds: usize, iters: usize, mut f: F) {
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
    let rounds = 5;

    // 1. CPA-SRK VdW (pure H2O and binary H2O+ETOH)
    let cpa_h2o = Cpa::new(&["H2O"]).unwrap();
    let prep_h2o = cpa_h2o.prepare(373.15, &[1.0]).unwrap();
    bench_ns("cpa_h2o_pressure_and_derivs", rounds, 50_000, |i| {
        let rho = 50.0 + (i % 100) as f64 * 10.0;
        let jet = black_box(&prep_h2o)
            .pressure_and_derivatives(black_box(rho))
            .unwrap();
        jet.pressure * 1e-6 + jet.dp_dt * 1e-3
    });
    bench_ns("cpa_h2o_state", rounds, 20_000, |i| {
        let rho = 50.0 + (i % 100) as f64 * 10.0;
        let st = black_box(&prep_h2o).state(black_box(rho)).unwrap();
        st.pressure * 1e-6 + st.ln_phi[0] + st.cp_total
    });

    let cpa_bin = Cpa::new(&["H2O", "ETOH"]).unwrap();
    let prep_bin = cpa_bin.prepare(360.0, &[0.4, 0.6]).unwrap();
    bench_ns("cpa_bin_tvn_prepared", rounds, 20_000, |i| {
        let v = 4.0e-5 + (i % 100) as f64 * 2.0e-6;
        let h = black_box(&prep_bin)
            .residual_helmholtz_tv_unit(black_box(v))
            .unwrap();
        h.f + h.f_tt + h.f_nn[0]
    });
    bench_ns("cpa_bin_state_tp", rounds, 2_000, |i| {
        let p = 1.0e6 + (i % 20) as f64 * 1.0e5;
        let st = black_box(&cpa_bin)
            .state_tp(
                black_box(340.0),
                black_box(p),
                black_box(&[0.55, 0.45]),
                Phase::Liquid,
            )
            .unwrap();
        st.volume * 1e5 + st.ln_fugacity.ln_phi[0] + st.cp
    });

    // 2. CPA-SRK HV2
    let mut cpa_hv2 = Cpa::with_cpa_options(
        &["H2O", "MEOH"],
        CpaOptions {
            mixing_rule: MixingRule::HuronVidal2,
            ..Default::default()
        },
    )
    .unwrap();
    cpa_hv2
        .set_hv_param(
            0,
            1,
            HvPairParams {
                alpha_ij: 0.25,
                alpha_ji: 0.30,
                a_ge_ij: 180.0,
                a_ge_ji: -90.0,
                b_ge_ij: 0.12,
                b_ge_ji: -0.08,
                c_ge_ij: 4.0e-4,
                c_ge_ji: -3.0e-4,
            },
        )
        .unwrap();
    let prep_hv2 = cpa_hv2.prepare(350.0, &[0.5, 0.5]).unwrap();
    bench_ns("cpa_hv2_tvn_prepared", rounds, 20_000, |i| {
        let v = 4.2e-5 + (i % 100) as f64 * 1.0e-6;
        let h = black_box(&prep_hv2)
            .residual_helmholtz_tv_unit(black_box(v))
            .unwrap();
        h.f + h.f_tv + h.f_nn[1]
    });

    // 3. Cubic HVWongSandler, VTPR, UMR
    let mut hvws = Cubic::with_options(
        CubicKind::PengRobinson,
        &["CO2", "C1"],
        CubicOptions {
            mixing_rule: MixingRule::HvWongSandler,
            ..Default::default()
        },
    )
    .unwrap();
    hvws.set_hv_param(
        0,
        1,
        HvPairParams {
            alpha_ij: 0.28,
            alpha_ji: 0.32,
            a_ge_ij: 220.0,
            a_ge_ji: -110.0,
            b_ge_ij: 0.15,
            b_ge_ji: -0.10,
            c_ge_ij: 3.0e-4,
            c_ge_ji: -2.0e-4,
        },
    )
    .unwrap();
    hvws.set_kij(0, 1, 0.06).unwrap();
    hvws.set_lij(0, 1, 0.02).unwrap();
    bench_ns("cubic_hvws_prepare_tvn", rounds, 20_000, |i| {
        let t = 310.0 + (i % 50) as f64 * 0.5;
        let h = black_box(&hvws)
            .residual_helmholtz_tv(black_box(t), black_box(2.5e-4), black_box(&[0.45, 0.55]))
            .unwrap();
        h.f + h.f_tt + h.f_nn[1]
    });

    let vtpr = init_vtpr(&["C1", "C3", "NC4"]).unwrap();
    bench_ns("cubic_vtpr_prepare_tvn", rounds, 20_000, |i| {
        let t = 310.0 + (i % 50) as f64 * 0.5;
        let h = black_box(&vtpr)
            .residual_helmholtz_tv(
                black_box(t),
                black_box(3.0e-4),
                black_box(&[0.4, 0.35, 0.25]),
            )
            .unwrap();
        h.f + h.f_tt + h.f_nn[4]
    });

    let umr = init_umr(&["C1", "C3", "H2O"]).unwrap();
    bench_ns("cubic_umr_prepare_tvn", rounds, 20_000, |i| {
        let t = 350.0 + (i % 50) as f64 * 0.5;
        let h = black_box(&umr)
            .residual_helmholtz_tv(
                black_box(t),
                black_box(4.0e-4),
                black_box(&[0.4, 0.35, 0.25]),
            )
            .unwrap();
        h.f + h.f_tt + h.f_nn[4]
    });

    // 4. Apparent composition transformation
    let mut app = ApparentComposition::with_ncsym(3, 5, 1).unwrap();
    app.set_v_stoich(1, 1, 1.0).unwrap();
    app.set_v_stoich(1, 2, 1.0).unwrap();
    app.set_v_stoich(2, 3, 1.0).unwrap();
    app.set_v_stoich(2, 2, 2.0).unwrap();
    let n = [0.85, 0.10, 0.05];
    let ne = app.apparent_to_real_mole_numbers(&n).unwrap();
    let lnf_r = [-0.20, -0.55, -0.65, -1.10, 0.0];
    let dt_r = [1.0e-3, -2.0e-3, -1.5e-3, -3.2e-3, 0.0];
    let dp_r = [2.5e-8, -3.5e-8, 4.0e-8, -6.0e-8, 0.0];
    let dn_r = [
        0.35, -0.08, 0.04, -0.02, 0.0, -0.08, 0.85, 0.20, 0.10, 0.0, 0.04, 0.20, 1.05, 0.18, 0.0,
        -0.02, 0.10, 0.18, 1.40, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ];
    bench_ns("apparent_tp_lnfug_3x5", rounds, 100_000, |_| {
        let r = black_box(&app)
            .tp_lnfug_apparent(
                black_box(&n),
                black_box(&ne),
                black_box(2.5e5),
                black_box(&lnf_r),
                Some(black_box(&dt_r)),
                Some(black_box(&dp_r)),
                Some(black_box(&dn_r)),
            )
            .unwrap();
        r.ln_fugacity[1] + r.dn.unwrap()[4]
    });
}
