//! Complete LJS (BH, WCA, UV, LJS-UF, LJ-UF) and WCA cavity/hard-sphere
//! exact-source parity, derivative, thermodynamic identity, and volume-solver tests.
//! Generated/validated by tools/rust_port/validate_ljs_complete.py.
//! Copyright thermotools and Contributors. Apache-2.0.

#![allow(clippy::excessive_precision)]

use std::fs;
use thermopack_rs::Phase;
use thermopack_rs::hardsphere_wca::{
    cavity_integral_eta_t, cavity_integral_lj_fres_tvn, hard_sphere_diameter,
    packing_fraction_vdw_no_segments, pure_hard_sphere_eta,
};
use thermopack_rs::helmholtz::PureHelmholtz;
use thermopack_rs::ljs::{
    ALPHA_LJS, LjUfPotential, LjsBh, LjsControl, LjsTpState, LjsTvState, LjsUf, LjsUfControl,
    LjsUv, LjsUvControl, a1_ljs_ex, a2_ljs_ex, a3_ljs_ex, ai_ljs_ex, calc_uf_wca, calc_uf_wca_tvn,
    calc_uv_wca_tvn, uv_a1_u_mult_t,
};
use thermopack_rs::ljs_wca::{LjsWca, LjsWcaControl, wca_term_tvn};

const EXPECTED_FIXTURE_CASES: usize = 72;

fn tvn_slots(h: PureHelmholtz) -> [f64; 10] {
    [
        h.value, h.dt, h.dv, h.dn, h.dtt, h.dtv, h.dtn, h.dvv, h.dvn, h.dnn,
    ]
}

fn tp_slots(tp: &LjsTpState) -> [f64; 20] {
    let t = tp.temperature;
    let p = tp.pressure;
    let r = thermopack_rs::R;
    // Residual TP derivatives follow directly from total minus ideal-gas values
    // (Cp_id = 2.5*R, dh_id/dP = 0, ds_id/dT = 2.5*R/T, ds_id/dP = -R/P).
    let dhr_dt = tp.enthalpy.dt - 2.5 * r;
    let dhr_dp = tp.enthalpy.dp;
    let dsr_dt = tp.entropy.dt - 2.5 * r / t;
    let dsr_dp = tp.entropy.dp + r / p;
    [
        tp.molar_volume,
        tp.compressibility.value,
        tp.compressibility.dt,
        tp.compressibility.dp,
        tp.tv.h_residual,
        dhr_dt,
        dhr_dp,
        tp.enthalpy.value,
        tp.enthalpy.dt,
        tp.enthalpy.dp,
        tp.tv.s_residual,
        dsr_dt,
        dsr_dp,
        tp.entropy.value,
        tp.entropy.dt,
        tp.entropy.dp,
        tp.ln_fugacity_coefficient.value,
        tp.ln_fugacity_coefficient.dt,
        tp.ln_fugacity_coefficient.dp,
        tp.speed_of_sound,
    ]
}

fn assert_close(label: &str, source: f64, native: f64, atol: f64, rtol: f64) {
    assert!(
        source.is_finite() && native.is_finite(),
        "{label}: non-finite value (source={source}, native={native})"
    );
    let diff = (source - native).abs();
    let tol = atol + rtol * source.abs();
    assert!(
        diff <= tol,
        "{label}: source={source:.16e}, native={native:.16e}, diff={diff:.3e}, tol={tol:.3e}"
    );
}

#[test]
fn exact_source_ljs_complete() {
    let text = match std::env::var("LJS_COMPLETE_FIXTURE") {
        Ok(path) => fs::read_to_string(path).expect("read LJS_COMPLETE_FIXTURE"),
        Err(_) => include_str!("fixtures/ljs_complete.csv").to_string(),
    };
    let mut checked = 0usize;
    let mut max_rel_tvn = 0.0_f64;
    let mut max_rel_tp = 0.0_f64;

    for line in text.lines().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split(';').collect();
        assert_eq!(cols.len(), 11, "malformed row: {line}");
        let model = cols[0];
        let variant = cols[1];
        let kind: usize = cols[2].parse().unwrap();
        let phase_flag: i32 = cols[3].parse().unwrap();
        let t: f64 = cols[4].parse().unwrap();
        let v_or_p: f64 = cols[5].parse().unwrap();
        let n1: f64 = cols[6].parse().unwrap();
        let sigma: f64 = cols[7].parse().unwrap();
        let eps: f64 = cols[8].parse().unwrap();
        let tag = cols[9];
        let source: Vec<f64> = cols[10].split(',').map(|s| s.parse().unwrap()).collect();
        let phase = if phase_flag == 1 {
            Phase::Liquid
        } else {
            Phase::Vapor
        };
        let label = format!("{model}/{variant}/kind={kind}/T={t}/v={v_or_p}");

        let native: Vec<f64> = match model {
            "CAVITY" => {
                let h = cavity_integral_lj_fres_tvn(sigma, eps, t, v_or_p, n1).unwrap();
                tvn_slots(h).to_vec()
            }
            "ZETAX2" => {
                let (d1, d1t, d1tt) = hard_sphere_diameter(sigma, eps, t);
                let (d2, d2t, d2tt) = hard_sphere_diameter(3.7e-10, 150.0, t);
                let d12 = 0.5 * (d1 + d2);
                let d12t = 0.5 * (d1t + d2t);
                let d12tt = 0.5 * (d1tt + d2tt);
                let z = packing_fraction_vdw_no_segments(
                    v_or_p,
                    &[n1, 0.7],
                    &[d1, d12, d12, d2],
                    &[d1t, d12t, d12t, d2t],
                    &[d1tt, d12tt, d12tt, d2tt],
                )
                .unwrap();
                vec![
                    z.value,
                    z.temperature,
                    z.temperature_second,
                    z.volume,
                    z.volume_second,
                    z.volume_third,
                    z.temperature_volume,
                    z.volume_temperature_second,
                    z.volume_temperature_third,
                    z.moles[0],
                    z.moles[1],
                    z.temperature_moles[0],
                    z.temperature_moles[1],
                    z.volume_moles[0],
                    z.volume_moles[1],
                    z.volume_temperature_moles[0],
                    z.volume_temperature_moles[1],
                    z.moles_second[0],
                    z.moles_second[1],
                    z.moles_second[2],
                    z.moles_second[3],
                    z.volume_moles_second[0],
                    z.volume_moles_second[1],
                    z.volume_moles_second[2],
                    z.volume_moles_second[3],
                ]
            }
            "UF-KERNEL" => {
                let is_lj = variant == "lj";
                if kind == 7 {
                    let u = calc_uf_wca(v_or_p, t, is_lj);
                    vec![u.0, u.1, u.2, u.3, u.4, u.5]
                } else {
                    let h = calc_uf_wca_tvn(sigma, eps, is_lj, t, v_or_p, n1).unwrap();
                    tvn_slots(h).to_vec()
                }
            }
            "LJS-BH" => {
                let mut m = LjsBh::new();
                if sigma > 0.0 && eps > 0.0 {
                    m.set_parameters(sigma, eps).unwrap();
                }
                match variant {
                    "default" => {}
                    "lafitte" => m.set_control(LjsControl {
                        use_lafitte_a3: true,
                        ..LjsControl::default()
                    }),
                    "no-chi" => m.set_control(LjsControl {
                        enable_chi_correction: false,
                        ..LjsControl::default()
                    }),
                    "bh-a1-only" => m.set_control(LjsControl {
                        enable_hard_sphere: false,
                        enable_a1: true,
                        enable_a2: false,
                        enable_a3: false,
                        ..LjsControl::default()
                    }),
                    "bh-a2-only" => m.set_control(LjsControl {
                        enable_hard_sphere: false,
                        enable_a1: false,
                        enable_a2: true,
                        enable_a3: false,
                        ..LjsControl::default()
                    }),
                    "bh-a3-only" => m.set_control(LjsControl {
                        enable_hard_sphere: false,
                        enable_a1: false,
                        enable_a2: false,
                        enable_a3: true,
                        ..LjsControl::default()
                    }),
                    "bh-hs-only" => m.set_control(LjsControl {
                        enable_hard_sphere: true,
                        enable_a1: false,
                        enable_a2: false,
                        enable_a3: false,
                        ..LjsControl::default()
                    }),
                    "kolafa" => m.set_kolafa_pure_reference(true),
                    other => panic!("unexpected BH variant: {other}"),
                }
                match kind {
                    0 => tvn_slots(m.residual_helmholtz_tvn(t, v_or_p, n1).unwrap()).to_vec(),
                    1 => tvn_slots(m.dispersion_tvn(t, v_or_p, n1).unwrap()).to_vec(),
                    2 => tvn_slots(m.hard_sphere_tvn(t, v_or_p, n1).unwrap()).to_vec(),
                    4 => tp_slots(&m.state_tp(t, v_or_p, phase).unwrap()).to_vec(),
                    5 => {
                        let (d, dt, dtt) = m.diameter(t).unwrap();
                        let (tc, vc) = m.estimate_critical_parameters().unwrap();
                        let (fmt_ok, na_en) = m.fmt_compatibility();
                        let (pot, force) = m.pair_potential(&[v_or_p * m.sigma()]).unwrap();
                        vec![
                            d,
                            dt,
                            dtt,
                            m.alpha_vdw(),
                            tc,
                            vc,
                            if fmt_ok { 1.0 } else { 0.0 },
                            if na_en { 1.0 } else { 0.0 },
                            pot[0],
                            force[0],
                        ]
                    }
                    other => panic!("unexpected BH kind: {other}"),
                }
            }
            "LJS-WCA" => {
                let mut m = LjsWca::new();
                if sigma > 0.0 && eps > 0.0 {
                    m.set_parameters(sigma, eps).unwrap();
                }
                match variant {
                    "default" => {}
                    "no-cavity" => m.set_control(LjsWcaControl {
                        enable_cavity: false,
                        ..LjsWcaControl::default()
                    }),
                    "no-hs" => m.set_control(LjsWcaControl {
                        enable_hard_sphere: false,
                        ..LjsWcaControl::default()
                    }),
                    "kolafa" => m.set_control(LjsWcaControl {
                        kolafa_pure_reference: true,
                        ..LjsWcaControl::default()
                    }),
                    "wca-a1-only" => m.set_control(LjsWcaControl {
                        enable_cavity: false,
                        enable_hard_sphere: false,
                        enable_a1: true,
                        enable_a2: false,
                        enable_a3: false,
                        enable_a4: false,
                        kolafa_pure_reference: false,
                    }),
                    "wca-a2-only" => m.set_control(LjsWcaControl {
                        enable_cavity: false,
                        enable_hard_sphere: false,
                        enable_a1: false,
                        enable_a2: true,
                        enable_a3: false,
                        enable_a4: false,
                        kolafa_pure_reference: false,
                    }),
                    "wca-a3-only" => m.set_control(LjsWcaControl {
                        enable_cavity: false,
                        enable_hard_sphere: false,
                        enable_a1: false,
                        enable_a2: false,
                        enable_a3: true,
                        enable_a4: false,
                        kolafa_pure_reference: false,
                    }),
                    "wca-a4-only" => m.set_control(LjsWcaControl {
                        enable_cavity: false,
                        enable_hard_sphere: false,
                        enable_a1: false,
                        enable_a2: false,
                        enable_a3: false,
                        enable_a4: true,
                        kolafa_pure_reference: false,
                    }),
                    other => panic!("unexpected WCA variant: {other}"),
                }
                match kind {
                    0 => tvn_slots(m.residual_helmholtz_tvn(t, v_or_p, n1).unwrap()).to_vec(),
                    1 => tvn_slots(m.dispersion_tvn(t, v_or_p, n1).unwrap()).to_vec(),
                    2 => tvn_slots(m.hard_sphere_tvn(t, v_or_p, n1).unwrap()).to_vec(),
                    3 => tvn_slots(m.soft_repulsion_tvn(t, v_or_p, n1).unwrap()).to_vec(),
                    4 => tp_slots(&m.state_tp(t, v_or_p, phase).unwrap()).to_vec(),
                    5 => {
                        let (d, dt, dtt) = m.diameter(t).unwrap();
                        let (tc, vc) = m.estimate_critical_parameters().unwrap();
                        let (fmt_ok, na_en) = m.fmt_compatibility();
                        let (pot, force) = m.pair_potential(&[v_or_p * m.sigma()]).unwrap();
                        vec![
                            d,
                            dt,
                            dtt,
                            m.alpha_vdw(),
                            tc,
                            vc,
                            if fmt_ok { 1.0 } else { 0.0 },
                            if na_en { 1.0 } else { 0.0 },
                            pot[0],
                            force[0],
                        ]
                    }
                    other => panic!("unexpected WCA kind: {other}"),
                }
            }
            "LJS-UV" => {
                let mut m = LjsUv::new();
                if sigma > 0.0 && eps > 0.0 {
                    m.set_parameters(sigma, eps).unwrap();
                }
                match variant {
                    "default" => {}
                    "temp-u" => m.set_control(LjsUvControl {
                        use_temperature_dependent_u_fraction: true,
                        ..LjsUvControl::default()
                    }),
                    "no-virial" => m.set_control(LjsUvControl {
                        enable_virial_term: false,
                        ..LjsUvControl::default()
                    }),
                    "no-cavity" => m.set_control(LjsUvControl {
                        enable_cavity: false,
                        ..LjsUvControl::default()
                    }),
                    "kolafa" => m.set_control(LjsUvControl {
                        kolafa_pure_reference: true,
                        ..LjsUvControl::default()
                    }),
                    other => panic!("unexpected UV variant: {other}"),
                }
                match kind {
                    0 => tvn_slots(m.residual_helmholtz_tvn(t, v_or_p, n1).unwrap()).to_vec(),
                    1 => tvn_slots(m.dispersion_tvn(t, v_or_p, n1).unwrap()).to_vec(),
                    2 => tvn_slots(m.hard_sphere_tvn(t, v_or_p, n1).unwrap()).to_vec(),
                    3 => tvn_slots(m.soft_repulsion_tvn(t, v_or_p, n1).unwrap()).to_vec(),
                    4 => tp_slots(&m.state_tp(t, v_or_p, phase).unwrap()).to_vec(),
                    5 => {
                        let (d, dt, dtt) = m.diameter(t).unwrap();
                        let (tc, vc) = m.estimate_critical_parameters().unwrap();
                        let (fmt_ok, na_en) = m.fmt_compatibility();
                        let (pot, force) = m.pair_potential(&[v_or_p * m.sigma()]).unwrap();
                        vec![
                            d,
                            dt,
                            dtt,
                            m.alpha_vdw(),
                            tc,
                            vc,
                            if fmt_ok { 1.0 } else { 0.0 },
                            if na_en { 1.0 } else { 0.0 },
                            pot[0],
                            force[0],
                        ]
                    }
                    other => panic!("unexpected UV kind: {other}"),
                }
            }
            "LJS-UF" | "LJ-UF" => {
                let pot_kind = if model == "LJ-UF" {
                    LjUfPotential::LennardJones
                } else {
                    LjUfPotential::LennardJonesSplined
                };
                let mut m = LjsUf::new(pot_kind);
                if sigma > 0.0 && eps > 0.0 {
                    m.set_parameters(sigma, eps).unwrap();
                }
                match variant {
                    "default" => {}
                    "no-cavity" => m.set_control(LjsUfControl {
                        enable_cavity: false,
                        ..LjsUfControl::default()
                    }),
                    "kolafa" => m.set_control(LjsUfControl {
                        kolafa_pure_reference: true,
                        ..LjsUfControl::default()
                    }),
                    other => panic!("unexpected UF variant: {other}"),
                }
                match kind {
                    0 => tvn_slots(m.residual_helmholtz_tvn(t, v_or_p, n1).unwrap()).to_vec(),
                    1 => tvn_slots(m.dispersion_tvn(t, v_or_p, n1).unwrap()).to_vec(),
                    2 => tvn_slots(m.hard_sphere_tvn(t, v_or_p, n1).unwrap()).to_vec(),
                    3 => tvn_slots(m.soft_repulsion_tvn(t, v_or_p, n1).unwrap()).to_vec(),
                    4 => tp_slots(&m.state_tp(t, v_or_p, phase).unwrap()).to_vec(),
                    other => panic!("unexpected UF kind: {other}"),
                }
            }
            other => panic!("unexpected model: {other}"),
        };

        assert_eq!(native.len(), source.len(), "{label}");
        let (atol, rtol) = match tag {
            "TP" => (1.0e-9, 1.0e-8),
            _ => (1.0e-12, 1.0e-10),
        };
        for (slot, (&s, &n)) in source.iter().zip(native.iter()).enumerate() {
            assert_close(&format!("{label}[{tag}:{slot}]"), s, n, atol, rtol);
            if s.abs() > 1.0e-12 {
                let rel = (s - n).abs() / s.abs();
                if tag == "TP" {
                    max_rel_tp = max_rel_tp.max(rel);
                } else {
                    max_rel_tvn = max_rel_tvn.max(rel);
                }
            }
        }
        checked += 1;
    }

    assert_eq!(checked, EXPECTED_FIXTURE_CASES);
    assert!(max_rel_tvn < 1.0e-10, "max_rel_tvn={max_rel_tvn:e}");
    assert!(max_rel_tp < 1.0e-8, "max_rel_tp={max_rel_tp:e}");
}

fn check_tvn_fd(
    label: &str,
    eval: impl Fn(f64, f64, f64) -> PureHelmholtz,
    t: f64,
    v: f64,
    n: f64,
) {
    let f = eval(t, v, n);
    let ht = 1.0e-6 * t;
    let hv = 1.0e-6 * v;
    let hn = 1.0e-6 * n;
    let fp_t = eval(t + ht, v, n);
    let fm_t = eval(t - ht, v, n);
    let fp_v = eval(t, v + hv, n);
    let fm_v = eval(t, v - hv, n);
    let fp_n = eval(t, v, n + hn);
    let fm_n = eval(t, v, n - hn);

    let checks = [
        ("dT", f.dt, (fp_t.value - fm_t.value) / (2.0 * ht)),
        ("dV", f.dv, (fp_v.value - fm_v.value) / (2.0 * hv)),
        ("dn", f.dn, (fp_n.value - fm_n.value) / (2.0 * hn)),
        ("dTT", f.dtt, (fp_t.dt - fm_t.dt) / (2.0 * ht)),
        ("dTV(V)", f.dtv, (fp_v.dt - fm_v.dt) / (2.0 * hv)),
        ("dTV(T)", f.dtv, (fp_t.dv - fm_t.dv) / (2.0 * ht)),
        ("dTn(n)", f.dtn, (fp_n.dt - fm_n.dt) / (2.0 * hn)),
        ("dTn(T)", f.dtn, (fp_t.dn - fm_t.dn) / (2.0 * ht)),
        ("dVV", f.dvv, (fp_v.dv - fm_v.dv) / (2.0 * hv)),
        ("dVn(n)", f.dvn, (fp_n.dv - fm_n.dv) / (2.0 * hn)),
        ("dVn(V)", f.dvn, (fp_v.dn - fm_v.dn) / (2.0 * hv)),
        ("dnn", f.dnn, (fp_n.dn - fm_n.dn) / (2.0 * hn)),
    ];
    for (name, ana, fd) in checks {
        assert_close(&format!("{label}.{name}"), fd, ana, 1.0e-6, 1.0e-5);
    }

    // Euler homogeneity of degree 1 for extensive F(T, V, n)
    assert_close(
        &format!("{label}.euler_F"),
        f.value,
        n * f.dn + v * f.dv,
        1.0e-11,
        1.0e-11,
    );
    assert_close(
        &format!("{label}.euler_FT"),
        f.dt,
        n * f.dtn + v * f.dtv,
        1.0e-11,
        1.0e-11,
    );
    assert_close(
        &format!("{label}.euler_0"),
        0.0,
        n * f.dnn + v * f.dvn,
        1.0e-10 * f.dnn.abs().max(1.0),
        1.0e-10,
    );
}

#[test]
fn finite_difference_tvn_and_euler_all_models() {
    let bh = LjsBh::new();
    check_tvn_fd(
        "BH-default",
        |t, v, n| bh.residual_helmholtz_tvn(t, v, n).unwrap(),
        120.0,
        1.0e-4,
        1.2,
    );

    let mut bh_laf = LjsBh::new();
    bh_laf.set_control(LjsControl {
        use_lafitte_a3: true,
        enable_chi_correction: false,
        ..LjsControl::default()
    });
    check_tvn_fd(
        "BH-lafitte-nochi",
        |t, v, n| bh_laf.residual_helmholtz_tvn(t, v, n).unwrap(),
        115.0,
        9.5e-5,
        1.3,
    );

    let wca = LjsWca::new();
    check_tvn_fd(
        "WCA-default",
        |t, v, n| wca.residual_helmholtz_tvn(t, v, n).unwrap(),
        110.0,
        9.0e-5,
        1.2,
    );
    // Above eta_c = 0.52 Taylor extrapolation branch
    check_tvn_fd(
        "WCA-high-eta",
        |t, v, n| wca.residual_helmholtz_tvn(t, v, n).unwrap(),
        80.0,
        4.3e-5,
        1.2,
    );

    let mut wca_csk = LjsWca::new();
    wca_csk.set_control(LjsWcaControl {
        kolafa_pure_reference: true,
        ..LjsWcaControl::default()
    });
    check_tvn_fd(
        "WCA-kolafa",
        |t, v, n| wca_csk.residual_helmholtz_tvn(t, v, n).unwrap(),
        115.0,
        9.5e-5,
        1.2,
    );

    let uv = LjsUv::new();
    check_tvn_fd(
        "UV-default",
        |t, v, n| uv.residual_helmholtz_tvn(t, v, n).unwrap(),
        125.0,
        9.5e-5,
        1.2,
    );

    let mut uv_t = LjsUv::new();
    uv_t.set_control(LjsUvControl {
        use_temperature_dependent_u_fraction: true,
        ..LjsUvControl::default()
    });
    check_tvn_fd(
        "UV-temp-u",
        |t, v, n| uv_t.residual_helmholtz_tvn(t, v, n).unwrap(),
        125.0,
        9.5e-5,
        1.2,
    );

    let ljs_uf = LjsUf::new(LjUfPotential::LennardJonesSplined);
    check_tvn_fd(
        "LJS-UF",
        |t, v, n| ljs_uf.residual_helmholtz_tvn(t, v, n).unwrap(),
        110.0,
        9.0e-5,
        1.2,
    );

    let lj_uf = LjsUf::new(LjUfPotential::LennardJones);
    check_tvn_fd(
        "LJ-UF",
        |t, v, n| lj_uf.residual_helmholtz_tvn(t, v, n).unwrap(),
        130.0,
        9.5e-5,
        1.2,
    );

    check_tvn_fd(
        "cavity-below-etac",
        |t, v, n| cavity_integral_lj_fres_tvn(3.42e-10, 124.0, t, v, n).unwrap(),
        120.0,
        9.0e-5,
        1.2,
    );
    check_tvn_fd(
        "cavity-above-etac",
        |t, v, n| cavity_integral_lj_fres_tvn(3.42e-10, 124.0, t, v, n).unwrap(),
        80.0,
        4.3e-5,
        1.2,
    );
}

#[test]
fn dft_additivity_decomposition_holds() {
    let (t, v, n) = (115.0, 9.5e-5, 1.3);

    // BH: molar(F) == hard_sphere + dispersion
    let bh = LjsBh::new();
    let bh_res = bh.residual_density(t, n / v).unwrap();
    let bh_molar = bh_res.value;
    let bh_hs = bh.hard_sphere_tvn(t, v, n).unwrap();
    let bh_disp = bh.dispersion_tvn(t, v, n).unwrap();
    assert_close(
        "BH additivity",
        bh_molar,
        bh_hs.value + bh_disp.value,
        1e-13,
        1e-13,
    );

    // WCA, UV, LJS-UF, LJ-UF: molar(F) == hard_sphere + soft_repulsion + dispersion
    let wca = LjsWca::new();
    let wca_tot = wca.residual_helmholtz_tvn(t, v, n).unwrap();
    let wca_hs = wca.hard_sphere_tvn(t, v, n).unwrap();
    let wca_sr = wca.soft_repulsion_tvn(t, v, n).unwrap();
    let wca_dp = wca.dispersion_tvn(t, v, n).unwrap();
    assert_close(
        "WCA F/n",
        wca_tot.value / n,
        wca_hs.value + wca_sr.value + wca_dp.value,
        1e-13,
        1e-13,
    );
    assert_close(
        "WCA F_T/n",
        wca_tot.dt / n,
        wca_hs.dt + wca_sr.dt + wca_dp.dt,
        1e-13,
        1e-13,
    );
    assert_close(
        "WCA F_V/n",
        wca_tot.dv / n,
        wca_hs.dv + wca_sr.dv + wca_dp.dv,
        1e-12,
        1e-12,
    );

    let uv = LjsUv::new();
    let uv_tot = uv.residual_helmholtz_tvn(t, v, n).unwrap();
    let uv_hs = uv.hard_sphere_tvn(t, v, n).unwrap();
    let uv_sr = uv.soft_repulsion_tvn(t, v, n).unwrap();
    let uv_dp = uv.dispersion_tvn(t, v, n).unwrap();
    assert_close(
        "UV F/n",
        uv_tot.value / n,
        uv_hs.value + uv_sr.value + uv_dp.value,
        1e-13,
        1e-13,
    );

    for pot in [
        LjUfPotential::LennardJonesSplined,
        LjUfPotential::LennardJones,
    ] {
        let uf = LjsUf::new(pot);
        let uf_tot = uf.residual_helmholtz_tvn(t, v, n).unwrap();
        let uf_hs = uf.hard_sphere_tvn(t, v, n).unwrap();
        let uf_sr = uf.soft_repulsion_tvn(t, v, n).unwrap();
        let uf_dp = uf.dispersion_tvn(t, v, n).unwrap();
        assert_close(
            "UF F/n",
            uf_tot.value / n,
            uf_hs.value + uf_sr.value + uf_dp.value,
            1e-13,
            1e-13,
        );
    }
}

fn check_tp_identities_and_fd(label: &str, eval: impl Fn(f64, f64) -> LjsTpState, t: f64, p: f64) {
    let s = eval(t, p);
    let ht = 1.0e-5 * t;
    let hp = 1.0e-5 * p;
    let st_p = eval(t + ht, p);
    let st_m = eval(t - ht, p);
    let sp_p = eval(t, p + hp);
    let sp_m = eval(t, p - hp);

    assert_close(
        &format!("{label}.dv_dt"),
        (st_p.molar_volume - st_m.molar_volume) / (2.0 * ht),
        s.dv_dt,
        1e-8,
        1e-5,
    );
    assert_close(
        &format!("{label}.dv_dp"),
        (sp_p.molar_volume - sp_m.molar_volume) / (2.0 * hp),
        s.dv_dp,
        1e-12,
        1e-5,
    );
    assert_close(
        &format!("{label}.dz_dt"),
        (st_p.compressibility.value - st_m.compressibility.value) / (2.0 * ht),
        s.compressibility.dt,
        1e-8,
        1e-5,
    );
    assert_close(
        &format!("{label}.dz_dp"),
        (sp_p.compressibility.value - sp_m.compressibility.value) / (2.0 * hp),
        s.compressibility.dp,
        1e-12,
        1e-5,
    );
    assert_close(
        &format!("{label}.dh_dt"),
        (st_p.enthalpy.value - st_m.enthalpy.value) / (2.0 * ht),
        s.enthalpy.dt,
        1e-6,
        1e-5,
    );
    assert_close(
        &format!("{label}.dh_dp"),
        (sp_p.enthalpy.value - sp_m.enthalpy.value) / (2.0 * hp),
        s.enthalpy.dp,
        1e-9,
        1e-5,
    );
    assert_close(
        &format!("{label}.ds_dt"),
        (st_p.entropy.value - st_m.entropy.value) / (2.0 * ht),
        s.entropy.dt,
        1e-8,
        1e-5,
    );
    assert_close(
        &format!("{label}.ds_dp"),
        (sp_p.entropy.value - sp_m.entropy.value) / (2.0 * hp),
        s.entropy.dp,
        1e-11,
        1e-5,
    );
    assert_close(
        &format!("{label}.dlnphi_dt"),
        (st_p.ln_fugacity_coefficient.value - st_m.ln_fugacity_coefficient.value) / (2.0 * ht),
        s.ln_fugacity_coefficient.dt,
        1e-8,
        1e-5,
    );
    assert_close(
        &format!("{label}.dlnphi_dp"),
        (sp_p.ln_fugacity_coefficient.value - sp_m.ln_fugacity_coefficient.value) / (2.0 * hp),
        s.ln_fugacity_coefficient.dp,
        1e-12,
        1e-5,
    );

    // Thermodynamic identities: Cp - Cv = -T * P_T^2 / P_V, dh/dP = v - T dv/dT
    let cp_minus_cv = -t * s.tv.pressure.dt * s.tv.pressure.dt / s.tv.pressure.dv;
    assert_close(
        &format!("{label}.Cp-Cv"),
        cp_minus_cv,
        s.cp_total - s.cv_total,
        1e-9,
        1e-10,
    );
    assert_close(
        &format!("{label}.Maxwell"),
        -s.dv_dt,
        s.entropy.dp,
        1e-12,
        1e-12,
    );
    assert!(s.speed_of_sound.is_finite() && s.speed_of_sound > 0.0);
}

#[test]
fn tp_derivatives_and_volume_round_trips_all_models() {
    let bh = LjsBh::new();
    let wca = LjsWca::new();
    let uv = LjsUv::new();
    let ljs_uf = LjsUf::new(LjUfPotential::LennardJonesSplined);
    let lj_uf = LjsUf::new(LjUfPotential::LennardJones);

    check_tp_identities_and_fd(
        "BH-liq",
        |t, p| bh.state_tp(t, p, Phase::Liquid).unwrap(),
        85.0,
        2.0e6,
    );
    check_tp_identities_and_fd(
        "WCA-liq",
        |t, p| wca.state_tp(t, p, Phase::Liquid).unwrap(),
        85.0,
        2.0e6,
    );
    check_tp_identities_and_fd(
        "UV-vap",
        |t, p| uv.state_tp(t, p, Phase::Vapor).unwrap(),
        135.0,
        1.2e6,
    );
    check_tp_identities_and_fd(
        "LJS-UF-liq",
        |t, p| ljs_uf.state_tp(t, p, Phase::Liquid).unwrap(),
        85.0,
        2.0e6,
    );
    check_tp_identities_and_fd(
        "LJ-UF-vap",
        |t, p| lj_uf.state_tp(t, p, Phase::Vapor).unwrap(),
        160.0,
        1.5e6,
    );

    // Verify LjsBh::volume matches ALL 8 source liquid/vapor volumes in tests/fixtures/ljs_state.csv
    // (previously broken for liquid roots due to eta > 1 pole in legacy solver).
    let mut checked = 0usize;
    for row in include_str!("fixtures/ljs_state.csv").lines().skip(1) {
        let fields: Vec<&str> = row.split(';').collect();
        let parts: Vec<&str> = fields[0].split('|').collect();
        let t: f64 = parts[0].parse().unwrap();
        let p: f64 = parts[1].parse().unwrap();
        let phase_idx: i32 = parts[2].parse().unwrap();
        let source_v: f64 = fields[1].parse().unwrap();
        let phase = if phase_idx == 1 {
            Phase::Liquid
        } else {
            Phase::Vapor
        };
        let native_v = bh.volume(t, p, phase).unwrap();
        assert_close(
            &format!("ljs_state[{t}|{p}|{phase_idx}]"),
            source_v,
            native_v,
            1e-12,
            1e-7,
        );
        let base_it = bh.baseline_volume_iteration(t, p, phase).unwrap();
        assert_close(
            &format!("ljs_state_base[{t}|{p}|{phase_idx}]"),
            source_v,
            base_it.volume,
            1e-12,
            1e-7,
        );
        assert!(base_it.volume_converged && base_it.pressure_converged);
        checked += 1;
    }
    assert_eq!(checked, 8);

    // Round trips across low-P vapor, compressed liquid, supercritical, and MinimumGibbs
    for (t, p, phase) in [
        (85.0, 2.0e6, Phase::Liquid),
        (90.0, 1.0e5, Phase::Vapor),
        (130.0, 1.0e6, Phase::Vapor),
        (250.0, 2.0e7, Phase::MinimumGibbs),
        (95.0, 6.0e5, Phase::MinimumGibbs),
        (200.0, 1.0e3, Phase::Vapor),
    ] {
        for (name, v, p_back, dp_dv) in [
            ("BH", bh.volume(t, p, phase).unwrap(), {
                let v = bh.volume(t, p, phase).unwrap();
                let pr = bh.pressure_tv(t, v, 1.0).unwrap();
                (pr.value, pr.dv)
            }),
            ("WCA", wca.volume(t, p, phase).unwrap(), {
                let v = wca.volume(t, p, phase).unwrap();
                let pr = wca.pressure_tv(t, v, 1.0).unwrap();
                (pr.value, pr.dv)
            }),
            ("UV", uv.volume(t, p, phase).unwrap(), {
                let v = uv.volume(t, p, phase).unwrap();
                let pr = uv.pressure_tv(t, v, 1.0).unwrap();
                (pr.value, pr.dv)
            }),
            ("LJS-UF", ljs_uf.volume(t, p, phase).unwrap(), {
                let v = ljs_uf.volume(t, p, phase).unwrap();
                let pr = ljs_uf.pressure_tv(t, v, 1.0).unwrap();
                (pr.value, pr.dv)
            }),
            ("LJ-UF", lj_uf.volume(t, p, phase).unwrap(), {
                let v = lj_uf.volume(t, p, phase).unwrap();
                let pr = lj_uf.pressure_tv(t, v, 1.0).unwrap();
                (pr.value, pr.dv)
            }),
        ]
        .into_iter()
        .map(|(n, v, (pb, dv))| (n, v, pb, dv))
        {
            assert_close(
                &format!("{name} round-trip T={t} P={p}"),
                p,
                p_back,
                1e-6,
                1e-9,
            );
            assert!(
                dp_dv < 0.0,
                "{name} unstable root at T={t}, P={p}, v={v}: dP/dV={dp_dv}"
            );
        }
    }
}

#[test]
fn edge_domains_standalone_kernels_and_invalid_inputs() {
    // Cavity eta_c = 0.52 continuity (C^2 in eta at eta_c)
    let (d, dt, dtt) = hard_sphere_diameter(3.42e-10, 124.0, 100.0);
    let below = cavity_integral_eta_t(3.42e-10, 124.0, d, dt, dtt, 0.52 - 1e-8, 100.0).unwrap();
    let at = cavity_integral_eta_t(3.42e-10, 124.0, d, dt, dtt, 0.52, 100.0).unwrap();
    let above = cavity_integral_eta_t(3.42e-10, 124.0, d, dt, dtt, 0.52 + 1e-8, 100.0).unwrap();
    assert_close("cavity C0 below", at.value, below.value, 1e-6, 1e-6);
    assert_close("cavity C0 above", at.value, above.value, 1e-6, 1e-6);
    assert_close("cavity C1 above", at.eta, above.eta, 1e-5, 1e-5);
    assert_close("cavity C2 above", at.eta_eta, above.eta_eta, 1e-4, 1e-4);

    // Pure HS CS and CSK invalid domain rejection (eta >= 1 or eta < 0)
    assert!(pure_hard_sphere_eta(-0.01, false).is_err());
    assert!(pure_hard_sphere_eta(1.0, false).is_err());
    assert!(pure_hard_sphere_eta(1.05, true).is_err());
    let (cs0, _, _) = pure_hard_sphere_eta(0.0, false).unwrap();
    let (csk0, _, _) = pure_hard_sphere_eta(0.0, true).unwrap();
    assert_close("CS(0)", 0.0, cs0, 1e-15, 1e-15);
    // Exact source CSK formula has constant offset 12.5/3 at eta=0.
    assert_close("CSK(0)", 12.5 / 3.0, csk0, 1e-14, 1e-14);

    // Standalone BH polynomial and WCA/UV TVN helpers agree with model methods
    let bh = LjsBh::new();
    let (t, v, n) = (130.0, 1.2e-4, 1.1);
    let (d_bh, _, _) = bh.diameter(t).unwrap();
    let x0 = bh.sigma() / d_bh;
    let rho_star = thermopack_rs::N_A * bh.sigma().powi(3) * n / v;
    assert_close(
        "a1_ljs_ex",
        bh.a1_tvn(t, v, n).unwrap().value,
        a1_ljs_ex(x0, rho_star, bh.eps()).unwrap(),
        1e-12,
        1e-12,
    );
    assert_close(
        "a2_ljs_ex",
        bh.a2_tvn(t, v, n).unwrap().value,
        a2_ljs_ex(x0, rho_star, bh.eps()).unwrap(),
        1e-12,
        1e-12,
    );
    assert_close(
        "a3_ljs_ex",
        bh.a3_tvn(t, v, n).unwrap().value,
        a3_ljs_ex(x0, rho_star, bh.eps()).unwrap(),
        1e-12,
        1e-12,
    );
    assert!(ai_ljs_ex(0, 0.05, 0.5, 1.0).is_err());
    assert!(ai_ljs_ex(5, 0.05, 0.5, 1.0).is_err());
    assert_eq!(bh.alpha_vdw(), ALPHA_LJS);

    let (dw, dwt, dwtt) = hard_sphere_diameter(3.42e-10, 124.0, t);
    let w1 = wca_term_tvn(1, 3.42e-10, 124.0, dw, dwt, dwtt, v, n).unwrap();
    let mut wca_a1 = LjsWca::new();
    wca_a1.set_control(LjsWcaControl {
        enable_cavity: false,
        enable_hard_sphere: false,
        enable_a1: true,
        enable_a2: false,
        enable_a3: false,
        enable_a4: false,
        kolafa_pure_reference: false,
    });
    assert_close(
        "wca_term_tvn/T",
        wca_a1.residual_helmholtz_tvn(t, v, n).unwrap().value,
        w1.value / t,
        1e-13,
        1e-13,
    );

    let uv_tvn = calc_uv_wca_tvn(3.42e-10, 124.0, t, v, n, dw, dwt, dwtt, true, false).unwrap();
    let mut uv_a1 = LjsUv::new();
    uv_a1.set_control(LjsUvControl {
        enable_hard_sphere: false,
        ..LjsUvControl::default()
    });
    assert_close(
        "calc_uv_wca_tvn",
        uv_a1.residual_helmholtz_tvn(t, v, n).unwrap().value,
        uv_tvn.value,
        1e-13,
        1e-13,
    );
    let u_mult = uv_a1_u_mult_t(124.0, dw / 3.42e-10, rho_star);
    assert!(u_mult.0.is_finite());

    // Extreme-density tanh saturation stays finite (no exp(2*x) NaN)
    let wca = LjsWca::new();
    let mut wca_disp_only = wca;
    wca_disp_only.set_control(LjsWcaControl {
        enable_hard_sphere: false,
        ..LjsWcaControl::default()
    });
    let extreme = wca_disp_only.residual_density(120.0, 2.0e5).unwrap();
    assert!(extreme.value.is_finite() && extreme.drho.is_finite());

    // Invalid inputs rejected cleanly
    assert!(bh.volume(-10.0, 1e5, Phase::Vapor).is_err());
    assert!(bh.volume(100.0, 0.0, Phase::Vapor).is_err());
    assert!(wca.residual_helmholtz_tvn(100.0, 0.0, 1.0).is_err());
    assert!(
        LjsUv::new()
            .residual_helmholtz_tvn(100.0, 1e-4, -1.0)
            .is_err()
    );
    assert!(
        LjsUf::new(LjUfPotential::LennardJones)
            .set_parameters(0.0, 124.0)
            .is_err()
    );

    // LjsTvState legacy conversion matches state()
    let tv: LjsTvState = wca.tv_state(120.0, 1.0e-4, 1.0).unwrap();
    let legacy = wca.state(120.0, 1.0e-4, 1.0).unwrap();
    assert_close("legacy h_total", tv.h_total, legacy.h_total, 1e-14, 1e-14);
}
