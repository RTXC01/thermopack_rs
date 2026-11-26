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
    calc_uv_wca_tvn, uv_a1_u_mult_t, uv_phi, uv_phi_temp,
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

#[test]
fn near_critical_both_spinodals_and_phase_selections() {
    // Sub-cell near-critical spinodal pairs where both grid endpoint derivatives
    // stay positive must still isolate both stable liquid and vapor roots.
    fn check_near_critical(
        label: &str,
        p: f64,
        roots: &[f64],
        v_liq: f64,
        v_vap: f64,
        v_min_g: f64,
        pressure_tv: impl Fn(f64) -> thermopack_rs::PurePressure,
    ) {
        assert_eq!(
            roots.len(),
            2,
            "{label}: expected 2 stable near-critical roots, got {roots:?}"
        );
        assert!(
            roots[0] < roots[1],
            "{label}: liquid root must be strictly smaller than vapor root ({roots:?})"
        );
        assert_close(&format!("{label}.liquid"), roots[0], v_liq, 1e-15, 1e-13);
        assert_close(&format!("{label}.vapor"), roots[1], v_vap, 1e-15, 1e-13);
        assert!(
            (v_min_g - roots[0]).abs() <= 1e-15 * roots[0]
                || (v_min_g - roots[1]).abs() <= 1e-15 * roots[1],
            "{label}: MinimumGibbs root {v_min_g:e} not among {roots:?}"
        );
        for (tag, v) in [("liquid", v_liq), ("vapor", v_vap)] {
            let pr = pressure_tv(v);
            assert!(
                pr.dv < 0.0,
                "{label}.{tag}: expected mechanically stable root (dP/dV < 0), got {}",
                pr.dv
            );
            let rel_err = (pr.value - p).abs() / p;
            assert!(
                rel_err < 1.0e-11,
                "{label}.{tag}: pressure closure rel_err={rel_err:.3e} (p={p:.16e}, got={:.16e})",
                pr.value
            );
        }
    }

    // WCA near critical point
    let wca = LjsWca::new();
    let (t, p) = (117.139920846527, 4_639_769.784019153);
    let roots = wca.volumes(t, p).unwrap();
    check_near_critical(
        "WCA-near-critical",
        p,
        &roots,
        wca.volume(t, p, Phase::Liquid).unwrap(),
        wca.volume(t, p, Phase::Vapor).unwrap(),
        wca.volume(t, p, Phase::MinimumGibbs).unwrap(),
        |v| wca.pressure_tv(t, v, 1.0).unwrap(),
    );
    assert_close("WCA-liq-val", 6.83918042e-5, roots[0], 1e-11, 1e-7);
    assert_close("WCA-vap-val", 6.86946043e-5, roots[1], 1e-11, 1e-7);

    // UV near critical point
    let uv = LjsUv::new();
    let (t, p) = (108.884564638390, 3_165_762.679804831);
    let roots = uv.volumes(t, p).unwrap();
    check_near_critical(
        "UV-near-critical",
        p,
        &roots,
        uv.volume(t, p, Phase::Liquid).unwrap(),
        uv.volume(t, p, Phase::Vapor).unwrap(),
        uv.volume(t, p, Phase::MinimumGibbs).unwrap(),
        |v| uv.pressure_tv(t, v, 1.0).unwrap(),
    );

    // LJS-UF near critical point
    let ljs_uf = LjsUf::new(LjUfPotential::LennardJonesSplined);
    let (t, p) = (114.282861338813, 4_377_015.446249527);
    let roots = ljs_uf.volumes(t, p).unwrap();
    check_near_critical(
        "LJS-UF-near-critical",
        p,
        &roots,
        ljs_uf.volume(t, p, Phase::Liquid).unwrap(),
        ljs_uf.volume(t, p, Phase::Vapor).unwrap(),
        ljs_uf.volume(t, p, Phase::MinimumGibbs).unwrap(),
        |v| ljs_uf.pressure_tv(t, v, 1.0).unwrap(),
    );

    // LJ-UF near critical point
    let lj_uf = LjsUf::new(LjUfPotential::LennardJones);
    let (t, p) = (163.366546486757, 5_652_595.510867001);
    let roots = lj_uf.volumes(t, p).unwrap();
    check_near_critical(
        "LJ-UF-near-critical",
        p,
        &roots,
        lj_uf.volume(t, p, Phase::Liquid).unwrap(),
        lj_uf.volume(t, p, Phase::Vapor).unwrap(),
        lj_uf.volume(t, p, Phase::MinimumGibbs).unwrap(),
        |v| lj_uf.pressure_tv(t, v, 1.0).unwrap(),
    );
}

#[test]
fn dilute_vapor_bracketing_and_pressure_relative_closure() {
    let bh = LjsBh::new();
    let t = 200.0;
    for p in [1.0e-13, 1.0e-9, 1.0e-8, 1.500_000_000_000_000_2e-8] {
        let v = bh.volume(t, p, Phase::Vapor).unwrap();
        let pr = bh.pressure_tv(t, v, 1.0).unwrap();
        let rel_err = (pr.value - p).abs() / p;
        assert!(
            rel_err < 1.0e-12,
            "BH dilute P={p:e}: got {:.16e} (rel_err={rel_err:.3e})",
            pr.value
        );
        assert!(pr.dv < 0.0);
    }

    // All remaining models close ultra-dilute vapor pressures cleanly
    let p = 1.0e-12;
    let wca = LjsWca::new();
    let v = wca.volume(t, p, Phase::Vapor).unwrap();
    assert!((wca.pressure_tv(t, v, 1.0).unwrap().value - p).abs() / p < 1.0e-12);

    let uv = LjsUv::new();
    let v = uv.volume(t, p, Phase::Vapor).unwrap();
    assert!((uv.pressure_tv(t, v, 1.0).unwrap().value - p).abs() / p < 1.0e-12);

    for pot in [
        LjUfPotential::LennardJonesSplined,
        LjUfPotential::LennardJones,
    ] {
        let uf = LjsUf::new(pot);
        let v = uf.volume(t, p, Phase::Vapor).unwrap();
        assert!((uf.pressure_tv(t, v, 1.0).unwrap().value - p).abs() / p < 1.0e-12);
    }
}

#[test]
fn extreme_tanh_second_derivatives_and_density_finiteness() {
    let p_uv = uv_phi(124.0, 120.0, 3.0);
    for (slot, val) in [p_uv.0, p_uv.1, p_uv.2, p_uv.3, p_uv.4, p_uv.5]
        .iter()
        .enumerate()
    {
        assert!(val.is_finite(), "uv_phi slot {slot} not finite: {val}");
    }
    let p_uv_t = uv_phi_temp(124.0, 120.0, 3.0);
    for (slot, val) in [p_uv_t.0, p_uv_t.1, p_uv_t.2, p_uv_t.3, p_uv_t.4, p_uv_t.5]
        .iter()
        .enumerate()
    {
        assert!(val.is_finite(), "uv_phi_temp slot {slot} not finite: {val}");
    }
    for is_lj in [false, true] {
        let uf = calc_uf_wca(10.0, 1.0, is_lj);
        for (slot, val) in [uf.0, uf.1, uf.2, uf.3, uf.4, uf.5].iter().enumerate() {
            assert!(
                val.is_finite(),
                "calc_uf_wca(is_lj={is_lj}) slot {slot} not finite: {val}"
            );
        }
    }

    // Extreme-temperature high-packing LJ-UF density evaluation stays finite
    let lj_uf = LjsUf::new(LjUfPotential::LennardJones);
    let prep = lj_uf.prepare(1.0e6).unwrap();
    let rho = 0.85 / prep.eta_per_rho();
    let dh = prep.residual_density(rho).unwrap();
    for (name, val) in [
        ("value", dh.value),
        ("dt", dh.dt),
        ("drho", dh.drho),
        ("dtt", dh.dtt),
        ("dtrho", dh.dtrho),
        ("drhorho", dh.drhorho),
    ] {
        assert!(
            val.is_finite(),
            "LJ-UF extreme density {name} not finite: {val}"
        );
    }
}

#[test]
fn legacy_nonpositive_pressure_state_error_and_unstable_sound_speed() {
    let bh = LjsBh::new();

    // Negative-pressure state: legacy state() returns Err(InvalidInput("pressure")),
    // while tv_state() exposes the negative pressure and NaN entropies explicitly.
    let (t_neg, v_neg) = (85.0, 7.904_664_559_912_66e-5);
    let err = bh.state(t_neg, v_neg, 1.0).unwrap_err();
    assert!(
        matches!(err, thermopack_rs::Error::InvalidInput("pressure")),
        "expected InvalidInput(\"pressure\"), got {err:?}"
    );
    let tv_neg = bh.tv_state(t_neg, v_neg, 1.0).unwrap();
    assert!(tv_neg.pressure.value < -1.0e6);
    assert!(tv_neg.s_residual.is_nan() && tv_neg.s_total.is_nan());

    // Consistent positive-pressure check on WCA, UV, and UF legacy state()
    assert!(matches!(
        LjsWca::new().state(t_neg, v_neg, 1.0),
        Err(thermopack_rs::Error::InvalidInput("pressure"))
    ));
    assert!(matches!(
        LjsUv::new().state(t_neg, v_neg, 1.0),
        Err(thermopack_rs::Error::InvalidInput("pressure"))
    ));
    assert!(matches!(
        LjsUf::new(LjUfPotential::LennardJonesSplined).state(t_neg, v_neg, 1.0),
        Err(thermopack_rs::Error::InvalidInput("pressure"))
    ));

    // Mechanically unstable TV state where both pressure.dv > 0 and cp_total < 0:
    // speed_of_sound must be NaN (never sqrt of positive product of two wrong signs).
    let prep = bh.prepare(85.0).unwrap();
    let v_unstable = prep.eta_per_rho() / 0.10;
    let tv_unstable = prep.tv_state(v_unstable, 1.0).unwrap();
    assert!(tv_unstable.pressure.dv > 0.0, "expected dP/dV > 0");
    assert!(tv_unstable.cp_total < 0.0, "expected cp_total < 0");
    assert!(
        tv_unstable.speed_of_sound.is_nan(),
        "unstable TV state must report NaN speed_of_sound, got {}",
        tv_unstable.speed_of_sound
    );
    assert!(
        bh.speed_of_sound_tv(85.0, v_unstable, 1.0)
            .unwrap()
            .is_nan()
    );
}

#[test]
fn fresh_and_prepared_evaluations_agree_across_models() {
    let (t, v, n, p) = (125.0, 1.1e-4, 1.2, 1.2e6);

    let bh = LjsBh::new();
    let bh_p = bh.prepare(t).unwrap();
    assert_close(
        "BH fresh vs prep F",
        bh.residual_helmholtz_tvn(t, v, n).unwrap().value,
        bh_p.residual_helmholtz_tvn(v, n).unwrap().value,
        0.0,
        1e-15,
    );
    assert_close(
        "BH fresh vs prep V",
        bh.volume(t, p, Phase::Vapor).unwrap(),
        bh_p.volume(p, Phase::Vapor).unwrap(),
        0.0,
        1e-15,
    );

    let wca = LjsWca::new();
    let wca_p = wca.prepare(t).unwrap();
    assert_close(
        "WCA fresh vs prep F",
        wca.residual_helmholtz_tvn(t, v, n).unwrap().value,
        wca_p.residual_helmholtz_tvn(v, n).unwrap().value,
        0.0,
        1e-15,
    );
    assert_close(
        "WCA fresh vs prep V",
        wca.volume(t, p, Phase::Vapor).unwrap(),
        wca_p.volume(p, Phase::Vapor).unwrap(),
        0.0,
        1e-15,
    );

    let uv = LjsUv::new();
    let uv_p = uv.prepare(t).unwrap();
    assert_close(
        "UV fresh vs prep F",
        uv.residual_helmholtz_tvn(t, v, n).unwrap().value,
        uv_p.residual_helmholtz_tvn(v, n).unwrap().value,
        0.0,
        1e-15,
    );
    assert_close(
        "UV fresh vs prep V",
        uv.volume(t, p, Phase::Vapor).unwrap(),
        uv_p.volume(p, Phase::Vapor).unwrap(),
        0.0,
        1e-15,
    );

    let uf = LjsUf::new(LjUfPotential::LennardJonesSplined);
    let uf_p = uf.prepare(t).unwrap();
    assert_close(
        "UF fresh vs prep F",
        uf.residual_helmholtz_tvn(t, v, n).unwrap().value,
        uf_p.residual_helmholtz_tvn(v, n).unwrap().value,
        0.0,
        1e-15,
    );
    assert_close(
        "UF fresh vs prep V",
        uf.volume(t, p, Phase::Vapor).unwrap(),
        uf_p.volume(p, Phase::Vapor).unwrap(),
        0.0,
        1e-15,
    );
}

#[test]
fn dense_low_pressure_cancellation_flat_root_close_gibbs_and_overflow_rejection() {
    // 1. Dense low-pressure liquid branch preserved under ideal/residual pressure cancellation
    // (BH T=85 K, P=0.1 Pa): both dense liquid and dilute vapor roots must be returned,
    // with the liquid root straddled by public forward pressure_tv().
    let bh = LjsBh::new();
    let (t_bh, p_bh) = (85.0, 0.1);
    let bh_roots = bh.volumes(t_bh, p_bh).unwrap();
    assert_eq!(
        bh_roots.len(),
        2,
        "BH T={t_bh}, P={p_bh}: expected both liquid and vapor roots, got {bh_roots:?}"
    );
    let v_liq = bh.volume(t_bh, p_bh, Phase::Liquid).unwrap();
    let v_vap = bh.volume(t_bh, p_bh, Phase::Vapor).unwrap();
    assert_close("BH dense low-P liquid", bh_roots[0], v_liq, 0.0, 1e-15);
    assert_close("BH dilute low-P vapor", bh_roots[1], v_vap, 0.0, 1e-15);
    assert!(
        v_liq > 3.336_851_694_343_900e-5 && v_liq < 3.336_851_694_343_950e-5,
        "BH liquid root {v_liq:.16e} outside straddled crossing"
    );
    let pr_lo = bh.pressure_tv(t_bh, 3.336_851_694_343_900e-5, 1.0).unwrap();
    let pr_hi = bh.pressure_tv(t_bh, 3.336_851_694_343_950e-5, 1.0).unwrap();
    assert!(pr_lo.value > p_bh && pr_hi.value < p_bh);
    assert!(pr_lo.dv < 0.0 && pr_hi.dv < 0.0);
    let pr_liq = bh.pressure_tv(t_bh, v_liq, 1.0).unwrap();
    assert!(pr_liq.dv < 0.0);
    assert!(
        (pr_liq.value - p_bh).abs() < 1.0e-6,
        "BH liquid residual too large: {:.3e} Pa",
        (pr_liq.value - p_bh).abs()
    );
    // Dilute vapor branch remains strictly pressure-relative
    let pr_vap = bh.pressure_tv(t_bh, v_vap, 1.0).unwrap();
    assert!(pr_vap.dv < 0.0);
    assert!(
        (pr_vap.value - p_bh).abs() / p_bh < 1.0e-12,
        "BH dilute vapor closure rel_err={:.3e}",
        (pr_vap.value - p_bh).abs() / p_bh
    );

    // 2. Flat near-critical resolvable pressure root (LJ-UF T=163.36664638716522 K,
    // P=5652616.747907021 Pa): returned liquid volume must converge to the local
    // pressure-evaluation resolution (~1e-8 Pa noise floor on P ≈ 5.65e6 Pa) rather
    // than stopping early at ~3.27e-7 Pa / 7.266782e-5.
    let lj_uf = LjsUf::new(LjUfPotential::LennardJones);
    let (t_flat, p_flat) = (163.366_646_387_165_22, 5_652_616.747_907_021);
    let flat_roots = lj_uf.volumes(t_flat, p_flat).unwrap();
    assert_eq!(flat_roots.len(), 2, "LJ-UF flat roots: {flat_roots:?}");
    let v_flat_liq = lj_uf.volume(t_flat, p_flat, Phase::Liquid).unwrap();
    assert_close("LJ-UF flat liq", flat_roots[0], v_flat_liq, 0.0, 1e-15);
    // Independent forward-pressure bracket at [7.266895e-5, 7.266902e-5] straddles
    // the crossing; the evaluation noise floor (~1.5e-8 Pa with dP/dV ≈ -50 Pa/(m^3/mol))
    // spans [7.266885e-5, 7.266915e-5] before the spinodal near 7.266917e-5.
    let (p_fast_lo, _) = lj_uf.pressure(t_flat, 7.266_895e-5, 1.0).unwrap();
    let (p_fast_hi, _) = lj_uf.pressure(t_flat, 7.266_902e-5, 1.0).unwrap();
    let p_br_lo = lj_uf.pressure_tv(t_flat, 7.266_895e-5, 1.0).unwrap();
    let p_br_hi = lj_uf.pressure_tv(t_flat, 7.266_902e-5, 1.0).unwrap();
    assert!(
        p_fast_lo > p_flat
            && p_fast_hi < p_flat
            && p_br_lo.value > p_flat
            && p_br_hi.value < p_flat
            && p_br_lo.dv < 0.0
            && p_br_hi.dv < 0.0
    );
    let (p_fast, _) = lj_uf.pressure(t_flat, v_flat_liq, 1.0).unwrap();
    let pr_flat = lj_uf.pressure_tv(t_flat, v_flat_liq, 1.0).unwrap();
    assert!(
        v_flat_liq > 7.266_885e-5 && v_flat_liq < 7.266_915e-5,
        "LJ-UF flat liquid root {v_flat_liq:.16e} fell outside evaluation-resolution band"
    );
    assert!(pr_flat.dv < 0.0);
    assert!(
        (p_fast - p_flat).abs() < 2.0e-8 && (pr_flat.value - p_flat).abs() < 2.0e-8,
        "LJ-UF flat residual not at evaluation resolution: fast={:.3e}, tv={:.3e}",
        p_fast - p_flat,
        pr_flat.value - p_flat
    );

    // 3. Close-phase MinimumGibbs selection via Maxwell pressure integral
    // (LJ-UF T=163.36664548723948 K, P=5652616.556625032 Pa):
    // independent trapezoid integration of (P(T, V) - P_spec) dV gives G_liq - G_vap < 0,
    // so MinimumGibbs must choose the liquid root; slightly lower pressure chooses vapor.
    let (t_gibbs, p_gibbs) = (163.366_645_487_239_48, 5_652_616.556_625_032);
    let g_roots = lj_uf.volumes(t_gibbs, p_gibbs).unwrap();
    assert_eq!(g_roots.len(), 2, "LJ-UF close Gibbs roots: {g_roots:?}");
    let independent_dg = |v1: f64, v2: f64, p_ref: f64, n: usize| -> f64 {
        let h = (v2 - v1) / n as f64;
        let mut s = 0.0;
        for i in 0..=n {
            let v = v1 + h * i as f64;
            let w = if i == 0 || i == n { 0.5 } else { 1.0 };
            let (pv, _) = lj_uf.pressure(t_gibbs, v, 1.0).unwrap();
            s += w * (pv - p_ref) * h;
        }
        s
    };
    let dg_256 = independent_dg(g_roots[0], g_roots[1], p_gibbs, 256);
    let dg_1024 = independent_dg(g_roots[0], g_roots[1], p_gibbs, 1024);
    assert!(
        dg_256 < -5.0e-13 && dg_1024 < -5.0e-13,
        "expected G_liq - G_vap ≈ -7.77e-13 J/mol, got dg_256={dg_256:.3e}, dg_1024={dg_1024:.3e}"
    );
    let v_min_g = lj_uf.volume(t_gibbs, p_gibbs, Phase::MinimumGibbs).unwrap();
    assert_close(
        "LJ-UF close MinimumGibbs liquid",
        g_roots[0],
        v_min_g,
        0.0,
        1e-15,
    );

    // Shift pressure slightly below coexistence: G_liq - G_vap > 0 so MinimumGibbs picks vapor
    let p_vap_favored = p_gibbs - 4.0e-5;
    let gv_roots = lj_uf.volumes(t_gibbs, p_vap_favored).unwrap();
    assert_eq!(gv_roots.len(), 2);
    let dg_vap = independent_dg(gv_roots[0], gv_roots[1], p_vap_favored, 256);
    assert!(
        dg_vap > 0.0,
        "expected vapor-favored dg > 0, got {dg_vap:.3e}"
    );
    let v_min_g_vap = lj_uf
        .volume(t_gibbs, p_vap_favored, Phase::MinimumGibbs)
        .unwrap();
    assert_close(
        "LJ-UF close MinimumGibbs vapor",
        gv_roots[1],
        v_min_g_vap,
        0.0,
        1e-15,
    );

    // 4. Reject unrepresentable molar volume (overflow to +inf) instead of returning Ok([inf])
    let err_vols = bh.volumes(200.0, 1.0e-306).unwrap_err();
    assert!(
        matches!(err_vols, thermopack_rs::Error::InvalidInput(_)),
        "expected InvalidInput on volume overflow, got {err_vols:?}"
    );
    let err_vap = bh.volume(200.0, 1.0e-306, Phase::Vapor).unwrap_err();
    assert!(
        matches!(err_vap, thermopack_rs::Error::InvalidInput(_)),
        "expected InvalidInput on vapor volume overflow, got {err_vap:?}"
    );
    assert!(LjsWca::new().volumes(200.0, 1.0e-306).is_err());
    assert!(LjsUv::new().volumes(200.0, 1.0e-306).is_err());
    assert!(lj_uf.volumes(200.0, 1.0e-306).is_err());
}

#[test]
fn vapor_volume_overflow_does_not_select_liquid() {
    let model = LjsBh::new();
    let temperature = 85.0;
    let pressure = 1.0e-306;
    let liquid = model.volume(temperature, pressure, Phase::Liquid).unwrap();
    assert_close(
        "dense liquid volume",
        liquid,
        3.336_851_696_116_317e-5,
        1e-12,
        0.0,
    );
    for result in [
        model.volumes(temperature, pressure).map(|_| ()),
        model
            .volume(temperature, pressure, Phase::Vapor)
            .map(|_| ()),
        model
            .volume(temperature, pressure, Phase::MinimumGibbs)
            .map(|_| ()),
    ] {
        assert!(
            matches!(result, Err(thermopack_rs::Error::InvalidInput(_))),
            "{result:?}"
        );
    }

    let pressure = 3.970_625_475_093_25e-306;
    let roots = model.volumes(temperature, pressure).unwrap();
    assert_eq!(roots.len(), 2);
    assert_close("dense liquid volume", roots[0], liquid, 1e-12, 0.0);
    let vapor = roots[1];
    assert!(vapor.is_finite() && vapor > 1.0e308);
    let round_trip = model.pressure_tv(temperature, vapor, 1.0).unwrap().value;
    assert!((round_trip / pressure - 1.0).abs() < 1e-8, "{round_trip:e}");
    for phase in [Phase::Vapor, Phase::MinimumGibbs] {
        assert_eq!(model.volume(temperature, pressure, phase).unwrap(), vapor);
    }
}
