//! Integration tests for associating and polar PC-SAFT (`PC-SAFT`, `sPC-SAFT`,
//! `PCP-SAFT`, `sPCP-SAFT`), shared Wertheim TPT1 association, Gross/Vrabec
//! multipolar terms, and simplified PC-SAFT RDF.

use thermopack_rs::association::{self, AssocCombRule, AssocScheme, DeltaCombRule};
use thermopack_rs::multipol::MultipolParams;
use thermopack_rs::pc_saft::{PcSaft, PcSaftKind, PcSaftOptions, RGAS};
use thermopack_rs::rdf::spcsaft_rdf;
use thermopack_rs::{Error, Phase};

fn approx_eq(actual: f64, expected: f64, rtol: f64, atol: f64, label: &str) {
    assert!(
        actual.is_finite() && expected.is_finite(),
        "{label}: non-finite actual={actual} expected={expected}"
    );
    let err = (actual - expected).abs();
    let tol = atol + rtol * expected.abs().max(actual.abs());
    assert!(
        err <= tol,
        "{label}: actual={actual:.15e} expected={expected:.15e} err={err:.3e} tol={tol:.3e}"
    );
}

fn physical_scales(nc: usize, nsites: usize, sumn: f64, t: f64, p: f64, v_tv: f64) -> Vec<f64> {
    let v_tp = v_tv / sumn;
    let rt = RGAS * t;
    let mut s = Vec::new();
    // F, F_T, F_V, F_n
    s.extend([sumn, sumn / t, sumn / v_tv]);
    s.extend(std::iter::repeat_n(1.0, nc));
    // F_TT, F_VV, F_TV, F_Tn, F_Vn, F_nn
    s.extend([sumn / (t * t), sumn / (v_tv * v_tv), sumn / (t * v_tv)]);
    s.extend(std::iter::repeat_n(1.0 / t, nc));
    s.extend(std::iter::repeat_n(1.0 / v_tv, nc));
    s.extend(std::iter::repeat_n(1.0 / sumn, nc * nc));
    // X_k
    s.extend(std::iter::repeat_n(1.0, nsites));
    // Z
    s.extend([1.0, 1.0 / t, 1.0 / p]);
    s.extend(std::iter::repeat_n(1.0, nc));
    // v_tp
    s.extend([v_tp, v_tp / t, v_tp / p]);
    s.extend(std::iter::repeat_n(v_tp, nc));
    // h_res
    s.extend([rt, RGAS, v_tp]);
    s.extend(std::iter::repeat_n(rt, nc));
    // h
    s.extend([rt, RGAS, v_tp]);
    s.extend(std::iter::repeat_n(rt, nc));
    // s_res
    s.extend([RGAS, RGAS / t, v_tp / t]);
    s.extend(std::iter::repeat_n(RGAS, nc));
    // s
    s.extend([RGAS, RGAS / t, v_tp / t]);
    s.extend(std::iter::repeat_n(RGAS, nc));
    // ln_phi, dt, dp, dn
    s.extend(std::iter::repeat_n(1.0, nc));
    s.extend(std::iter::repeat_n(1.0 / t, nc));
    s.extend(std::iter::repeat_n(1.0 / p, nc));
    s.extend(std::iter::repeat_n(1.0, nc * nc));
    // speed of sound
    s.push(100.0);
    s
}

#[test]
fn offline_exact_source_fixture_parity_all_35_cases() {
    let mut lines = include_str!("fixtures/pcsaft_full.csv").lines();
    assert_eq!(lines.next().unwrap(), "case;nsites;values");
    let mut compared = 0usize;
    let rtol = 1.0e-8;
    let atol_scale = 1.0e-11;

    for row in lines {
        let cols: Vec<&str> = row.split(';').collect();
        let fields: Vec<&str> = cols[0].split('|').collect();
        assert_eq!(fields[0], "PCSAFT_FULL");
        let kind = match fields[1] {
            "PC-SAFT" => PcSaftKind::PcSaft,
            "sPC-SAFT" => PcSaftKind::SPcSaft,
            "PCP-SAFT" => PcSaftKind::PcpSaft,
            "sPCP-SAFT" => PcSaftKind::SPcpSaft,
            other => panic!("unexpected kind {other}"),
        };
        let param_ref = fields[2];
        let comps: Vec<&str> = fields[3].split(',').collect();
        let t: f64 = fields[4].parse().unwrap();
        let p: f64 = fields[5].parse().unwrap();
        let v_tv: f64 = fields[6].parse().unwrap();
        let moles: Vec<f64> = fields[7].split(',').map(|v| v.parse().unwrap()).collect();
        let phase = if fields[8] == "L" {
            Phase::Liquid
        } else {
            Phase::Vapor
        };
        let nsites: usize = cols[1].parse().unwrap();
        let expected: Vec<f64> = cols[2..].iter().map(|v| v.parse().unwrap()).collect();

        let nc = comps.len();
        let sumn: f64 = moles.iter().sum();
        let x: Vec<f64> = moles.iter().map(|&v| v / sumn).collect();

        let model = PcSaft::new_with_options(kind, &comps, param_ref).unwrap();
        let f = model.residual_helmholtz_tv(t, v_tv, &moles).unwrap();
        assert_eq!(f.x_k.len(), nsites);
        let s = model.state_tp(t, p, &x, phase).unwrap();

        let mut actual = Vec::with_capacity(expected.len());
        actual.extend([f.f, f.f_t, f.f_v]);
        actual.extend(f.f_n);
        actual.extend([f.f_tt, f.f_vv, f.f_tv]);
        actual.extend(f.f_tn);
        actual.extend(f.f_vn);
        actual.extend(f.f_nn);
        actual.extend(f.x_k);
        actual.extend([
            s.compressibility_factor.value,
            s.compressibility_factor.dt,
            s.compressibility_factor.dp,
        ]);
        actual.extend(s.compressibility_factor.dn);
        actual.extend([
            s.specific_volume.value,
            s.specific_volume.dt,
            s.specific_volume.dp,
        ]);
        actual.extend(s.specific_volume.dn);
        actual.extend([
            s.residual_enthalpy.value,
            s.residual_enthalpy.dt,
            s.residual_enthalpy.dp,
        ]);
        actual.extend(s.residual_enthalpy.dn);
        actual.extend([s.enthalpy.value, s.enthalpy.dt, s.enthalpy.dp]);
        actual.extend(s.enthalpy.dn);
        actual.extend([
            s.residual_entropy.value,
            s.residual_entropy.dt,
            s.residual_entropy.dp,
        ]);
        actual.extend(s.residual_entropy.dn);
        actual.extend([s.entropy.value, s.entropy.dt, s.entropy.dp]);
        actual.extend(s.entropy.dn);
        actual.extend(s.ln_fugacity.ln_phi);
        actual.extend(s.ln_fugacity.dt);
        actual.extend(s.ln_fugacity.dp);
        actual.extend(s.ln_fugacity.dn);
        actual.push(s.speed_of_sound);

        assert_eq!(actual.len(), expected.len());
        let scales = physical_scales(nc, nsites, sumn, t, p, v_tv);
        for (idx, ((&a, &e), &sc)) in actual
            .iter()
            .zip(expected.iter())
            .zip(scales.iter())
            .enumerate()
        {
            approx_eq(
                a,
                e,
                rtol,
                atol_scale * sc,
                &format!("{}[{}]", cols[0], idx),
            );
        }
        compared += 1;
    }
    assert_eq!(compared, 35);
}

fn check_tvn_10_blocks_fd(model: &PcSaft, t: f64, v: f64, n: &[f64], label: &str) {
    let nc = n.len();
    let base = model.residual_helmholtz_tv(t, v, n).unwrap();

    let ht = 1.0e-5 * t;
    let hv = 1.0e-5 * v;
    let fp_t = model.residual_helmholtz_tv(t + ht, v, n).unwrap();
    let fm_t = model.residual_helmholtz_tv(t - ht, v, n).unwrap();
    let fp_v = model.residual_helmholtz_tv(t, v + hv, n).unwrap();
    let fm_v = model.residual_helmholtz_tv(t, v - hv, n).unwrap();

    approx_eq(
        base.f_t,
        (fp_t.f - fm_t.f) / (2.0 * ht),
        1e-6,
        1e-9,
        &format!("{label} F_T"),
    );
    approx_eq(
        base.f_v,
        (fp_v.f - fm_v.f) / (2.0 * hv),
        1e-6,
        1e-9 / v,
        &format!("{label} F_V"),
    );
    approx_eq(
        base.f_tt,
        (fp_t.f_t - fm_t.f_t) / (2.0 * ht),
        1e-6,
        1e-9 / (t * t),
        &format!("{label} F_TT"),
    );
    approx_eq(
        base.f_vv,
        (fp_v.f_v - fm_v.f_v) / (2.0 * hv),
        1e-6,
        1e-9 / (v * v),
        &format!("{label} F_VV"),
    );
    approx_eq(
        base.f_tv,
        (fp_t.f_v - fm_t.f_v) / (2.0 * ht),
        1e-6,
        1e-9 / (t * v),
        &format!("{label} F_TV"),
    );

    let sumn: f64 = n.iter().sum();
    for i in 0..nc {
        let hn = 1.0e-5 * n[i].max(1.0e-2 * sumn);
        let mut np = n.to_vec();
        let mut nm = n.to_vec();
        np[i] += hn;
        nm[i] -= hn;
        let fp = model.residual_helmholtz_tv(t, v, &np).unwrap();
        let fm = model.residual_helmholtz_tv(t, v, &nm).unwrap();

        approx_eq(
            base.f_n[i],
            (fp.f - fm.f) / (2.0 * hn),
            1e-6,
            1e-9,
            &format!("{label} F_n[{i}]"),
        );
        approx_eq(
            base.f_tn[i],
            (fp.f_t - fm.f_t) / (2.0 * hn),
            1e-6,
            1e-9 / t,
            &format!("{label} F_Tn[{i}]"),
        );
        approx_eq(
            base.f_vn[i],
            (fp.f_v - fm.f_v) / (2.0 * hn),
            1e-6,
            1e-9 / v,
            &format!("{label} F_Vn[{i}]"),
        );
        for j in 0..nc {
            approx_eq(
                base.f_nn[i * nc + j],
                (fp.f_n[j] - fm.f_n[j]) / (2.0 * hn),
                1e-6,
                1e-9 / sumn,
                &format!("{label} F_nn[{i},{j}]"),
            );
            approx_eq(
                base.f_nn[i * nc + j],
                base.f_nn[j * nc + i],
                1e-12,
                1e-12 / sumn,
                &format!("{label} F_nn sym [{i},{j}]"),
            );
        }
    }

    // Euler homogeneity identity: V * F_V + sum_i n_i * F_n[i] == F
    let euler = v * base.f_v
        + n.iter()
            .zip(&base.f_n)
            .map(|(ni, fni)| ni * fni)
            .sum::<f64>();
    approx_eq(euler, base.f, 1e-11, 1e-11, &format!("{label} Euler"));
}

#[test]
fn finite_difference_all_10_blocks_across_assoc_polar_and_rules() {
    // 1. Original PC-SAFT associating mixture (4C + 2B)
    let m1 = PcSaft::new(&["H2O", "MEOH"]).unwrap();
    check_tvn_10_blocks_fd(&m1, 340.0, 3.2e-5, &[0.6, 0.4], "PC-SAFT H2O+MEOH");

    // 2. Simplified sPC-SAFT with ariComb + geoComb (NH3 + H2O)
    let m2 = PcSaft::new_simplified(&["NH3", "H2O"]).unwrap();
    check_tvn_10_blocks_fd(&m2, 310.0, 2.8e-5, &[0.4, 0.6], "sPC-SAFT NH3+H2O");

    // 3. Polar PCP-SAFT D-D + Q-Q + D-Q ternary
    let m3 = PcSaft::new_polar(&["CO2", "C1", "ACETONE"]).unwrap();
    check_tvn_10_blocks_fd(
        &m3,
        330.0,
        9.0e-5,
        &[0.3, 0.3, 0.4],
        "PCP-SAFT CO2+C1+ACETONE",
    );

    // 4. Polar + associating PCP-SAFT (original BMCSL pair RDF!)
    let m4 = PcSaft::new_polar(&["NH3", "H2O"]).unwrap();
    check_tvn_10_blocks_fd(&m4, 320.0, 2.8e-5, &[0.4, 0.6], "PCP-SAFT NH3+H2O");

    // 5. Elliot cross-association rule
    let mut m5 = PcSaft::new(&["H2O", "MEOH"]).unwrap();
    m5.set_delta_comb_rule(DeltaCombRule::Elliot).unwrap();
    check_tvn_10_blocks_fd(&m5, 340.0, 3.2e-5, &[0.6, 0.4], "PC-SAFT Elliot H2O+MEOH");

    // 6. Induced association (Scheme1ea on CO2 cross-associating with H2O 4C)
    let mut m6 = PcSaft::new(&["H2O", "CO2"]).unwrap();
    m6.set_assoc_params(1, 0.0, 0.05, AssocScheme::Scheme1ea)
        .unwrap();
    m6.set_assoc_comb_rules(0, 1, AssocCombRule::AriComb, AssocCombRule::GeoComb)
        .unwrap();
    check_tvn_10_blocks_fd(&m6, 350.0, 3.0e-5, &[0.7, 0.3], "PC-SAFT Scheme1ea H2O+CO2");
}

#[test]
fn pcp_saft_associating_mixture_is_component_order_invariant() {
    // Upstream Fortran `saft_association.f90:136,147` omits `eosPCP_SAFT` and
    // breaks component-order invariance by reusing g(1, 1) for every pair.
    // Rust evaluates the physical per-pair g(i, j) and preserves symmetry.
    let m_ab = PcSaft::new_polar(&["CO2", "H2S"]).unwrap();
    let m_ba = PcSaft::new_polar(&["H2S", "CO2"]).unwrap();
    let h_ab = m_ab
        .residual_helmholtz_tv(290.0, 5.744_117_894_993_951e-5, &[0.5, 0.5])
        .unwrap();
    let h_ba = m_ba
        .residual_helmholtz_tv(290.0, 5.744_117_894_993_951e-5, &[0.5, 0.5])
        .unwrap();
    approx_eq(h_ab.f, h_ba.f, 1e-13, 1e-14, "PCP-SAFT CO2/H2S order F");
    approx_eq(
        h_ab.f_v,
        h_ba.f_v,
        1e-13,
        1e-10,
        "PCP-SAFT CO2/H2S order F_V",
    );
    approx_eq(h_ab.f_n[0], h_ba.f_n[1], 1e-13, 1e-13, "PCP-SAFT CO2 F_n");
    approx_eq(h_ab.f_n[1], h_ba.f_n[0], 1e-13, 1e-13, "PCP-SAFT H2S F_n");

    let n_ab = PcSaft::new_polar(&["NH3", "H2O"]).unwrap();
    let n_ba = PcSaft::new_polar(&["H2O", "NH3"]).unwrap();
    let fa = n_ab
        .residual_helmholtz_tv(320.0, 2.5e-5, &[0.4, 0.6])
        .unwrap();
    let fb = n_ba
        .residual_helmholtz_tv(320.0, 2.5e-5, &[0.6, 0.4])
        .unwrap();
    approx_eq(fa.f, fb.f, 1e-13, 1e-14, "PCP-SAFT NH3/H2O order F");
    approx_eq(fa.f_v, fb.f_v, 1e-13, 1e-9, "PCP-SAFT NH3/H2O order F_V");
}

#[test]
fn total_mole_scaling_invariance_down_to_1e_140_and_overflow_guard() {
    for model in [
        PcSaft::new(&["H2O", "MEOH"]).unwrap(),
        PcSaft::new_polar(&["CO2", "ACETONE"]).unwrap(),
        PcSaft::new_simplified_polar(&["NH3", "H2O"]).unwrap(),
    ] {
        let t = 340.0;
        let v_unit = 4.0e-5;
        let x = [0.6, 0.4];
        let ref_h = model.residual_helmholtz_tv(t, v_unit, &x).unwrap();

        for &scale in &[1.0e-12_f64, 1.0e-140_f64] {
            let v = v_unit * scale;
            let n = [x[0] * scale, x[1] * scale];
            let h = model.residual_helmholtz_tv(t, v, &n).unwrap();
            approx_eq(h.f / scale, ref_h.f, 1e-13, 1e-14, "F/N");
            approx_eq(h.f_t / scale, ref_h.f_t, 1e-13, 1e-14, "F_T/N");
            approx_eq(h.f_v, ref_h.f_v, 1e-13, 1e-10, "F_V");
            approx_eq(h.f_tt / scale, ref_h.f_tt, 1e-13, 1e-14, "F_TT/N");
            approx_eq(h.f_tv, ref_h.f_tv, 1e-13, 1e-10, "F_TV");
            approx_eq(h.f_vv * scale, ref_h.f_vv, 1e-13, 1e-6, "N*F_VV");
            for i in 0..2 {
                approx_eq(h.f_n[i], ref_h.f_n[i], 1e-13, 1e-13, "F_n");
                approx_eq(h.f_tn[i], ref_h.f_tn[i], 1e-13, 1e-13, "F_Tn");
                approx_eq(h.f_vn[i] * scale, ref_h.f_vn[i], 1e-13, 1e-9, "N*F_Vn");
                for j in 0..2 {
                    approx_eq(
                        h.f_nn[i * 2 + j] * scale,
                        ref_h.f_nn[i * 2 + j],
                        1e-13,
                        1e-12,
                        "N*F_nn",
                    );
                }
            }
            for (xk, xk_ref) in h.x_k.iter().zip(&ref_h.x_k) {
                approx_eq(*xk, *xk_ref, 1e-13, 1e-14, "X_k");
            }
        }

        // Genuine overflow scale (`1 / N` overflows `f64`) must return an explicit Error.
        let tiny = 1.0e-310_f64;
        assert!(
            model
                .residual_helmholtz_tv(t, v_unit * tiny, &[x[0] * tiny, x[1] * tiny])
                .is_err()
        );
    }
}

#[test]
fn zero_mole_cross_associating_infinite_dilution_matches_fd() {
    let model = PcSaft::new(&["H2O", "MEOH"]).unwrap();
    let t = 340.0;
    let v = model
        .volume_tp(t, 2.0e6, &[1.0, 0.0], Phase::Liquid)
        .unwrap();
    let h0 = model.residual_helmholtz_tv(t, v, &[1.0, 0.0]).unwrap();

    // Infinite-dilution methanol site fractions must be strictly below 1 (bonded to water).
    assert_eq!(h0.x_k.len(), 6);
    assert!(h0.x_k[4] < 0.05 && h0.x_k[4] > 0.0);
    assert!(h0.x_k[5] < 0.05 && h0.x_k[5] > 0.0);

    // Match one-sided finite difference in n_MEOH (resolves upstream Fortran 9.03-unit bug).
    let eps = 1.0e-7;
    let hp = model.residual_helmholtz_tv(t, v, &[1.0, eps]).unwrap();
    let hp2 = model
        .residual_helmholtz_tv(t, v, &[1.0, 2.0 * eps])
        .unwrap();
    let fd_fn_meoh = (-3.0 * h0.f + 4.0 * hp.f - hp2.f) / (2.0 * eps);
    approx_eq(
        h0.f_n[1],
        fd_fn_meoh,
        1e-6,
        1e-6,
        "infinite-dilution F_n(MEOH)",
    );
    // And site fractions are continuous at n_MEOH -> 0+.
    approx_eq(h0.x_k[4], hp.x_k[4], 1e-4, 1e-6, "X_k[4] continuity");
    approx_eq(h0.x_k[5], hp.x_k[5], 1e-4, 1e-6, "X_k[5] continuity");
}

#[test]
fn tp_state_derivatives_and_thermodynamic_identities() {
    let model = PcSaft::new_simplified_polar(&["NH3", "H2O"]).unwrap();
    let t = 320.0;
    let p = 3.0e6;
    let x = [0.4, 0.6];
    let st = model.state_tp(t, p, &x, Phase::Liquid).unwrap();

    // Central FD in T and P
    let ht = 1.0e-4;
    let hp = 1.0e2;
    let st_tp = model.state_tp(t + ht, p, &x, Phase::Liquid).unwrap();
    let st_tm = model.state_tp(t - ht, p, &x, Phase::Liquid).unwrap();
    let st_pp = model.state_tp(t, p + hp, &x, Phase::Liquid).unwrap();
    let st_pm = model.state_tp(t, p - hp, &x, Phase::Liquid).unwrap();

    approx_eq(
        st.specific_volume.dt,
        (st_tp.volume - st_tm.volume) / (2.0 * ht),
        1e-6,
        1e-12,
        "dv/dT",
    );
    approx_eq(
        st.specific_volume.dp,
        (st_pp.volume - st_pm.volume) / (2.0 * hp),
        1e-6,
        1e-15,
        "dv/dP",
    );
    approx_eq(
        st.enthalpy.dt,
        (st_tp.enthalpy.value - st_tm.enthalpy.value) / (2.0 * ht),
        1e-6,
        1e-6,
        "dH/dT == Cp",
    );
    approx_eq(
        st.entropy.dt,
        (st_tp.entropy.value - st_tm.entropy.value) / (2.0 * ht),
        1e-6,
        1e-8,
        "dS/dT == Cp/T",
    );
    approx_eq(
        st.entropy.dp,
        -st.specific_volume.dt,
        1e-12,
        1e-15,
        "Maxwell dS/dP == -dV/dT",
    );
    approx_eq(
        st.residual_enthalpy.dp,
        st.volume - t * st.specific_volume.dt,
        1e-12,
        1e-12,
        "H_P == V - T V_T",
    );
    approx_eq(
        st.ln_fugacity.dn[1],
        st.ln_fugacity.dn[2],
        1e-12,
        1e-12,
        "Maxwell dlnphi_0/dn_1 == dlnphi_1/dn_0",
    );
    // Gibbs-Duhem: sum_i x_i * dlnphi_i/dn_j == 0
    for j in 0..2 {
        let gd = x[0] * st.ln_fugacity.dn[j] + x[1] * st.ln_fugacity.dn[2 + j];
        approx_eq(gd, 0.0, 1e-10, 1e-10, &format!("Gibbs-Duhem j={j}"));
    }
    assert!(
        st.specific_volume.dp < 0.0,
        "mechanical stability dV/dP < 0"
    );
    assert!(st.cp > st.cv && st.cv > 0.0 && st.speed_of_sound > 500.0);
}

#[test]
fn spcsaft_rdf_matches_finite_differences() {
    let model = PcSaft::new_simplified(&["H2O", "MEOH"]).unwrap();
    let t = 340.0;
    let v = 3.5e-5;
    let n = [0.6, 0.4];
    let r = spcsaft_rdf(&model, t, v, &n).unwrap();
    let ht = 1.0e-5 * t;
    let hv = 1.0e-5 * v;
    let rp_t = spcsaft_rdf(&model, t + ht, v, &n).unwrap();
    let rm_t = spcsaft_rdf(&model, t - ht, v, &n).unwrap();
    let rp_v = spcsaft_rdf(&model, t, v + hv, &n).unwrap();
    let rm_v = spcsaft_rdf(&model, t, v - hv, &n).unwrap();

    approx_eq(
        r.temperature,
        (rp_t.g - rm_t.g) / (2.0 * ht),
        1e-6,
        1e-10,
        "g_T",
    );
    approx_eq(r.volume, (rp_v.g - rm_v.g) / (2.0 * hv), 1e-6, 1e-6, "g_V");
    approx_eq(
        r.temperature_second,
        (rp_t.temperature - rm_t.temperature) / (2.0 * ht),
        1e-6,
        1e-10,
        "g_TT",
    );
    approx_eq(
        r.temperature_volume,
        (rp_t.volume - rm_t.volume) / (2.0 * ht),
        1e-6,
        1e-6,
        "g_TV",
    );
    approx_eq(
        r.volume_second,
        (rp_v.volume - rm_v.volume) / (2.0 * hv),
        1e-6,
        1e-2,
        "g_VV",
    );
    for i in 0..2 {
        let hn = 1.0e-5 * n[i];
        let mut np = n;
        let mut nm = n;
        np[i] += hn;
        nm[i] -= hn;
        let rp = spcsaft_rdf(&model, t, v, &np).unwrap();
        let rm = spcsaft_rdf(&model, t, v, &nm).unwrap();
        approx_eq(r.moles[i], (rp.g - rm.g) / (2.0 * hn), 1e-6, 1e-9, "g_n");
        approx_eq(
            r.temperature_moles[i],
            (rp.temperature - rm.temperature) / (2.0 * hn),
            1e-6,
            1e-9,
            "g_Tn",
        );
        approx_eq(
            r.volume_moles[i],
            (rp.volume - rm.volume) / (2.0 * hn),
            1e-6,
            1e-6,
            "g_Vn",
        );
        for j in 0..2 {
            approx_eq(
                r.moles_second[i * 2 + j],
                (rp.moles[j] - rm.moles[j]) / (2.0 * hn),
                1e-6,
                1e-9,
                "g_nn",
            );
        }
    }
}

#[test]
fn multipol_subterm_controls_and_fast_skip() {
    let mut model = PcSaft::new_polar(&["CO2", "ACETONE"]).unwrap();
    let t = 320.0;
    let v = 6.0e-5;
    let n = [0.4, 0.6];

    let full = model.multipol_residual_helmholtz_tv(t, v, &n).unwrap();
    assert!(full.f.abs() > 1.0e-2);

    model.set_multipol_control(true, false, false).unwrap();
    let dd = model.multipol_residual_helmholtz_tv(t, v, &n).unwrap();
    model.set_multipol_control(false, true, false).unwrap();
    let qq = model.multipol_residual_helmholtz_tv(t, v, &n).unwrap();
    model.set_multipol_control(false, false, true).unwrap();
    let dq = model.multipol_residual_helmholtz_tv(t, v, &n).unwrap();

    assert!(dd.f.abs() > 1.0e-3);
    assert!(qq.f.abs() > 1.0e-3);
    assert!(dq.f.abs() > 1.0e-4);
    approx_eq(dd.f + qq.f + dq.f, full.f, 1e-13, 1e-14, "DD+QQ+DQ == full");

    // Fast-skip when all polar terms are disabled or moments are zero.
    model.set_multipol_control(false, false, false).unwrap();
    assert!(!model.has_polar());
    let none = model.multipol_residual_helmholtz_tv(t, v, &n).unwrap();
    approx_eq(none.f, 0.0, 0.0, 1e-15, "disabled multipol == 0");

    let non_polar = PcSaft::with_options(
        &["C1", "C3"],
        PcSaftOptions {
            kind: PcSaftKind::PcpSaft,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(!non_polar.has_polar());
    let zero_mp = non_polar
        .multipol_residual_helmholtz_tv(300.0, 1.0e-4, &[0.5, 0.5])
        .unwrap();
    approx_eq(zero_mp.f, 0.0, 0.0, 1e-15, "non-polar fluid multipol == 0");
}

#[test]
fn pure_fluid_saturation_for_associating_water_and_polar_co2() {
    for (model, t) in [
        (PcSaft::new(&["H2O"]).unwrap(), 373.15),
        (PcSaft::new(&["MEOH"]).unwrap(), 337.85),
        (PcSaft::new_polar(&["CO2"]).unwrap(), 260.0),
    ] {
        let prep = model.prepare(t, &[1.0]).unwrap();
        let sat = prep.saturation_pressure(t).unwrap();
        assert!(sat.pressure > 0.0 && sat.liquid_volume < sat.vapor_volume);
        assert!(sat.enthalpy_vaporization > 0.0 && sat.entropy_vaporization > 0.0);

        let sl = prep.state(1.0 / sat.liquid_volume).unwrap();
        let sv = prep.state(1.0 / sat.vapor_volume).unwrap();
        approx_eq(sl.pressure, sat.pressure, 1e-7, 1.0, "sat liquid P");
        approx_eq(sv.pressure, sat.pressure, 1e-7, 1.0, "sat vapor P");
        approx_eq(sl.ln_phi[0], sv.ln_phi[0], 1e-8, 1e-9, "sat fugacity eq");
    }
}

#[test]
fn association_site_scheme_utilities() {
    assert_eq!(AssocScheme::NoAssoc.num_sites(), 0);
    assert_eq!(AssocScheme::Scheme1.num_sites(), 1);
    assert_eq!(AssocScheme::Scheme1ea.num_sites(), 1);
    assert_eq!(AssocScheme::Scheme2B.num_sites(), 2);
    assert_eq!(AssocScheme::Scheme3B.num_sites(), 3);
    assert_eq!(AssocScheme::Scheme4C.num_sites(), 4);
    let empty = association::AssociationParams::empty(2);
    assert!(!empty.has_sites());
}

#[test]
fn induced_association_zero_mole_absent_partner_matches_limit_and_fd() {
    let mut model = PcSaft::new(&["CO2", "H2O"]).unwrap();
    model
        .set_assoc_params(0, 10_000.0, 0.03, AssocScheme::Scheme1ea)
        .unwrap();

    let t = 340.0;
    let v = 1.0e-4;
    let f0 = model.residual_helmholtz_tv(t, v, &[1.0, 0.0]).unwrap();
    let f_lim = model.residual_helmholtz_tv(t, v, &[1.0, 1.0e-8]).unwrap();

    // Present CO2 has no self-association (`Scheme1ea` cannot bond to itself),
    // while absent H2O sites still cross-associate with present CO2.
    approx_eq(f0.x_k[0], 1.0, 1e-12, 1e-13, "present CO2 X_k");
    assert!(
        f0.x_k[1..].iter().any(|&xk| xk < 0.8),
        "expected induced cross-association on absent H2O sites, got {:?}",
        f0.x_k
    );
    for (k, (&xk0, &xkl)) in f0.x_k.iter().zip(f_lim.x_k.iter()).enumerate() {
        approx_eq(xk0, xkl, 1e-5, 1e-6, &format!("X_k[{k}] limit"));
    }

    // Verify F_n(H2O), F_Tn(H2O), F_Vn(H2O), and F_nn(*, H2O) against one-sided FD.
    let h = 1.0e-7;
    let f_p = model.residual_helmholtz_tv(t, v, &[1.0, h]).unwrap();
    let f_2p = model.residual_helmholtz_tv(t, v, &[1.0, 2.0 * h]).unwrap();

    let fn1_fd = (-3.0 * f0.f + 4.0 * f_p.f - f_2p.f) / (2.0 * h);
    let ftn1_fd = (-3.0 * f0.f_t + 4.0 * f_p.f_t - f_2p.f_t) / (2.0 * h);
    let fvn1_fd = (-3.0 * f0.f_v + 4.0 * f_p.f_v - f_2p.f_v) / (2.0 * h);
    let fnn01_fd = (-3.0 * f0.f_n[0] + 4.0 * f_p.f_n[0] - f_2p.f_n[0]) / (2.0 * h);
    let fnn11_fd = (-3.0 * f0.f_n[1] + 4.0 * f_p.f_n[1] - f_2p.f_n[1]) / (2.0 * h);

    approx_eq(f0.f_n[1], fn1_fd, 1e-6, 1e-7, "induced assoc F_n(H2O)");
    approx_eq(f0.f_tn[1], ftn1_fd, 1e-5, 1e-7, "induced assoc F_Tn(H2O)");
    approx_eq(f0.f_vn[1], fvn1_fd, 1e-5, 1e-4, "induced assoc F_Vn(H2O)");
    approx_eq(
        f0.f_nn[1],
        fnn01_fd,
        1e-5,
        1e-6,
        "induced assoc F_nn(CO2,H2O)",
    );
    approx_eq(
        f0.f_nn[3],
        fnn11_fd,
        1e-5,
        1e-6,
        "induced assoc F_nn(H2O,H2O)",
    );
}

#[test]
fn multipol_dq_zero_mole_boundary_matches_fd_and_rejects_ambiguous_mixed_boundary() {
    let mut model = PcSaft::new_polar(&["ACETONE", "CO2"]).unwrap();
    model.set_multipol_control(false, false, true).unwrap();

    let t = 340.0;
    let v = 1.0e-4;
    let h = 1.0e-7;

    // Test both boundary faces: absent quadrupole [1, 0] and absent dipole [0, 1].
    for (n0, absent, present, label) in [
        ([1.0, 0.0], 1usize, 0usize, "absent_CO2"),
        ([0.0, 1.0], 0usize, 1usize, "absent_ACETONE"),
    ] {
        let f0 = model.multipol_residual_helmholtz_tv(t, v, &n0).unwrap();
        let mut n1 = n0;
        let mut n2 = n0;
        n1[absent] = h;
        n2[absent] = 2.0 * h;
        let f1 = model.multipol_residual_helmholtz_tv(t, v, &n1).unwrap();
        let f2 = model.multipol_residual_helmholtz_tv(t, v, &n2).unwrap();

        let fn_abs_fd = (-3.0 * f0.f + 4.0 * f1.f - f2.f) / (2.0 * h);
        let ftn_abs_fd = (-3.0 * f0.f_t + 4.0 * f1.f_t - f2.f_t) / (2.0 * h);
        let fvn_abs_fd = (-3.0 * f0.f_v + 4.0 * f1.f_v - f2.f_v) / (2.0 * h);
        let fnn_pa_fd =
            (-3.0 * f0.f_n[present] + 4.0 * f1.f_n[present] - f2.f_n[present]) / (2.0 * h);
        let fnn_aa_fd = (-3.0 * f0.f_n[absent] + 4.0 * f1.f_n[absent] - f2.f_n[absent]) / (2.0 * h);

        approx_eq(f0.f, 0.0, 0.0, 1e-15, &format!("{label} F == 0"));
        assert!(
            f0.f_n[absent].abs() > 1.0e-3,
            "{label}: expected non-trivial DQ Pade boundary derivative"
        );
        approx_eq(
            f0.f_n[absent],
            fn_abs_fd,
            1e-6,
            1e-8,
            &format!("{label} F_n"),
        );
        approx_eq(
            f0.f_tn[absent],
            ftn_abs_fd,
            1e-5,
            1e-8,
            &format!("{label} F_Tn"),
        );
        approx_eq(
            f0.f_vn[absent],
            fvn_abs_fd,
            1e-5,
            1e-4,
            &format!("{label} F_Vn"),
        );
        approx_eq(
            f0.f_nn[present * 2 + absent],
            fnn_pa_fd,
            1e-5,
            1e-7,
            &format!("{label} F_nn(present,absent)"),
        );
        approx_eq(
            f0.f_nn[absent * 2 + absent],
            fnn_aa_fd,
            1e-5,
            1e-6,
            &format!("{label} F_nn(absent,absent)"),
        );
    }

    // Explicit reference magnitude check for ACETONE + absent CO2 (proves Pade
    // factor a2'^2 / (a2' - a3') is preserved rather than collapsing to a2').
    let f_ref = model
        .multipol_residual_helmholtz_tv(t, v, &[1.0, 0.0])
        .unwrap();
    approx_eq(
        f_ref.f_n[1],
        -0.861_517_502_5,
        1e-6,
        1e-6,
        "ACETONE/absent-CO2 DQ F_n reference",
    );

    // When two species in the absent polar family are simultaneously zero, the
    // mixed DQ Pade boundary second derivative is path-dependent and must error.
    let mut ambig_q = PcSaft::new_polar(&["ACETONE", "CO2", "N2"]).unwrap();
    ambig_q.set_polar_moments(2, 0.0, 1.5).unwrap();
    ambig_q.set_multipol_control(false, false, true).unwrap();
    let err_q = ambig_q
        .multipol_residual_helmholtz_tv(t, v, &[1.0, 0.0, 0.0])
        .unwrap_err();
    assert!(matches!(err_q, Error::InvalidInput(_)));

    let mut ambig_d = PcSaft::new_polar(&["CO2", "ACETONE", "NH3"]).unwrap();
    ambig_d.set_polar_moments(2, 1.47, 0.0).unwrap();
    ambig_d.set_multipol_control(false, false, true).unwrap();
    let err_d = ambig_d
        .multipol_residual_helmholtz_tv(t, v, &[1.0, 0.0, 0.0])
        .unwrap_err();
    assert!(matches!(err_d, Error::InvalidInput(_)));

    // When both polar families are absent in a non-polar solvent, a_2 = O(eps^2)
    // and a_3 = O(eps^3), so the boundary is smooth and matches FD.
    let mut both_absent = PcSaft::new_polar(&["ACETONE", "CO2", "C1"]).unwrap();
    both_absent
        .set_multipol_control(false, false, true)
        .unwrap();
    let f_ba = both_absent
        .multipol_residual_helmholtz_tv(t, v, &[0.0, 0.0, 1.0])
        .unwrap();
    let f_ba_h = both_absent
        .multipol_residual_helmholtz_tv(t, v, &[h, h, 1.0])
        .unwrap();
    approx_eq(
        f_ba.f_nn[1],
        f_ba_h.f / (h * h),
        1e-4,
        1e-6,
        "both polar families absent F_nn(ACETONE,CO2)",
    );
}

#[test]
fn zero_mole_state_tp_partial_molar_ideal_entropy_uses_thermopack_cutoff() {
    let mix = PcSaft::new(&["C1", "C2"]).unwrap();
    let pure = PcSaft::new(&["C1"]).unwrap();

    let t = 340.0;
    let p = 1.0e5;
    let s_mix = mix.state_tp(t, p, &[1.0, 0.0], Phase::Vapor).unwrap();
    let s_pure = pure.state_tp(t, p, &[1.0], Phase::Vapor).unwrap();

    // Mixture entropy preserves exact 0 * ln(0) = 0 matching pure C1.
    approx_eq(
        s_mix.entropy.value,
        s_pure.entropy.value,
        1e-12,
        1e-12,
        "mixture entropy at [1,0] == pure C1",
    );
    // Partial molar entropy of absent C2 uses ThermoPack's 1e-100 floor
    // (`s_id_i = s_id_pure_i - R * ln(1e-100)`), staying finite (~2151.24 J/(mol K)).
    assert!(s_mix.entropy.dn[1].is_finite());
    approx_eq(
        s_mix.entropy.dn[1],
        2_151.242_215_167,
        1e-6,
        1e-4,
        "zero-mole partial molar entropy matches ThermoPack cutoff",
    );

    // Fugacity coefficients remain finite and exact (`ln_phi_i = F_n_i - ln(Z)`).
    let f = mix
        .residual_helmholtz_tv(t, s_mix.specific_volume.value, &[1.0, 0.0])
        .unwrap();
    let ln_z = s_mix.compressibility_factor.value.ln();
    for i in 0..2 {
        approx_eq(
            s_mix.ln_fugacity.ln_phi[i],
            f.f_n[i] - ln_z,
            1e-12,
            1e-13,
            &format!("zero-mole ln_phi[{i}]"),
        );
    }
}

#[test]
fn set_kij_validates_inputs_and_is_atomic_on_error() {
    let mut model = PcSaft::new_polar(&["CO2", "ACETONE"]).unwrap();
    let orig_kij = [
        model.get_kij(0, 0),
        model.get_kij(0, 1),
        model.get_kij(1, 0),
        model.get_kij(1, 1),
    ];
    let orig_eps = model.multipol().unwrap().eps_divk_ij.clone();

    for (i, j, val) in [
        (0, 0, 1.0),
        (0, 0, 0.1),
        (0, 1, 1.0),
        (0, 1, 1.5),
        (0, 2, 0.0),
        (2, 0, 0.0),
        (0, 1, f64::NAN),
        (0, 1, f64::INFINITY),
        (0, 1, f64::NEG_INFINITY),
    ] {
        let err = model.set_kij(i, j, val).unwrap_err();
        assert!(
            matches!(err, Error::InvalidInput(_)),
            "expected InvalidInput for set_kij({i}, {j}, {val}), got {err:?}"
        );
        assert_eq!(model.get_kij(0, 0), orig_kij[0]);
        assert_eq!(model.get_kij(0, 1), orig_kij[1]);
        assert_eq!(model.get_kij(1, 0), orig_kij[2]);
        assert_eq!(model.get_kij(1, 1), orig_kij[3]);
        assert_eq!(model.multipol().unwrap().eps_divk_ij, orig_eps);
    }

    model.set_kij(0, 1, 0.05).unwrap();
    assert_eq!(model.get_kij(0, 1), 0.05);
    assert_eq!(model.get_kij(1, 0), 0.05);
    assert_ne!(model.multipol().unwrap().eps_divk_ij[1], orig_eps[1]);
    assert_eq!(
        model.multipol().unwrap().eps_divk_ij[1],
        model.multipol().unwrap().eps_divk_ij[2]
    );
}

#[test]
fn pure_fluid_saturation_converges_for_polar_meoh_low_t_water_and_rejects_supercritical() {
    for (comp, t, expected_p, rtol_p) in [
        ("MEOH", 450.0, 444_367.074_811, 1e-5),
        ("MEOH", 500.0, 1_227_151.027, 1e-4),
        ("H2O", 280.0, 991.643, 1e-3),
    ] {
        let model = PcSaft::new_polar(&[comp]).unwrap();
        let prep = model.prepare(t, &[1.0]).unwrap();
        let sat = prep.saturation_pressure(t).unwrap();

        assert!(sat.pressure > 0.0 && sat.liquid_volume < sat.vapor_volume);
        assert!(sat.enthalpy_vaporization > 0.0 && sat.entropy_vaporization > 0.0);
        approx_eq(
            sat.pressure,
            expected_p,
            rtol_p,
            1.0,
            &format!("{comp}@{t}K P_sat"),
        );

        let sl = prep.state(1.0 / sat.liquid_volume).unwrap();
        let sv = prep.state(1.0 / sat.vapor_volume).unwrap();
        approx_eq(
            sl.pressure,
            sat.pressure,
            1e-6,
            0.0,
            &format!("{comp}@{t}K liquid P"),
        );
        approx_eq(
            sv.pressure,
            sat.pressure,
            1e-6,
            0.0,
            &format!("{comp}@{t}K vapor P"),
        );
        let rt = RGAS * t;
        let zl = sat.pressure * sat.liquid_volume / rt;
        let zv = sat.pressure * sat.vapor_volume / rt;
        let ln_phi_l_gibbs = sl.residual_helmholtz + (zl - 1.0) - zl.ln();
        let ln_phi_v_gibbs = sv.residual_helmholtz + (zv - 1.0) - zv.ln();
        approx_eq(
            ln_phi_l_gibbs,
            ln_phi_v_gibbs,
            1e-10,
            1e-10,
            &format!("{comp}@{t}K Gibbs fugacity eq"),
        );
        approx_eq(
            sl.ln_phi[0],
            sv.ln_phi[0],
            1e-8,
            1e-8,
            &format!("{comp}@{t}K state fugacity eq"),
        );
    }

    // Supercritical CO2 has a monotone isotherm and must return NotConverged.
    let co2 = PcSaft::new_polar(&["CO2"]).unwrap();
    let prep = co2.prepare(400.0, &[1.0]).unwrap();
    let err = prep.saturation_pressure(400.0).unwrap_err();
    assert!(matches!(err, Error::NotConverged(_)));
}

#[test]
fn state_tp_raw_non_unit_moles_scaling_and_fd_derivatives() {
    type ScalingCase<'a> = (
        &'a str,
        PcSaftKind,
        &'a [&'a str],
        f64,
        f64,
        &'a [f64],
        Phase,
    );
    let cases: &[ScalingCase<'_>] = &[
        (
            "sPC-SAFT_C1_C2_C3_N2",
            PcSaftKind::SPcSaft,
            &["C1", "C2", "C3"],
            320.0,
            4.0e6,
            &[1.2, 0.5, 0.3],
            Phase::Vapor,
        ),
        (
            "sPC-SAFT_H2O_MEOH_C1_N2",
            PcSaftKind::SPcSaft,
            &["H2O", "MEOH", "C1"],
            360.0,
            6.0e6,
            &[0.8, 0.6, 0.6],
            Phase::Liquid,
        ),
        (
            "sPCP-SAFT_NH3_H2O_1e-140",
            PcSaftKind::SPcpSaft,
            &["NH3", "H2O"],
            320.0,
            3.0e6,
            &[0.4e-140, 0.6e-140],
            Phase::Liquid,
        ),
    ];

    for &(label, kind, comps, t, p, moles, phase) in cases {
        let model = PcSaft::new_with_options(kind, comps, "default").unwrap();
        let nc = comps.len();
        let total: f64 = moles.iter().sum();
        let x: Vec<f64> = moles.iter().map(|&v| v / total).collect();

        let s_raw = model.state_tp(t, p, moles, phase).unwrap();
        let s_unit = model.state_tp(t, p, &x, phase).unwrap();

        // Intensive quantities and partial molar properties are invariant in total moles.
        approx_eq(
            s_raw.specific_volume.value,
            s_unit.specific_volume.value,
            1e-12,
            1e-15,
            &format!("{label} v"),
        );
        approx_eq(
            s_raw.compressibility_factor.value,
            s_unit.compressibility_factor.value,
            1e-12,
            1e-15,
            &format!("{label} Z"),
        );
        approx_eq(
            s_raw.enthalpy.value,
            s_unit.enthalpy.value,
            1e-12,
            1e-12,
            &format!("{label} h"),
        );
        approx_eq(
            s_raw.entropy.value,
            s_unit.entropy.value,
            1e-12,
            1e-12,
            &format!("{label} s"),
        );
        approx_eq(
            s_raw.speed_of_sound,
            s_unit.speed_of_sound,
            1e-12,
            1e-12,
            &format!("{label} sos"),
        );

        for i in 0..nc {
            approx_eq(
                s_raw.specific_volume.dn[i],
                s_unit.specific_volume.dn[i],
                1e-12,
                1e-15,
                &format!("{label} v_bar[{i}]"),
            );
            approx_eq(
                s_raw.enthalpy.dn[i],
                s_unit.enthalpy.dn[i],
                1e-12,
                1e-12,
                &format!("{label} h_bar[{i}]"),
            );
            approx_eq(
                s_raw.entropy.dn[i],
                s_unit.entropy.dn[i],
                1e-12,
                1e-12,
                &format!("{label} s_bar[{i}]"),
            );
            approx_eq(
                s_raw.ln_fugacity.ln_phi[i],
                s_unit.ln_fugacity.ln_phi[i],
                1e-12,
                1e-13,
                &format!("{label} ln_phi[{i}]"),
            );
            // Intensive mole derivatives scale strictly as 1/N.
            approx_eq(
                s_raw.compressibility_factor.dn[i] * total,
                s_unit.compressibility_factor.dn[i],
                1e-12,
                1e-13,
                &format!("{label} Z_dn[{i}] * N"),
            );
            for j in 0..nc {
                approx_eq(
                    s_raw.ln_fugacity.dn[i * nc + j] * total,
                    s_unit.ln_fugacity.dn[i * nc + j],
                    1e-12,
                    1e-13,
                    &format!("{label} ln_phi_dn[{i},{j}] * N"),
                );
            }
        }

        // Direct central FD verification at raw non-unit moles for N = 2.0 cases.
        if total > 0.5 {
            for j in 0..nc {
                let h = 1.0e-6 * moles[j].max(1.0);
                let mut np = moles.to_vec();
                let mut nm = moles.to_vec();
                np[j] += h;
                nm[j] -= h;
                let sp = model.state_tp(t, p, &np, phase).unwrap();
                let sm = model.state_tp(t, p, &nm, phase).unwrap();

                let z_dn_fd =
                    (sp.compressibility_factor.value - sm.compressibility_factor.value) / (2.0 * h);
                approx_eq(
                    s_raw.compressibility_factor.dn[j],
                    z_dn_fd,
                    1e-6,
                    1e-8,
                    &format!("{label} FD Z_dn[{j}]"),
                );

                // Partial molar volume relation: dv_molar/dn_j = (V_bar_j - v_molar) / N.
                let v_dn_fd = (sp.specific_volume.value - sm.specific_volume.value) / (2.0 * h);
                let v_dn_analytic =
                    (s_raw.specific_volume.dn[j] - s_raw.specific_volume.value) / total;
                approx_eq(
                    v_dn_analytic,
                    v_dn_fd,
                    1e-6,
                    1e-10,
                    &format!("{label} FD dv_molar/dn[{j}]"),
                );

                for i in 0..nc {
                    let lnp_dn_fd =
                        (sp.ln_fugacity.ln_phi[i] - sm.ln_fugacity.ln_phi[i]) / (2.0 * h);
                    approx_eq(
                        s_raw.ln_fugacity.dn[i * nc + j],
                        lnp_dn_fd,
                        1e-6,
                        1e-8,
                        &format!("{label} FD ln_phi_dn[{i},{j}]"),
                    );
                }
            }
        }
    }
}

#[test]
fn multipol_dq_positive_trace_sweep_is_smooth_across_fres_fn_and_full() {
    let mut model = PcSaft::new_polar(&["ACETONE", "CO2"]).unwrap();
    model.set_multipol_control(false, false, true).unwrap();
    let mp = model.multipol().unwrap();

    let t = 340.0;
    let v = 1.0e-4;
    let f_zero = model
        .multipol_residual_helmholtz_tv(t, v, &[1.0, 0.0])
        .unwrap();
    approx_eq(
        f_zero.f_n[1],
        -0.861_517_495,
        1e-6,
        1e-6,
        "DQ F_n(CO2) at n_CO2=0",
    );
    approx_eq(
        f_zero.f_nn[3],
        -1.417_187,
        1e-4,
        1e-4,
        "DQ F_nn(CO2,CO2) at n_CO2=0",
    );

    for n_co2 in [0.0, 1.0e-30, 1.0e-16, 1.0e-15, 1.0e-12, 1.0e-10, 1.0e-8] {
        let n = [1.0, n_co2];
        let val = mp.fres_value(t, v, &n).unwrap();
        let jet = mp.fres_pressure_jet(t, v, &n).unwrap();
        let fn_only = mp.fres_fn(t, v, &n).unwrap();
        let full = model.multipol_residual_helmholtz_tv(t, v, &n).unwrap();

        approx_eq(val, full.f, 1e-13, 1e-15, &format!("DQ val n={n_co2:e}"));
        approx_eq(
            jet.f,
            full.f,
            1e-13,
            1e-15,
            &format!("DQ jet.f n={n_co2:e}"),
        );
        approx_eq(
            jet.f_v,
            full.f_v,
            1e-13,
            1e-15,
            &format!("DQ jet.f_v n={n_co2:e}"),
        );
        for i in 0..2 {
            approx_eq(
                fn_only[i],
                full.f_n[i],
                1e-13,
                1e-15,
                &format!("DQ fres_fn[{i}] vs full n={n_co2:e}"),
            );
        }

        // Continuity relative to the exact n_CO2 = 0 boundary (no 1e-15 clamp jump).
        let tol = if n_co2 <= 1.0e-12 { 1e-10 } else { 1e-5 };
        approx_eq(
            full.f_n[1],
            f_zero.f_n[1],
            tol,
            tol,
            &format!("DQ smooth F_n(CO2) n={n_co2:e}"),
        );
        approx_eq(
            full.f_tn[1],
            f_zero.f_tn[1],
            tol,
            tol,
            &format!("DQ smooth F_Tn(CO2) n={n_co2:e}"),
        );
        approx_eq(
            full.f_vn[1],
            f_zero.f_vn[1],
            tol,
            tol * 1.0e4,
            &format!("DQ smooth F_Vn(CO2) n={n_co2:e}"),
        );
        approx_eq(
            full.f_nn[1],
            f_zero.f_nn[1],
            tol,
            tol,
            &format!("DQ smooth F_nn(0,1) n={n_co2:e}"),
        );
        approx_eq(
            full.f_nn[3],
            f_zero.f_nn[3],
            tol,
            tol,
            &format!("DQ smooth F_nn(1,1) n={n_co2:e}"),
        );
    }

    // DD and QQ ultra-trace sweep (including 1e-16 and 1e-140) stays finite and smooth.
    for (comps, dd, qq) in [
        (&["C1", "ACETONE"][..], true, false),
        (&["C1", "CO2"][..], false, true),
    ] {
        let mut m = PcSaft::new_polar(comps).unwrap();
        m.set_multipol_control(dd, qq, false).unwrap();
        let mp_sub = m.multipol().unwrap();
        let f0 = m.multipol_residual_helmholtz_tv(t, v, &[1.0, 0.0]).unwrap();
        for n_tr in [1.0e-140, 1.0e-30, 1.0e-16, 1.0e-15, 1.0e-12] {
            let n = [1.0, n_tr];
            let fn_only = mp_sub.fres_fn(t, v, &n).unwrap();
            let full = m.multipol_residual_helmholtz_tv(t, v, &n).unwrap();
            approx_eq(
                fn_only[1],
                full.f_n[1],
                1e-12,
                1e-15,
                &format!("DD/QQ fres_fn n={n_tr:e}"),
            );
            approx_eq(
                full.f_nn[3],
                f0.f_nn[3],
                1e-8,
                1e-10,
                &format!("DD/QQ smooth F_nn(1,1) n={n_tr:e}"),
            );
        }
    }
}

#[test]
fn multipol_dq_identical_trace_species_boundary_matches_binary_diagonal() {
    let t = 340.0;
    let v = 1.0e-4;

    let mut bin = PcSaft::new_polar(&["ACETONE", "CO2"]).unwrap();
    bin.set_multipol_control(false, false, true).unwrap();
    let f_bin = bin
        .multipol_residual_helmholtz_tv(t, v, &[1.0, 0.0])
        .unwrap();

    let mut dup_q = PcSaft::new_polar(&["ACETONE", "CO2", "CO2"]).unwrap();
    dup_q.set_multipol_control(false, false, true).unwrap();

    for n in [[1.0, 0.0, 0.0], [1.0, 1.0e-16, 0.0]] {
        let f_dup = dup_q.multipol_residual_helmholtz_tv(t, v, &n).unwrap();
        approx_eq(
            f_dup.f_n[1],
            f_bin.f_n[1],
            1e-12,
            1e-13,
            &format!("dup Q F_n[1] at {n:?}"),
        );
        approx_eq(
            f_dup.f_n[2],
            f_bin.f_n[1],
            1e-12,
            1e-13,
            &format!("dup Q F_n[2] at {n:?}"),
        );
        approx_eq(
            f_dup.f_nn[4],
            f_bin.f_nn[3],
            1e-10,
            1e-10,
            &format!("dup Q F_nn(1,1) at {n:?}"),
        );
        approx_eq(
            f_dup.f_nn[5],
            f_bin.f_nn[3],
            1e-10,
            1e-10,
            &format!("dup Q F_nn(1,2) at {n:?}"),
        );
        approx_eq(
            f_dup.f_nn[8],
            f_bin.f_nn[3],
            1e-10,
            1e-10,
            &format!("dup Q F_nn(2,2) at {n:?}"),
        );
        approx_eq(
            f_dup.f_nn[5],
            -1.417_187,
            1e-4,
            1e-4,
            &format!("dup Q F_nn(1,2) reference at {n:?}"),
        );
    }

    // Symmetric check: duplicate absent dipolar species with present quadrupole.
    let f_bin_d = bin
        .multipol_residual_helmholtz_tv(t, v, &[0.0, 1.0])
        .unwrap();
    let mut dup_d = PcSaft::new_polar(&["CO2", "ACETONE", "ACETONE"]).unwrap();
    dup_d.set_multipol_control(false, false, true).unwrap();
    for n in [[1.0, 0.0, 0.0], [1.0, 1.0e-16, 0.0]] {
        let f_dup = dup_d.multipol_residual_helmholtz_tv(t, v, &n).unwrap();
        approx_eq(
            f_dup.f_n[1],
            f_bin_d.f_n[0],
            1e-12,
            1e-13,
            &format!("dup D F_n[1] at {n:?}"),
        );
        approx_eq(
            f_dup.f_nn[5],
            f_bin_d.f_nn[0],
            1e-10,
            1e-10,
            &format!("dup D F_nn(1,2) at {n:?}"),
        );
    }
}

#[test]
fn near_critical_methane_saturation_and_ultra_low_t_water_rejection() {
    // Near-critical methane: narrow van der Waals loop stepped over by coarse grid.
    let c1 = PcSaft::new(&["C1"]).unwrap();
    let t_nc = 191.399_581_288_3;
    let prep = c1.prepare(t_nc, &[1.0]).unwrap();
    let sat = prep.saturation_pressure(t_nc).unwrap();

    approx_eq(sat.pressure, 4_674_933.015_24, 1e-6, 1.0, "near-crit C1 P");
    approx_eq(
        sat.liquid_volume,
        1.076_017_99e-4,
        1e-4,
        1e-9,
        "near-crit C1 Vl",
    );
    approx_eq(
        sat.vapor_volume,
        1.091_267_66e-4,
        1e-4,
        1e-9,
        "near-crit C1 Vv",
    );
    assert!(sat.liquid_volume < sat.vapor_volume);

    let sl = prep.state(1.0 / sat.liquid_volume).unwrap();
    let sv = prep.state(1.0 / sat.vapor_volume).unwrap();
    approx_eq(
        sl.pressure,
        sat.pressure,
        1e-6,
        0.0,
        "near-crit C1 liquid P",
    );
    approx_eq(sv.pressure, sat.pressure, 1e-6, 0.0, "near-crit C1 vapor P");
    approx_eq(
        sl.ln_phi[0],
        sv.ln_phi[0],
        1e-9,
        1e-9,
        "near-crit C1 fugacity",
    );

    // Ultra-low-T PCP-SAFT water below f64 liquid pressure resolution must return NotConverged.
    let h2o = PcSaft::new_polar(&["H2O"]).unwrap();
    for t_low in [150.0, 100.0] {
        let p_low = h2o.prepare(t_low, &[1.0]).unwrap();
        let err = p_low.saturation_pressure(t_low).unwrap_err();
        assert!(
            matches!(err, Error::NotConverged(_)),
            "expected NotConverged for PCP-SAFT H2O at {t_low} K, got {err:?}"
        );
    }
}

#[test]
fn elliot_precondition_enforced_at_build_time_and_atomic_on_error() {
    let t = 340.0;
    let v = 1.0e-4;
    let n = [0.4, 0.6];

    // Order (a): Scheme1ea already set, then switching to Elliot must fail atomically.
    let mut co2_h2o = PcSaft::new(&["CO2", "H2O"]).unwrap();
    co2_h2o
        .set_assoc_params(0, 10_000.0, 0.03, AssocScheme::Scheme1ea)
        .unwrap();
    let before_a = co2_h2o.residual_helmholtz_tv(t, v, &n).unwrap();
    let err_a = co2_h2o
        .set_delta_comb_rule(DeltaCombRule::Elliot)
        .unwrap_err();
    assert!(matches!(err_a, Error::Unsupported(_)));
    assert_eq!(
        co2_h2o.association().delta_comb_rule,
        DeltaCombRule::Standard
    );
    let after_a = co2_h2o.residual_helmholtz_tv(t, v, &n).unwrap();
    approx_eq(
        after_a.f,
        before_a.f,
        0.0,
        1e-15,
        "Order (a) model preserved",
    );

    // Order (b): Elliot already active, then switching component 0 to Scheme1ea must fail atomically.
    let mut nh3_h2o = PcSaft::new(&["NH3", "H2O"]).unwrap();
    nh3_h2o.set_delta_comb_rule(DeltaCombRule::Elliot).unwrap();
    let before_b = nh3_h2o.residual_helmholtz_tv(t, v, &n).unwrap();
    let err_b = nh3_h2o
        .set_assoc_params(0, 10_000.0, 0.03, AssocScheme::Scheme1ea)
        .unwrap_err();
    assert!(matches!(err_b, Error::Unsupported(_)));
    assert_eq!(nh3_h2o.association().delta_comb_rule, DeltaCombRule::Elliot);
    let after_b = nh3_h2o.residual_helmholtz_tv(t, v, &n).unwrap();
    approx_eq(
        after_b.f,
        before_b.f,
        0.0,
        1e-15,
        "Order (b) model preserved",
    );
}

#[test]
fn set_kij_rejects_overflowing_cross_dispersion_energy_on_polar_and_nonpolar() {
    // Polar CO2/H2O
    let mut polar = PcSaft::new_polar(&["CO2", "H2O"]).unwrap();
    let orig_polar_kij = polar.get_kij(0, 1);
    let orig_polar_eps = polar.multipol().unwrap().eps_divk_ij.clone();
    for bad_kij in [-f64::MAX, -1.0e160] {
        let err = polar.set_kij(0, 1, bad_kij).unwrap_err();
        assert!(matches!(err, Error::InvalidInput(_)));
        assert_eq!(polar.get_kij(0, 1), orig_polar_kij);
        assert_eq!(polar.get_kij(1, 0), orig_polar_kij);
        assert_eq!(polar.multipol().unwrap().eps_divk_ij, orig_polar_eps);
    }
    let f_polar = polar
        .residual_helmholtz_tv(340.0, 1.0e-4, &[0.5, 0.5])
        .unwrap();
    assert!(f_polar.f.is_finite());

    // Non-polar C1/C2
    let mut nonpolar = PcSaft::new(&["C1", "C2"]).unwrap();
    let orig_np_kij = nonpolar.get_kij(0, 1);
    for bad_kij in [-f64::MAX, -1.0e160] {
        let err = nonpolar.set_kij(0, 1, bad_kij).unwrap_err();
        assert!(matches!(err, Error::InvalidInput(_)));
        assert_eq!(nonpolar.get_kij(0, 1), orig_np_kij);
        assert_eq!(nonpolar.get_kij(1, 0), orig_np_kij);
    }
    let f_np = nonpolar
        .residual_helmholtz_tv(300.0, 1.0e-4, &[0.5, 0.5])
        .unwrap();
    assert!(f_np.f.is_finite());

    // Direct MultipolParams::new validation of sigma_ij, eps_divk_ij, m, mu, q.
    assert!(matches!(
        MultipolParams::new(&[1.0], &[0.0], &[1.0], &[3.0e-10], &[f64::INFINITY]),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        MultipolParams::new(&[1.0], &[0.0], &[1.0], &[3.0e-10], &[0.0]),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        MultipolParams::new(&[1.0], &[0.0], &[1.0], &[0.0], &[200.0]),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        MultipolParams::new(&[1.0], &[0.0], &[1.0], &[1.0e100], &[200.0]),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        MultipolParams::new(&[1.0], &[0.0], &[0.0], &[3.0e-10], &[200.0]),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        MultipolParams::new(&[f64::NAN], &[0.0], &[1.0], &[3.0e-10], &[200.0]),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn association_combining_rules_are_overflow_free_and_validate_finiteness() {
    // Direct AriComb and GeoComb overflow-free pure helper checks.
    let ari_max = AssocCombRule::AriComb.apply(f64::MAX, f64::MAX).unwrap();
    assert!(ari_max.is_finite());
    let geo_big = AssocCombRule::GeoComb.apply(1.0e200, 1.0e200).unwrap();
    approx_eq(geo_big, 1.0e200, 1e-14, 0.0, "GeoComb(1e200, 1e200)");

    let t = 340.0;
    let v = 1.0e-4;
    let n = [0.4, 0.6];

    // set_assoc_params rejects unphysical beta > MAX_ASSOC_BETA transactionally
    // across AriComb, GeoComb, and DefaultComb without corrupting the prior model.
    for rule in [
        AssocCombRule::AriComb,
        AssocCombRule::GeoComb,
        AssocCombRule::DefaultComb,
    ] {
        let mut m = PcSaft::new(&["NH3", "H2O"]).unwrap();
        m.set_assoc_comb_rules(0, 1, rule, rule).unwrap();
        let before = m.assoc_residual_helmholtz_tv(t, v, &n).unwrap();
        for bad_beta in [1.0e200, 1.0e308] {
            let err = m
                .set_assoc_params(0, 10_000.0, bad_beta, AssocScheme::Scheme4C)
                .unwrap_err();
            assert!(matches!(err, Error::InvalidInput(_)));
        }
        let after = m.assoc_residual_helmholtz_tv(t, v, &n).unwrap();
        approx_eq(
            after.f,
            before.f,
            0.0,
            1e-15,
            &format!("{rule:?} preserved after rejected beta"),
        );
    }

    // Supported upper bound beta = MAX_ASSOC_BETA (1e12) succeeds across all
    // combining rules and Delta rules and matches central finite differences.
    let max_beta = association::MAX_ASSOC_BETA;
    for comb in [
        AssocCombRule::AriComb,
        AssocCombRule::GeoComb,
        AssocCombRule::DefaultComb,
    ] {
        for delta_rule in [DeltaCombRule::Standard, DeltaCombRule::Elliot] {
            let mut m = PcSaft::new(&["NH3", "H2O"]).unwrap();
            m.set_assoc_comb_rules(0, 1, comb, comb).unwrap();
            m.set_delta_comb_rule(delta_rule).unwrap();
            m.set_assoc_params(0, 10_000.0, max_beta, AssocScheme::Scheme4C)
                .unwrap();
            m.set_assoc_params(1, 10_000.0, max_beta, AssocScheme::Scheme4C)
                .unwrap();
            assert!(m.association().beta_kl.iter().all(|v| v.is_finite()));

            let f0 = m.assoc_residual_helmholtz_tv(t, v, &n).unwrap();
            assert!(f0.f.is_finite() && f0.f_t.is_finite() && f0.f_v.is_finite());
            assert!(f0.f_tt.is_finite() && f0.f_tv.is_finite() && f0.f_vv.is_finite());
            assert!(f0.f_n.iter().all(|v| v.is_finite()));
            assert!(f0.f_tn.iter().all(|v| v.is_finite()));
            assert!(f0.f_vn.iter().all(|v| v.is_finite()));
            assert!(f0.f_nn.iter().all(|v| v.is_finite()));

            // Central FD check for Standard and Elliot at MAX_ASSOC_BETA.
            if comb == AssocCombRule::DefaultComb {
                let ht = 1.0e-4;
                let hv = 1.0e-10;
                let hn = 1.0e-6;
                let ftp = m.assoc_residual_helmholtz_tv(t + ht, v, &n).unwrap();
                let ftm = m.assoc_residual_helmholtz_tv(t - ht, v, &n).unwrap();
                let fvp = m.assoc_residual_helmholtz_tv(t, v + hv, &n).unwrap();
                let fvm = m.assoc_residual_helmholtz_tv(t, v - hv, &n).unwrap();
                let fnp = m
                    .assoc_residual_helmholtz_tv(t, v, &[n[0] + hn, n[1]])
                    .unwrap();
                let fnm = m
                    .assoc_residual_helmholtz_tv(t, v, &[n[0] - hn, n[1]])
                    .unwrap();

                approx_eq(
                    f0.f_t,
                    (ftp.f - ftm.f) / (2.0 * ht),
                    1e-6,
                    1e-8,
                    &format!("MAX_ASSOC_BETA {delta_rule:?} F_T"),
                );
                approx_eq(
                    f0.f_v,
                    (fvp.f - fvm.f) / (2.0 * hv),
                    1e-6,
                    1e-4,
                    &format!("MAX_ASSOC_BETA {delta_rule:?} F_V"),
                );
                approx_eq(
                    f0.f_n[0],
                    (fnp.f - fnm.f) / (2.0 * hn),
                    1e-6,
                    1e-7,
                    &format!("MAX_ASSOC_BETA {delta_rule:?} F_n[0]"),
                );
                approx_eq(
                    f0.f_tt,
                    (ftp.f_t - ftm.f_t) / (2.0 * ht),
                    1e-5,
                    1e-8,
                    &format!("MAX_ASSOC_BETA {delta_rule:?} F_TT"),
                );
                approx_eq(
                    f0.f_tv,
                    (ftp.f_v - ftm.f_v) / (2.0 * ht),
                    1e-5,
                    1e-5,
                    &format!("MAX_ASSOC_BETA {delta_rule:?} F_TV"),
                );
                approx_eq(
                    f0.f_vv,
                    (fvp.f_v - fvm.f_v) / (2.0 * hv),
                    1e-5,
                    1.0,
                    &format!("MAX_ASSOC_BETA {delta_rule:?} F_VV"),
                );
                approx_eq(
                    f0.f_nn[0],
                    (fnp.f_n[0] - fnm.f_n[0]) / (2.0 * hn),
                    1e-5,
                    1e-6,
                    &format!("MAX_ASSOC_BETA {delta_rule:?} F_nn[0,0]"),
                );
                approx_eq(
                    f0.f_nn[1],
                    (fnp.f_n[1] - fnm.f_n[1]) / (2.0 * hn),
                    1e-5,
                    1e-6,
                    &format!("MAX_ASSOC_BETA {delta_rule:?} F_nn[0,1]"),
                );
            }
        }
    }
}

#[test]
fn multipol_dq_structural_proportionality_rejects_kij_perturbation_and_allows_scaled_moments() {
    let t = 340.0;
    let v = 1.0e-4;

    // 1. Unmodified duplicate CO2 species match [1, 1e-30, 1e-30] exactly.
    let mut dup = PcSaft::new_polar(&["ACETONE", "CO2", "CO2"]).unwrap();
    dup.set_multipol_control(false, false, true).unwrap();
    let f_zero = dup
        .multipol_residual_helmholtz_tv(t, v, &[1.0, 0.0, 0.0])
        .unwrap();
    let f_tiny = dup
        .multipol_residual_helmholtz_tv(t, v, &[1.0, 1.0e-30, 1.0e-30])
        .unwrap();
    approx_eq(
        f_zero.f_nn[5],
        -1.417_186_792_393,
        1e-9,
        1e-9,
        "unmodified dup CO2 F_nn(1,2)",
    );
    approx_eq(
        f_zero.f_nn[5],
        f_tiny.f_nn[5],
        1e-12,
        1e-12,
        "unmodified dup CO2 boundary vs 1e-30",
    );

    // 2. Any nonzero kij perturbation breaks structural proportionality and must error at [1, 0, 0].
    for dk in [1.0e-12, 1.0e-15] {
        let mut pert = PcSaft::new_polar(&["ACETONE", "CO2", "CO2"]).unwrap();
        pert.set_multipol_control(false, false, true).unwrap();
        pert.set_kij(0, 2, dk).unwrap();
        let err = pert
            .multipol_residual_helmholtz_tv(t, v, &[1.0, 0.0, 0.0])
            .unwrap_err();
        assert!(
            matches!(err, Error::InvalidInput(_)),
            "expected InvalidInput for kij perturbation {dk:e}, got {err:?}"
        );
    }

    // 3. Structurally proportional species with scaled quadrupole moment (q_2 = 2 * q_1,
    //    identical m, sigma, eps_ij) remain smooth at [1, 0, 0] and match [1, 1e-30, 1e-30].
    let mut scaled = PcSaft::new_polar(&["ACETONE", "CO2", "CO2"]).unwrap();
    let q_co2 = scaled.components()[1].q;
    scaled.set_polar_moments(2, 0.0, 2.0 * q_co2).unwrap();
    scaled.set_multipol_control(false, false, true).unwrap();
    let fs_zero = scaled
        .multipol_residual_helmholtz_tv(t, v, &[1.0, 0.0, 0.0])
        .unwrap();
    let fs_tiny = scaled
        .multipol_residual_helmholtz_tv(t, v, &[1.0, 1.0e-30, 1.0e-30])
        .unwrap();
    approx_eq(
        fs_zero.f_n[2],
        4.0 * fs_zero.f_n[1],
        1e-12,
        1e-13,
        "scaled Q F_n scales as q^2",
    );
    approx_eq(
        fs_zero.f_nn[5],
        fs_tiny.f_nn[5],
        1e-11,
        1e-11,
        "scaled Q boundary vs 1e-30",
    );

    // 4. Both polar families absent in non-polar solvent: smooth zero mixed derivative
    //    even for distinct absent dipole species.
    let mut all_absent = PcSaft::new_polar(&["CO2", "ACETONE", "NH3", "C1"]).unwrap();
    all_absent.set_polar_moments(2, 1.47, 0.0).unwrap();
    all_absent.set_multipol_control(false, false, true).unwrap();
    let faa = all_absent
        .multipol_residual_helmholtz_tv(t, v, &[0.0, 0.0, 0.0, 1.0])
        .unwrap();
    approx_eq(
        faa.f_nn[6],
        0.0,
        0.0,
        1e-15,
        "both polar families absent F_nn(ACETONE, NH3) == 0",
    );
}

#[test]
fn multipol_eps_square_and_active_coeff_validation_and_ultra_trace_dd_qq_precision() {
    // Direct MultipolParams::new rejects eps^2 overflow (1e155) and active coeff overflow (mu=1e80),
    // while accepting eps=1e150 and evaluating cleanly.
    assert!(matches!(
        MultipolParams::new(&[1.0], &[0.0], &[1.0], &[3.0e-10], &[1.0e155]),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        MultipolParams::new(&[1.0e80], &[0.0], &[1.0], &[3.0e-10], &[200.0]),
        Err(Error::InvalidInput(_))
    ));
    let mp_big = MultipolParams::new(&[1.0], &[0.0], &[1.0], &[3.0e-10], &[1.0e150]).unwrap();
    let val_big = mp_big.fres_value(340.0, 1.0e-4, &[1.0]).unwrap();
    let all_big = mp_big.fres_all(340.0, 1.0e-4, &[1.0]).unwrap();
    assert!(val_big.is_finite() && all_big.f.is_finite() && all_big.f_nn[0].is_finite());

    let t = 340.0;
    let v = 1.0e-4;
    let eps1 = 1.0e-300;
    let eps2 = 2.0e-300;

    // DD ultra-trace at 1e-300: no subnormal precision loss in C1/ACETONE/MEOH gradient.
    let mut dd_model = PcSaft::new_polar(&["C1", "ACETONE", "MEOH"]).unwrap();
    dd_model.set_multipol_control(true, false, false).unwrap();
    let dd0 = dd_model
        .multipol_residual_helmholtz_tv(t, v, &[1.0, 0.0, 0.0])
        .unwrap();
    let dd_tr = dd_model
        .multipol_residual_helmholtz_tv(t, v, &[1.0, eps1, eps2])
        .unwrap();
    let dd_expected_fn1 = dd0.f_nn[4] * eps1 + dd0.f_nn[5] * eps2;
    assert!(dd_expected_fn1.abs() > 0.0);
    approx_eq(
        dd_tr.f_n[1],
        dd_expected_fn1,
        1e-12,
        0.0,
        "DD 1e-300 trace F_n(ACETONE) without subnormal loss",
    );

    // QQ ultra-trace at 1e-300: no subnormal precision loss in C1/CO2/N2 gradient.
    let mut qq_model = PcSaft::new_polar(&["C1", "CO2", "N2"]).unwrap();
    qq_model.set_polar_moments(2, 0.0, 1.5).unwrap();
    qq_model.set_multipol_control(false, true, false).unwrap();
    let qq0 = qq_model
        .multipol_residual_helmholtz_tv(t, v, &[1.0, 0.0, 0.0])
        .unwrap();
    let qq_tr = qq_model
        .multipol_residual_helmholtz_tv(t, v, &[1.0, eps1, eps2])
        .unwrap();
    let qq_expected_fn1 = qq0.f_nn[4] * eps1 + qq0.f_nn[5] * eps2;
    assert!(qq_expected_fn1.abs() > 0.0);
    approx_eq(
        qq_tr.f_n[1],
        qq_expected_fn1,
        1e-12,
        0.0,
        "QQ 1e-300 trace F_n(CO2) without subnormal loss",
    );
}

#[test]
fn multipol_dq_custom_diameters_have_a_smooth_mixed_boundary() {
    let s = 2.328_306_436_538_696_3e-10;
    let sigma = [s, s, 2.0 * s, s, s, 1.5 * s, 2.0 * s, 1.5 * s, 2.0 * s];
    let mut mp = MultipolParams::new(
        &[2.0, 0.0, 0.0],
        &[0.0, 3.0, 3.0],
        &[1.0; 3],
        &sigma,
        &[200.0; 9],
    )
    .unwrap();
    mp.enable_dd = false;
    mp.enable_qq = false;
    let (t, v) = (340.0, 1e-4);
    for (a, b) in [(0.0, 0.0), (1e-30, 1e-30), (2e-30, 1e-30), (1e-30, 2e-30)] {
        let full = mp.fres_all(t, v, &[1.0, a, b]).unwrap();
        approx_eq(
            full.f_nn[5],
            -14.054_142_104_610_15,
            0.0,
            2e-14,
            "DQ custom diameter boundary",
        );
    }

    // Away from the boundary the two-quadrupole terms have different diameter ratios.
    let n = [1.0, 0.08, 0.03];
    let h = 1e-6;
    let plus = mp.fres_fn(t, v, &[n[0], n[1], n[2] + h]).unwrap();
    let minus = mp.fres_fn(t, v, &[n[0], n[1], n[2] - h]).unwrap();
    approx_eq(
        mp.fres_all(t, v, &n).unwrap().f_nn[5],
        (plus[1] - minus[1]) / (2.0 * h),
        1e-8,
        0.0,
        "DQ interior mixed derivative",
    );
}

#[test]
fn multipol_dq_absent_dipoles_require_equal_diameter_products() {
    let s: f64 = 3.0 * 2.0_f64.powi(-33);
    let make_dq = |cross, pure| {
        let sigma = [s, s, cross, s, s, s, cross, s, pure];
        let mut mp = MultipolParams::new(
            &[0.0, 2.0, 2.0],
            &[3.0, 0.0, 0.0],
            &[1.0; 3],
            &sigma,
            &[200.0; 9],
        )
        .unwrap();
        mp.enable_dd = false;
        mp.enable_qq = false;
        mp
    };
    let (t, v) = (340.0, 1e-4);
    let mp = make_dq(s / 2.0, 2.0 * s);
    let boundary = mp.fres_all(t, v, &[1.0, 0.0, 0.0]).unwrap();
    for (a, b) in [(1e-30, 2e-30), (2e-30, 1e-30)] {
        let trace = mp.fres_all(t, v, &[1.0, a, b]).unwrap();
        approx_eq(
            trace.f_nn[5],
            boundary.f_nn[5],
            1e-13,
            0.0,
            "DQ dipole limit",
        );
    }

    // (s + ulp)*(s - ulp) rounds to s*s but differs by ulp squared.
    assert_eq!(s.next_up() * s.next_down(), s * s);
    let perturbed = make_dq(s.next_up(), s.next_down());
    assert!(matches!(
        perturbed.fres_all(t, v, &[1.0, 0.0, 0.0]),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn multipol_large_dipole_hessian_matches_gradient_differences() {
    let mp = MultipolParams::new(&[1e32], &[0.0], &[1.0], &[3e-10], &[200.0]).unwrap();
    let (t, v) = (340.0, 1e-4);
    let h = 1e-5;
    let plus = mp.fres_fn(t, v, &[1.0 + h]).unwrap();
    let minus = mp.fres_fn(t, v, &[1.0 - h]).unwrap();
    let fd = (plus[0] - minus[0]) / (2.0 * h);
    let full = mp.fres_all(t, v, &[1.0]).unwrap();
    approx_eq(
        full.f,
        -1.591_022_486_043_386e64,
        1e-14,
        0.0,
        "large DD value",
    );
    approx_eq(
        full.f_nn[0],
        2.863_395_002e63,
        1e-8,
        0.0,
        "large DD Hessian",
    );
    approx_eq(full.f_nn[0], fd, 1e-8, 0.0, "large DD gradient difference");
}

#[test]
fn multipol_dq_hessian_with_a_subnormal_pade_denominator() {
    let mut mp = MultipolParams::new(
        &[2e-12, 0.0, 0.0],
        &[0.0, 2e-12, 3e-12],
        &[1.0; 3],
        &[2.328_306_436_538_696_3e-10; 9],
        &[200.0; 9],
    )
    .unwrap();
    mp.enable_dd = false;
    mp.enable_qq = false;

    let (t, v, n) = (1e60, 1e145, [0.5, 0.25, 0.25]);
    let full = mp.fres_all(t, v, &n).unwrap();
    approx_eq(
        full.f,
        -1.478_840_498_596e-312,
        0.0,
        1e-320,
        "subnormal DQ value",
    );
    approx_eq(
        full.f_nn[1],
        -3.640_222_765_773e-312,
        0.0,
        1e-320,
        "subnormal DQ Hessian",
    );
}

#[test]
fn multipol_dq_hessian_correction_stays_finite_at_subnormal_traces() {
    let mut eps = [200.0; 9];
    eps[2] = 199.8;
    eps[6] = 199.8;
    let mut mp = MultipolParams::new(
        &[2.0, 0.0, 0.0],
        &[0.0, 3.0, 3.0],
        &[1.0; 3],
        &[2.328_306_436_538_696_3e-10; 9],
        &eps,
    )
    .unwrap();
    mp.enable_dd = false;
    mp.enable_qq = false;

    let (t, v, trace) = (340.0, 1e-4, 1e-312);
    let h = 0.05 * trace;
    let plus = mp.fres_fn(t, v, &[1.0, trace + h, trace]).unwrap();
    let minus = mp.fres_fn(t, v, &[1.0, trace - h, trace]).unwrap();
    let finite_difference = (plus[1] - minus[1]) / (2.0 * h);
    let hessian = mp.fres_all(t, v, &[1.0, trace, trace]).unwrap().f_nn[4];
    approx_eq(
        hessian,
        finite_difference,
        2e-3,
        0.0,
        "ultra-trace DQ Hessian",
    );

    let larger_trace = mp.fres_all(t, v, &[1.0, 1e-310, 1e-310]).unwrap().f_nn[4];
    approx_eq(
        hessian * trace,
        larger_trace * 1e-310,
        1e-8,
        0.0,
        "DQ trace-scaled Hessian",
    );
}
