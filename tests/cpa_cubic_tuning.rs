//! Comprehensive qualification for CPA second derivatives & mixing,
//! UMR/VTPR/HV-WS cubic mixing, apparent_composition, tuning, and complexmodelinit.
//!
//! Includes:
//! 1. Exact-source Fortran fixture / live parity (`exact_source_cpa_cubic_tuning`).
//! 2. Independent central finite-difference verification of all 10 TVN Helmholtz
//!    blocks and all TP `(T, P, n)` derivatives across CPA (`VdW`, `HV1`, `HV2`,
//!    `WS`, `HVWS`, original/simplified RDF, `Elliot`, `CPA-PR`), `UMR`, `VTPR`,
//!    `init_vtpr`, `init_umr`, `HVWongSandler`, and power-mean `b_exponent`.
//! 3. Thermodynamic identities (Maxwell symmetry, Gibbs-Duhem, Euler scaling
//!    invariance from `1e-140` to `1e12` mol, and exact infinite-dilution trace limits).
//! 4. Apparent-composition mathematical consistency (including multi-salt absent-ion
//!    stability fixing Fortran's `0 * ln(0) = NaN` defect) and tuning round trips.

use std::env;
use std::fs;
use thermopack_rs::wong_sandler::{Fraction, WongSandlerData, fidel_alpha};
use thermopack_rs::{
    AlphaCorrelation, ApparentComposition, Cpa, CpaCubicKind, CpaOptions, Cubic, CubicHelmholtz,
    CubicKind, CubicOptions, DeltaCombRule, HvPairParams, MixingRule, Phase, R, TranslationKind,
    VolumeShiftParams, WsPairParams, init_umr, init_vtpr,
};

fn assert_close(label: &str, actual: f64, expected: f64, atol: f64, rtol: f64) {
    assert!(
        actual.is_finite() && expected.is_finite(),
        "{label}: non-finite actual={actual:?} expected={expected:?}"
    );
    let diff = (actual - expected).abs();
    let tol = atol + rtol * expected.abs().max(actual.abs());
    assert!(
        diff <= tol,
        "{label}: actual={actual:.16e}, expected={expected:.16e}, diff={diff:.3e}, tol={tol:.3e}"
    );
}

fn parse_f64_list(s: &str) -> Vec<f64> {
    s.split(',')
        .map(|v| v.trim().parse::<f64>().expect("valid f64"))
        .collect()
}

fn parse_cubic_kind(s: &str) -> CubicKind {
    match s {
        "SRK" => CubicKind::SoaveRedlichKwong,
        "PR" => CubicKind::PengRobinson,
        "RK" => CubicKind::RedlichKwong,
        "VDW" => CubicKind::VanDerWaals,
        other => panic!("unknown cubic kind {other}"),
    }
}

fn parse_cpa_kind(s: &str) -> CpaCubicKind {
    match s {
        "SRK" => CpaCubicKind::Srk,
        "PR" => CpaCubicKind::Pr,
        other => panic!("unknown CPA kind {other}"),
    }
}

fn flatten_cubic_tvn(h: &CubicHelmholtz) -> Vec<f64> {
    let mut out = vec![h.f, h.f_t, h.f_v];
    out.extend_from_slice(&h.f_n);
    out.push(h.f_tt);
    out.push(h.f_tv);
    out.extend_from_slice(&h.f_tn);
    out.push(h.f_vv);
    out.extend_from_slice(&h.f_vn);
    out.extend_from_slice(&h.f_nn);
    out
}

fn flatten_cpa_tp(st: &thermopack_rs::CpaTpState) -> Vec<f64> {
    let mut out = vec![
        st.volume,
        st.compressibility_factor.value,
        st.compressibility_factor.dt,
        st.compressibility_factor.dp,
    ];
    out.extend_from_slice(&st.compressibility_factor.dn);
    out.push(st.residual_enthalpy.value);
    out.push(st.residual_enthalpy.dt);
    out.push(st.residual_enthalpy.dp);
    out.extend_from_slice(&st.residual_enthalpy.dn);
    out.push(st.enthalpy.value);
    out.push(st.enthalpy.dt);
    out.push(st.enthalpy.dp);
    out.extend_from_slice(&st.enthalpy.dn);
    out.push(st.residual_entropy.value);
    out.push(st.residual_entropy.dt);
    out.push(st.residual_entropy.dp);
    out.extend_from_slice(&st.residual_entropy.dn);
    out.push(st.entropy.value);
    out.push(st.entropy.dt);
    out.push(st.entropy.dp);
    out.extend_from_slice(&st.entropy.dn);
    out.extend_from_slice(&st.ln_fugacity.ln_phi);
    out.extend_from_slice(&st.ln_fugacity.dt);
    out.extend_from_slice(&st.ln_fugacity.dp);
    out.extend_from_slice(&st.ln_fugacity.dn);
    out.push(st.speed_of_sound);
    out
}

fn evaluate_cubic_tp_flat(model: &Cubic, t: f64, p: f64, moles: &[f64], phase: Phase) -> Vec<f64> {
    let nc = moles.len();
    let total: f64 = moles.iter().sum();
    let x: Vec<f64> = moles.iter().map(|&v| v / total).collect();
    let state = model.prepare(t, &x).unwrap();
    let v = state.volume(p, phase).unwrap();
    let h = state.residual_helmholtz_tv_unit(v).unwrap();
    let inv_n = 1.0 / total;
    let p_eval = R * t / v - R * t * h.f_v;
    let p_v = -R * t / (v * v) - R * t * h.f_vv;
    let p_t = p_eval / t - R * t * h.f_tv;
    let p_n: Vec<f64> = h.f_vn.iter().map(|&fvn| R * t / v - R * t * fvn).collect();

    let v_t = -p_t / p_v;
    let v_p = 1.0 / p_v;
    let v_n: Vec<f64> = p_n.iter().map(|&pni| -pni / p_v).collect();
    let z = p_eval * v / (R * t);
    let ln_z = z.ln();
    let z_t = (p_eval / (R * t)) * (v_t - v / t);
    let z_p = (v + p_eval * v_p) / (R * t);
    let z_dn: Vec<f64> = v_n
        .iter()
        .map(|&vni| (p_eval / (R * t)) * (vni - v) * inv_n)
        .collect();

    let h_res = -R * t * t * h.f_t + p_eval * v - R * t;
    let cv_res = -R * (2.0 * t * h.f_t + t * t * h.f_tt);
    let cp_res = cv_res - t * p_t * p_t / p_v - R;
    let h_res_p = v - t * v_t;
    let h_res_dn: Vec<f64> = (0..nc)
        .map(|i| -R * t * t * h.f_tn[i] + t * p_t * v_n[i] - R * t)
        .collect();

    let s_res = -R * (h.f + t * h.f_t - ln_z);
    let s_res_t = cp_res / t;
    let s_res_p = -v_t + R / p_eval;
    let s_res_dn: Vec<f64> = (0..nc)
        .map(|i| -R * (h.f_n[i] + t * h.f_tn[i] - ln_z) + p_t * v_n[i] - R)
        .collect();

    let mut h_id = 0.0;
    let mut cp_id = 0.0;
    let mut s_id = 0.0;
    let mut h_tot_dn = Vec::with_capacity(nc);
    let mut s_tot_dn = Vec::with_capacity(nc);
    for i in 0..nc {
        let id = model.ideal_gas(i).unwrap().properties(t, p_eval).unwrap();
        h_id += x[i] * id.enthalpy;
        cp_id += x[i] * id.cp;
        if x[i] > 0.0 {
            s_id += x[i] * (id.entropy - R * x[i].ln());
        }
        h_tot_dn.push(h_res_dn[i] + id.enthalpy);
        s_tot_dn.push(s_res_dn[i] + id.entropy - R * x[i].max(1.0e-100).ln());
    }

    let mut ln_phi = Vec::with_capacity(nc);
    let mut dlnphi_dt = Vec::with_capacity(nc);
    let mut dlnphi_dp = Vec::with_capacity(nc);
    let mut dlnphi_dn = vec![0.0; nc * nc];
    for i in 0..nc {
        ln_phi.push(h.f_n[i] - ln_z);
        dlnphi_dt.push(h.f_tn[i] + 1.0 / t - v_n[i] * p_t / (R * t));
        dlnphi_dp.push(v_n[i] / (R * t) - 1.0 / p_eval);
        for j in 0..nc {
            dlnphi_dn[i * nc + j] =
                (h.f_nn[i * nc + j] + 1.0 + p_n[i] * p_n[j] / (R * t * p_v)) * inv_n;
        }
    }
    let cv = cp_id - R + cv_res;
    let cp = cp_id + cp_res;
    let dp_drho = -p_v * v * v;
    let mw_kg: f64 = x
        .iter()
        .zip(model.components())
        .map(|(&xi, c)| xi * c.molar_mass * 1.0e-3)
        .sum();
    let sos = ((cp / cv) * dp_drho / mw_kg).sqrt();

    let mut out = vec![v, z, z_t, z_p];
    out.extend_from_slice(&z_dn);
    out.push(h_res);
    out.push(cp_res);
    out.push(h_res_p);
    out.extend_from_slice(&h_res_dn);
    out.push(h_id + h_res);
    out.push(cp);
    out.push(h_res_p);
    out.extend_from_slice(&h_tot_dn);
    out.push(s_res);
    out.push(s_res_t);
    out.push(s_res_p);
    out.extend_from_slice(&s_res_dn);
    out.push(s_id + s_res);
    out.push(cp / t);
    out.push(-v_t);
    out.extend_from_slice(&s_tot_dn);
    out.extend_from_slice(&ln_phi);
    out.extend_from_slice(&dlnphi_dt);
    out.extend_from_slice(&dlnphi_dp);
    out.extend_from_slice(&dlnphi_dn);
    out.push(sos);
    out
}

fn build_cpa_case(comps: &[&str], eos: &str, mix: &str, pref: &str, variant: &str) -> Cpa {
    let mut opts = CpaOptions {
        cubic_kind: parse_cpa_kind(eos),
        mixing_rule: MixingRule::from_name(mix).unwrap(),
        parameter_reference: Some(pref),
        simplified_rdf: if variant == "orig_rdf" {
            Some(false)
        } else {
            None
        },
        delta_comb_rule: if variant == "elliot" {
            DeltaCombRule::Elliot
        } else {
            DeltaCombRule::Standard
        },
        ..Default::default()
    };
    let ws_data;
    if variant == "ws_tuned" {
        let mut d = WongSandlerData::new(comps.len());
        d.set_alphaij(0, 1, 0.3).unwrap();
        d.set_alphaij(1, 0, 0.3).unwrap();
        d.set_kij(0, 1, Fraction::constant(0.08)).unwrap();
        d.set_kij(1, 0, Fraction::constant(0.08)).unwrap();
        d.set_tauij(0, 1, Fraction::constant(0.22)).unwrap();
        d.set_tauij(1, 0, Fraction::constant(-0.14)).unwrap();
        ws_data = d;
        opts.wong_sandler = Some(&ws_data);
    }
    let mut model = Cpa::with_cpa_options(comps, opts).unwrap();
    if variant == "hv1_tuned" {
        model
            .set_hv_param(
                0,
                1,
                HvPairParams {
                    alpha_ij: 0.25,
                    alpha_ji: 0.25,
                    a_ge_ij: 180.0,
                    a_ge_ji: -90.0,
                    b_ge_ij: 0.12,
                    b_ge_ji: -0.08,
                    c_ge_ij: 0.0,
                    c_ge_ji: 0.0,
                },
            )
            .unwrap();
    } else if variant == "hv2_tuned" {
        model
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
    }
    model
}

fn build_cubic_case(comps: &[&str], eos: &str, mix: &str, pref: &str, variant: &str) -> Cubic {
    if variant == "init_vtpr" {
        // Fortran cbCalcBmix ignores lowcase_bij; compare against VTPR with b_exponent = 1.0
        let mut m = init_vtpr(comps).unwrap();
        m.set_b_exponent(1.0).unwrap();
        return m;
    }
    let mut m = Cubic::with_options(
        parse_cubic_kind(eos),
        comps,
        CubicOptions {
            mixing_rule: MixingRule::from_name(mix).unwrap(),
            parameter_reference: pref,
            ..Default::default()
        },
    )
    .unwrap();
    if variant == "hv_ws_tuned" {
        if comps.len() == 2 {
            m.set_hv_param(
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
            m.set_kij(0, 1, 0.06).unwrap();
            m.set_lij(0, 1, 0.02).unwrap();
        } else {
            m.set_hv_param(
                0,
                1,
                HvPairParams {
                    alpha_ij: 0.20,
                    alpha_ji: 0.25,
                    a_ge_ij: 160.0,
                    a_ge_ji: -80.0,
                    b_ge_ij: 0.10,
                    b_ge_ji: -0.05,
                    c_ge_ij: 2.0e-4,
                    c_ge_ji: -1.5e-4,
                },
            )
            .unwrap();
            m.set_kij(0, 1, 0.04).unwrap();
            m.set_lij(0, 1, 0.015).unwrap();
        }
    }
    m
}

fn eval_apparent_case(case_name: &str) -> Vec<f64> {
    let (app, t, p, n, fe_n, fe_tn, fe_vn, fe_nn, lnf_r, dt_r, dp_r, dn_r) = match case_name {
        "h2o_nacl" => {
            let mut app = ApparentComposition::with_ncsym(2, 3, 1).unwrap();
            app.set_v_stoich(1, 1, 1.0).unwrap();
            app.set_v_stoich(1, 2, 1.0).unwrap();
            (
                app,
                310.0,
                1.0e5,
                vec![0.92, 0.08],
                vec![-0.4, -1.2, -1.5],
                vec![0.002, -0.005, 0.004],
                vec![-120.0, -310.0, -290.0],
                vec![0.8, 0.2, -0.1, 0.2, 1.5, 0.4, -0.1, 0.4, 1.7],
                vec![-0.25, -0.60, -0.70],
                vec![1.2e-3, -2.3e-3, -1.8e-3],
                vec![3.1e-8, -4.2e-8, 5.0e-8],
                vec![0.4, -0.1, 0.05, -0.1, 0.9, 0.25, 0.05, 0.25, 1.1],
            )
        }
        "h2o_nacl_cacl2" => {
            let mut app = ApparentComposition::with_ncsym(3, 4, 1).unwrap();
            app.set_v_stoich(1, 1, 1.0).unwrap();
            app.set_v_stoich(1, 2, 1.0).unwrap();
            app.set_v_stoich(2, 3, 1.0).unwrap();
            app.set_v_stoich(2, 2, 2.0).unwrap();
            (
                app,
                320.0,
                2.5e5,
                vec![0.85, 0.10, 0.05],
                vec![-0.3, -1.1, -1.4, -2.2],
                vec![0.001, -0.004, 0.003, -0.006],
                vec![-100.0, -280.0, -250.0, -410.0],
                vec![
                    0.7, 0.1, -0.05, 0.02, //
                    0.1, 1.4, 0.3, 0.15, //
                    -0.05, 0.3, 1.6, 0.25, //
                    0.02, 0.15, 0.25, 2.1,
                ],
                vec![-0.20, -0.55, -0.65, -1.10],
                vec![1.0e-3, -2.0e-3, -1.5e-3, -3.2e-3],
                vec![2.5e-8, -3.5e-8, 4.0e-8, -6.0e-8],
                vec![
                    0.35, -0.08, 0.04, -0.02, //
                    -0.08, 0.85, 0.20, 0.10, //
                    0.04, 0.20, 1.05, 0.18, //
                    -0.02, 0.10, 0.18, 1.40,
                ],
            )
        }
        "identity" => (
            ApparentComposition::new(2, 2).unwrap(),
            300.0,
            1.0e5,
            vec![0.4, 0.6],
            vec![-0.5, -0.8],
            vec![0.001, -0.002],
            vec![-90.0, -110.0],
            vec![0.6, 0.1, 0.1, 0.9],
            vec![-0.15, -0.22],
            vec![1.1e-3, -0.9e-3],
            vec![2.0e-8, -1.5e-8],
            vec![0.3, -0.05, -0.05, 0.5],
        ),
        other => panic!("unknown apparent case {other}"),
    };
    let ne = app.apparent_to_real_mole_numbers(&n).unwrap();
    let diffs = app
        .real_to_apparent_differentials(Some(&fe_n), Some(&fe_tn), Some(&fe_vn), Some(&fe_nn))
        .unwrap();
    let fug = app
        .tp_lnfug_apparent(&n, &ne, p, &lnf_r, Some(&dt_r), Some(&dp_r), Some(&dn_r))
        .unwrap();
    let sumn: f64 = n.iter().sum();
    let x: Vec<f64> = n.iter().map(|&v| v / sumn).collect();
    let mut mod_lnf = fug.ln_fugacity.clone();
    let sumne = app.mod_fugacity(t, p, &x, &mut mod_lnf).unwrap();

    let mut out = Vec::new();
    out.extend_from_slice(&ne);
    out.extend_from_slice(diffs.f_n.as_ref().unwrap());
    out.extend_from_slice(diffs.f_tn.as_ref().unwrap());
    out.extend_from_slice(diffs.f_vn.as_ref().unwrap());
    out.extend_from_slice(diffs.f_nn.as_ref().unwrap());
    out.extend_from_slice(&fug.ln_fugacity);
    out.extend_from_slice(fug.dt.as_ref().unwrap());
    out.extend_from_slice(fug.dp.as_ref().unwrap());
    out.extend_from_slice(fug.dn.as_ref().unwrap());
    out.extend_from_slice(&mod_lnf);
    out.push(sumne);
    out
}

fn eval_tuning_case(cid: &str) -> Vec<f64> {
    let parts: Vec<&str> = cid.split('|').collect();
    let kind = parts[1];
    let comps: Vec<&str> = parts[2].split(',').collect();
    let eos = parse_cubic_kind(parts[3]);
    let mut out = Vec::new();
    match kind {
        "hv2" => {
            let mut m = Cubic::with_options(
                eos,
                &comps,
                CubicOptions {
                    mixing_rule: MixingRule::HuronVidal2,
                    ..Default::default()
                },
            )
            .unwrap();
            let p = HvPairParams {
                alpha_ij: 0.22,
                alpha_ji: 0.27,
                a_ge_ij: 310.0,
                a_ge_ji: -140.0,
                b_ge_ij: 0.20,
                b_ge_ji: -0.15,
                c_ge_ij: 5.0e-4,
                c_ge_ji: -4.0e-4,
            };
            m.set_hv_param(0, 1, p).unwrap();
            let g = m.hv_param(0, 1).unwrap();
            out.extend([
                g.alpha_ij, g.alpha_ji, g.a_ge_ij, g.a_ge_ji, g.b_ge_ij, g.b_ge_ji, g.c_ge_ij,
                g.c_ge_ji,
            ]);
            out.extend(flatten_cubic_tvn(
                &m.residual_helmholtz_tv(320.0, 2.2e-4, &[0.4, 0.6]).unwrap(),
            ));
        }
        "ws" => {
            let mut m = Cubic::with_options(
                eos,
                &comps,
                CubicOptions {
                    mixing_rule: MixingRule::WongSandler,
                    ..Default::default()
                },
            )
            .unwrap();
            let p = WsPairParams {
                alpha_ij: 0.30,
                alpha_ji: 0.35,
                k_ij: 0.11,
                k_ji: 0.11,
                tau_ij: 0.25,
                tau_ji: -0.18,
            };
            m.set_ws_param(0, 1, p).unwrap();
            let g = m.ws_param(0, 1).unwrap();
            out.extend([g.alpha_ij, g.alpha_ji, g.k_ij, g.k_ji, g.tau_ij, g.tau_ji]);
            out.extend(flatten_cubic_tvn(
                &m.residual_helmholtz_tv(320.0, 2.2e-4, &[0.4, 0.6]).unwrap(),
            ));
        }
        "twu_mc_beta" => {
            let mut m = Cubic::new(eos, &comps).unwrap();
            m.set_twu_param(0, 0.21, 0.88, 1.95).unwrap();
            m.set_mc_param(1, 0.72, 0.18, -0.09).unwrap();
            m.set_beta_corr(1, "QUANTUM", &[3.0696, 12.682]).unwrap();
            out.extend(m.twu_param(0).unwrap());
            out.extend(m.mc_param(1).unwrap());
            out.extend(flatten_cubic_tvn(
                &m.residual_helmholtz_tv(180.0, 1.8e-4, &[0.5, 0.5]).unwrap(),
            ));
        }
        "alpha_name" => {
            let mut m = Cubic::new(eos, &comps).unwrap();
            m.set_alpha_corr(0, "SW", &[0.42]).unwrap();
            m.set_alpha_corr(1, "PR78", &[0.61]).unwrap();
            out.extend([0.42, 0.61]);
            out.extend(flatten_cubic_tvn(
                &m.residual_helmholtz_tv(280.0, 2.4e-4, &[0.45, 0.55])
                    .unwrap(),
            ));
        }
        "kij_lij" => {
            let mut m = Cubic::new(eos, &comps).unwrap();
            m.set_kij(0, 1, 0.12).unwrap();
            m.set_kij_entry(1, 0, 0.05).unwrap();
            m.set_lij(0, 1, 0.03).unwrap();
            out.extend([
                m.kij(0, 1).unwrap(),
                m.kij(1, 0).unwrap(),
                m.lij(0, 1).unwrap(),
                m.lij(1, 0).unwrap(),
            ]);
            out.extend(flatten_cubic_tvn(
                &m.residual_helmholtz_tv(300.0, 2.0e-4, &[0.4, 0.6]).unwrap(),
            ));
        }
        "vshift_quad" => {
            let mut m = Cubic::new(eos, &comps).unwrap();
            m.set_volume_shift_parameters(
                0,
                VolumeShiftParams {
                    coefficients: [1.2e-6, -2.5e-9, 1.5e-12, 0.0, 0.0, 0.0],
                    ci_type: TranslationKind::Quadratic,
                },
            )
            .unwrap();
            m.set_volume_translation_enabled(true);
            let g = m.volume_shift_parameters(0).unwrap();
            out.extend(g.coefficients);
            out.push(g.ci_type as u8 as f64);
            out.extend(flatten_cubic_tvn(
                &m.residual_helmholtz_tv(310.0, 2.1e-4, &[0.4, 0.6]).unwrap(),
            ));
        }
        "unifac_get" => {
            let m = Cubic::with_options(
                eos,
                &comps,
                CubicOptions {
                    mixing_rule: MixingRule::Unifac,
                    ..Default::default()
                },
            )
            .unwrap();
            let (a, b, c) = m.unifac_group_interaction(1, 57).unwrap();
            out.extend([a, b, c]);
            out.extend(flatten_cubic_tvn(
                &m.residual_helmholtz_tv(310.0, 2.0e-4, &[0.4, 0.6]).unwrap(),
            ));
        }
        other => panic!("unknown tuning case {other}"),
    }
    out
}

#[test]
fn exact_source_cpa_cubic_tuning() {
    let content = match env::var("CPA_CUBIC_TUNING_FIXTURE") {
        Ok(path) => fs::read_to_string(path).expect("read live fixture"),
        Err(_) => include_str!("fixtures/cpa_cubic_tuning.csv").to_string(),
    };
    let mut rows = content.lines();
    let header = rows.next().expect("csv header");
    assert_eq!(header, "case_id;category;tag;values");

    let mut count = 0usize;
    for line in rows {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split(';').collect();
        assert_eq!(cols.len(), 4, "malformed row: {line}");
        let (cid, category, _tag, val_str) = (cols[0], cols[1], cols[2], cols[3]);
        let expected = parse_f64_list(val_str);
        let parts: Vec<&str> = cid.split('|').collect();

        let (actual, atol, rtol) = match category {
            "CPA_TVN" => {
                let comps: Vec<&str> = parts[1].split(',').collect();
                let t: f64 = parts[6].parse().unwrap();
                let v: f64 = parts[7].parse().unwrap();
                let moles = parse_f64_list(parts[8]);
                let model = build_cpa_case(&comps, parts[2], parts[3], parts[4], parts[5]);
                let h = model.residual_helmholtz_tv(t, v, &moles).unwrap();
                let mut flat = vec![h.f, h.f_t, h.f_v];
                flat.extend_from_slice(&h.f_n);
                flat.push(h.f_tt);
                flat.push(h.f_tv);
                flat.extend_from_slice(&h.f_tn);
                flat.push(h.f_vv);
                flat.extend_from_slice(&h.f_vn);
                flat.extend_from_slice(&h.f_nn);
                flat.extend_from_slice(&h.x_k);
                (flat, 1.0e-11, 1.0e-10)
            }
            "CPA_TP" => {
                let comps: Vec<&str> = parts[1].split(',').collect();
                let phase = if parts[6] == "1" {
                    Phase::Liquid
                } else {
                    Phase::Vapor
                };
                let t: f64 = parts[7].parse().unwrap();
                let p: f64 = parts[8].parse().unwrap();
                let moles = parse_f64_list(parts[9]);
                let model = build_cpa_case(&comps, parts[2], parts[3], parts[4], parts[5]);
                let st = model.state_tp(t, p, &moles, phase).unwrap();
                (flatten_cpa_tp(&st), 1.0e-9, 1.0e-8)
            }
            "CUBIC_TVN" => {
                let comps: Vec<&str> = parts[1].split(',').collect();
                let t: f64 = parts[6].parse().unwrap();
                let v: f64 = parts[7].parse().unwrap();
                let moles = parse_f64_list(parts[8]);
                let model = build_cubic_case(&comps, parts[2], parts[3], parts[4], parts[5]);
                let h = model.residual_helmholtz_tv(t, v, &moles).unwrap();
                (flatten_cubic_tvn(&h), 1.0e-11, 1.0e-10)
            }
            "CUBIC_TP" => {
                let comps: Vec<&str> = parts[1].split(',').collect();
                let phase = if parts[6] == "1" {
                    Phase::Liquid
                } else {
                    Phase::Vapor
                };
                let t: f64 = parts[7].parse().unwrap();
                let p: f64 = parts[8].parse().unwrap();
                let moles = parse_f64_list(parts[9]);
                let model = build_cubic_case(&comps, parts[2], parts[3], parts[4], parts[5]);
                (
                    evaluate_cubic_tp_flat(&model, t, p, &moles, phase),
                    1.0e-9,
                    1.0e-8,
                )
            }
            "FIDEL" => {
                let tau12: f64 = parts[1].parse().unwrap();
                let tau21: f64 = parts[2].parse().unwrap();
                (vec![fidel_alpha(tau12, tau21).unwrap()], 1.0e-14, 1.0e-12)
            }
            "BMIX" => {
                let comps: Vec<&str> = parts[1].split(',').collect();
                let b_exp: f64 = parts[3].parse().unwrap();
                let moles = parse_f64_list(parts[4]);
                let model = if parts[2] == "init_vtpr" {
                    init_vtpr(&comps).unwrap()
                } else {
                    let mut m = Cubic::new(CubicKind::PengRobinson, &comps).unwrap();
                    m.set_b_exponent(b_exp).unwrap();
                    m
                };
                let mut out = vec![model.calc_bmix_cubic(&moles).unwrap()];
                for i in 0..comps.len() {
                    for j in 0..comps.len() {
                        out.push(model.lowcase_bij(i, j).unwrap());
                    }
                }
                (out, 1.0e-15, 1.0e-12)
            }
            "APPARENT" => (eval_apparent_case(parts[1]), 1.0e-13, 1.0e-12),
            "TUNING" => (eval_tuning_case(cid), 1.0e-11, 1.0e-10),
            other => panic!("unknown category {other}"),
        };

        assert_eq!(
            actual.len(),
            expected.len(),
            "{cid}: length mismatch actual={} expected={}",
            actual.len(),
            expected.len()
        );
        for (idx, (&a, &e)) in actual.iter().zip(expected.iter()).enumerate() {
            assert_close(&format!("{cid}[{idx}]"), a, e, atol, rtol);
        }
        count += 1;
    }
    assert_eq!(count, 57, "expected all 57 pinned parity cases");
}

#[test]
fn independent_fd_all_cpa_mixing_and_derivatives() {
    let configs: Vec<(&str, Cpa, f64, f64, Vec<f64>)> = vec![
        (
            "CPA-SRK VdW H2O+MEOH+C1",
            Cpa::new(&["H2O", "MEOH", "C1"]).unwrap(),
            370.0,
            1.8e-4,
            vec![0.7, 0.5, 0.3],
        ),
        (
            "CPA-SRK orig_rdf H2O+ETOH",
            Cpa::with_cpa_options(
                &["H2O", "ETOH"],
                CpaOptions {
                    simplified_rdf: Some(false),
                    ..Default::default()
                },
            )
            .unwrap(),
            350.0,
            5.5e-5,
            vec![0.6, 0.8],
        ),
        (
            "CPA-SRK Elliot H2O+MEOH",
            Cpa::with_cpa_options(
                &["H2O", "MEOH"],
                CpaOptions {
                    delta_comb_rule: DeltaCombRule::Elliot,
                    ..Default::default()
                },
            )
            .unwrap(),
            345.0,
            4.5e-5,
            vec![0.55, 0.65],
        ),
        (
            "CPA-SRK HV2 H2O+MEOH",
            build_cpa_case(&["H2O", "MEOH"], "SRK", "HV2", "DEFAULT", "hv2_tuned"),
            355.0,
            4.6e-5,
            vec![0.5, 0.7],
        ),
        (
            "CPA-SRK WS H2O+CO2",
            build_cpa_case(&["H2O", "CO2"], "SRK", "WS", "DEFAULT", "ws_tuned"),
            385.0,
            1.7e-4,
            vec![0.6, 0.4],
        ),
        (
            "CPA-PR HVWS H2O+CO2",
            {
                let mut m = Cpa::with_cpa_options(
                    &["H2O", "CO2"],
                    CpaOptions {
                        cubic_kind: CpaCubicKind::Pr,
                        mixing_rule: MixingRule::HvWongSandler,
                        ..Default::default()
                    },
                )
                .unwrap();
                m.set_hv_param(
                    0,
                    1,
                    HvPairParams {
                        alpha_ij: 0.25,
                        alpha_ji: 0.30,
                        a_ge_ij: 210.0,
                        a_ge_ji: -95.0,
                        b_ge_ij: 0.10,
                        b_ge_ji: -0.05,
                        c_ge_ij: 2.0e-4,
                        c_ge_ji: -1.0e-4,
                    },
                )
                .unwrap();
                m.try_set_kij(0, 1, 0.05).unwrap();
                m.set_lij(0, 1, 0.02).unwrap();
                m
            },
            390.0,
            1.9e-4,
            vec![0.55, 0.45],
        ),
    ];

    for (label, model, t, v, n) in configs {
        let nc = n.len();
        let h = model.residual_helmholtz_tv(t, v, &n).unwrap();
        let dt = t * 1.0e-5;
        let dv = v * 1.0e-5;
        let hp_t = model.residual_helmholtz_tv(t + dt, v, &n).unwrap();
        let hm_t = model.residual_helmholtz_tv(t - dt, v, &n).unwrap();
        let hp_v = model.residual_helmholtz_tv(t, v + dv, &n).unwrap();
        let hm_v = model.residual_helmholtz_tv(t, v - dv, &n).unwrap();

        assert_close(
            &format!("{label} F_T"),
            h.f_t,
            (hp_t.f - hm_t.f) / (2.0 * dt),
            1e-8,
            1e-6,
        );
        assert_close(
            &format!("{label} F_TT"),
            h.f_tt,
            (hp_t.f_t - hm_t.f_t) / (2.0 * dt),
            1e-8,
            1e-6,
        );
        assert_close(
            &format!("{label} F_V"),
            h.f_v,
            (hp_v.f - hm_v.f) / (2.0 * dv),
            1e-6,
            1e-6,
        );
        assert_close(
            &format!("{label} F_VV"),
            h.f_vv,
            (hp_v.f_v - hm_v.f_v) / (2.0 * dv),
            1e-4,
            1e-6,
        );
        assert_close(
            &format!("{label} F_TV"),
            h.f_tv,
            (hp_t.f_v - hm_t.f_v) / (2.0 * dt),
            1e-6,
            1e-6,
        );

        for i in 0..nc {
            let dn = n[i] * 1.0e-5;
            let mut np = n.clone();
            let mut nm = n.clone();
            np[i] += dn;
            nm[i] -= dn;
            let hp = model.residual_helmholtz_tv(t, v, &np).unwrap();
            let hm = model.residual_helmholtz_tv(t, v, &nm).unwrap();
            assert_close(
                &format!("{label} F_n[{i}]"),
                h.f_n[i],
                (hp.f - hm.f) / (2.0 * dn),
                1e-8,
                1e-6,
            );
            assert_close(
                &format!("{label} F_Tn[{i}]"),
                h.f_tn[i],
                (hp.f_t - hm.f_t) / (2.0 * dn),
                1e-8,
                1e-6,
            );
            assert_close(
                &format!("{label} F_Vn[{i}]"),
                h.f_vn[i],
                (hp.f_v - hm.f_v) / (2.0 * dn),
                1e-6,
                1e-6,
            );
            for j in 0..nc {
                assert_close(
                    &format!("{label} F_nn[{i},{j}]"),
                    h.f_nn[i * nc + j],
                    (hp.f_n[j] - hm.f_n[j]) / (2.0 * dn),
                    1e-7,
                    1e-6,
                );
                assert_close(
                    &format!("{label} symmetry F_nn[{i},{j}]"),
                    h.f_nn[i * nc + j],
                    h.f_nn[j * nc + i],
                    1e-12,
                    1e-12,
                );
            }
        }

        // Verify Euler homogeneity: F = sum_i n_i F_n[i] + V F_V
        let euler: f64 = n
            .iter()
            .zip(&h.f_n)
            .map(|(&ni, &fni)| ni * fni)
            .sum::<f64>()
            + v * h.f_v;
        assert_close(&format!("{label} Euler"), h.f, euler, 1e-11, 1e-11);
    }
}

#[test]
fn cpa_scaling_invariance_and_zero_mole_infinite_dilution() {
    let model = Cpa::new(&["H2O", "MEOH", "C1"]).unwrap();
    let t = 360.0;
    let v0 = 1.5e-4;
    let n0 = [0.5, 0.3, 0.2];
    let base = model.residual_helmholtz_tv(t, v0, &n0).unwrap();

    for scale in [1.0e-140, 1.0e-40, 1.0e-6, 1.0e6, 1.0e12] {
        let ns: Vec<f64> = n0.iter().map(|&x| x * scale).collect();
        let hs = model.residual_helmholtz_tv(t, v0 * scale, &ns).unwrap();
        assert_close(
            &format!("scale={scale:e} F/s"),
            hs.f / scale,
            base.f,
            1e-12,
            1e-11,
        );
        assert_close(
            &format!("scale={scale:e} F_V"),
            hs.f_v,
            base.f_v,
            1e-10,
            1e-11,
        );
        assert_close(
            &format!("scale={scale:e} F_n[0]"),
            hs.f_n[0],
            base.f_n[0],
            1e-12,
            1e-11,
        );
        assert_close(
            &format!("scale={scale:e} X_k[0]"),
            hs.x_k[0],
            base.x_k[0],
            1e-13,
            1e-12,
        );
    }

    // Zero-mole infinite-dilution cross-association check: n_MEOH = 0
    // Analytic F_n[MEOH] must match forward finite difference (F(1, eps) - F(1, 0)) / eps.
    let bin = Cpa::new(&["H2O", "MEOH"]).unwrap();
    let v = 2.0e-5;
    let h0 = bin.residual_helmholtz_tv(340.0, v, &[1.0, 0.0]).unwrap();
    let eps = 1.0e-7;
    let hp = bin.residual_helmholtz_tv(340.0, v, &[1.0, eps]).unwrap();
    let fd = (hp.f - h0.f) / eps;
    // MEOH sites at index 4, 5 must be physically < 1.0 (bonded to water!), not Fortran's buggy 1.0.
    assert!(
        h0.x_k[4] < 0.5 && h0.x_k[5] < 0.5,
        "infinite dilution X_k={:?}",
        h0.x_k
    );
    assert_close("zero-mole CPA F_n[MEOH]", h0.f_n[1], fd, 1e-5, 1e-5);
}

#[test]
fn complex_models_vtpr_umr_and_b_exponent_fd_consistency() {
    let vtpr = init_vtpr(&["C1", "C3", "NC4"]).unwrap();
    assert_close(
        "VTPR b_exponent",
        vtpr.b_exponent(),
        4.0 / 3.0,
        1e-15,
        1e-15,
    );
    assert!(vtpr.volume_translation_enabled());

    let umr = init_umr(&["C1", "C3", "H2O"]).unwrap();
    assert_close("UMR b_exponent", umr.b_exponent(), 2.0, 1e-15, 1e-15);
    assert_eq!(
        umr.alpha(0).unwrap().correlation,
        AlphaCorrelation::ClassicFit
    );
    assert!(umr.volume_translation_enabled());

    let umr_mc =
        thermopack_rs::init_umr_with_refs(&["C1", "C3", "H2O"], "DEFAULT", "Chapoy2005").unwrap();
    assert_eq!(
        umr_mc.alpha(2).unwrap().correlation,
        AlphaCorrelation::MathiasCopeman
    );

    // Verify full 10-block Helmholtz FD for both physical VTPR (s=4/3) and UMR (s=2)
    // AND quantum-beta + b_exponent=2.0 combination!
    let mut q_pow = Cubic::quantum(&["H2", "HE"]).unwrap();
    q_pow.set_b_exponent(2.0).unwrap();
    q_pow.set_lij(0, 1, 0.04).unwrap();

    for (label, model, t, v, n) in [
        ("VTPR s=4/3", &vtpr, 330.0, 3.0e-4, vec![0.4, 0.35, 0.25]),
        ("UMR s=2.0", &umr, 380.0, 4.0e-4, vec![0.4, 0.35, 0.25]),
        ("Quantum b_exp=2", &q_pow, 60.0, 2.0e-4, vec![0.6, 0.4]),
    ] {
        let nc = n.len();
        let h = model.residual_helmholtz_tv(t, v, &n).unwrap();
        let dt = t * 1.0e-5;
        let dv = v * 1.0e-5;
        let hp_t = model.residual_helmholtz_tv(t + dt, v, &n).unwrap();
        let hm_t = model.residual_helmholtz_tv(t - dt, v, &n).unwrap();
        let hp_v = model.residual_helmholtz_tv(t, v + dv, &n).unwrap();
        let hm_v = model.residual_helmholtz_tv(t, v - dv, &n).unwrap();
        assert_close(
            &format!("{label} F_T"),
            h.f_t,
            (hp_t.f - hm_t.f) / (2.0 * dt),
            1e-8,
            1e-6,
        );
        assert_close(
            &format!("{label} F_TT"),
            h.f_tt,
            (hp_t.f_t - hm_t.f_t) / (2.0 * dt),
            1e-8,
            1e-6,
        );
        assert_close(
            &format!("{label} F_V"),
            h.f_v,
            (hp_v.f - hm_v.f) / (2.0 * dv),
            1e-6,
            1e-6,
        );
        assert_close(
            &format!("{label} F_VV"),
            h.f_vv,
            (hp_v.f_v - hm_v.f_v) / (2.0 * dv),
            1e-4,
            1e-6,
        );
        assert_close(
            &format!("{label} F_TV"),
            h.f_tv,
            (hp_t.f_v - hm_t.f_v) / (2.0 * dt),
            1e-6,
            1e-6,
        );
        for i in 0..nc {
            let dn = n[i] * 1.0e-5;
            let mut np = n.clone();
            let mut nm = n.clone();
            np[i] += dn;
            nm[i] -= dn;
            let hp = model.residual_helmholtz_tv(t, v, &np).unwrap();
            let hm = model.residual_helmholtz_tv(t, v, &nm).unwrap();
            assert_close(
                &format!("{label} F_n[{i}]"),
                h.f_n[i],
                (hp.f - hm.f) / (2.0 * dn),
                1e-8,
                1e-6,
            );
            assert_close(
                &format!("{label} F_Tn[{i}]"),
                h.f_tn[i],
                (hp.f_t - hm.f_t) / (2.0 * dn),
                1e-8,
                1e-6,
            );
            assert_close(
                &format!("{label} F_Vn[{i}]"),
                h.f_vn[i],
                (hp.f_v - hm.f_v) / (2.0 * dn),
                1e-6,
                1e-6,
            );
            for j in 0..nc {
                assert_close(
                    &format!("{label} F_nn[{i},{j}]"),
                    h.f_nn[i * nc + j],
                    (hp.f_n[j] - hm.f_n[j]) / (2.0 * dn),
                    1e-7,
                    1e-6,
                );
            }
        }
        // Check fugacity_derivatives_into consistency with TP state at the same root
        let prep = model.prepare(t, &n).unwrap();
        let p = 2.0e6;
        let v_root = prep.volume(p, Phase::Vapor).unwrap();
        let (mut d_t, mut d_p, mut d_n) = (vec![0.0; nc], vec![0.0; nc], vec![0.0; nc * nc]);
        prep.fugacity_derivatives_into(v_root, &mut d_t, &mut d_p, &mut d_n)
            .unwrap();
        let tp = evaluate_cubic_tp_flat(model, t, p, &n, Phase::Vapor);
        let off_dt = 16 + 6 * nc;
        for i in 0..nc {
            assert_close(
                &format!("{label} dlnphi_dt[{i}]"),
                d_t[i],
                tp[off_dt + i],
                1e-9,
                1e-8,
            );
            assert_close(
                &format!("{label} dlnphi_dp[{i}]"),
                d_p[i],
                tp[off_dt + nc + i],
                1e-11,
                1e-8,
            );
            for j in 0..nc {
                assert_close(
                    &format!("{label} dlnphi_dn[{i},{j}]"),
                    d_n[i * nc + j],
                    tp[off_dt + 2 * nc + i * nc + j],
                    1e-9,
                    1e-8,
                );
            }
        }
    }
}

#[test]
fn apparent_composition_fixes_absent_salt_nan_and_matches_fd() {
    // Multi-salt H2O + NaCl + CaCl2 with n_CaCl2 = 0 -> Fortran produced NaN for NaCl!
    let mut app = ApparentComposition::with_ncsym(3, 5, 1).unwrap();
    app.set_v_stoich(1, 1, 1.0).unwrap();
    app.set_v_stoich(1, 2, 1.0).unwrap();
    app.set_v_stoich(2, 3, 1.0).unwrap();
    app.set_v_stoich(2, 2, 2.0).unwrap();

    let n = [0.9, 0.1, 0.0];
    let ne = app.apparent_to_real_mole_numbers(&n).unwrap();
    assert_eq!(ne, vec![0.9, 0.1, 0.1, 0.0, 0.0]);
    let lnf_r = [-0.2, -0.5, -0.6, -0.8, 0.0];
    let res = app
        .tp_lnfug_apparent(&n, &ne, 1.0e5, &lnf_r, None, None, None)
        .unwrap();
    assert!(
        res.ln_fugacity.iter().all(|v| v.is_finite()),
        "{:?}",
        res.ln_fugacity
    );
    assert_eq!(res.ln_fugacity[2], 0.0);

    // Independent FD check of apparent dlnfug/dn on H2O + NaCl + CaCl2 using a real Cubic EOS
    let eos = Cubic::new(
        CubicKind::SoaveRedlichKwong,
        &["H2O", "C1", "C2", "C3", "NC4"],
    )
    .unwrap();
    let t = 350.0;
    let p = 2.0e6;
    let n_base = [0.80, 0.12, 0.08];
    let eval = |nv: &[f64]| {
        let nev = app.apparent_to_real_mole_numbers(nv).unwrap();
        let sumne: f64 = nev.iter().sum();
        let xe: Vec<f64> = nev.iter().map(|&v| v / sumne).collect();
        let st = eos.prepare(t, &xe).unwrap();
        let v = st.volume(p, Phase::Vapor).unwrap();
        let mut lnf = vec![0.0; 5];
        let mut dt = vec![0.0; 5];
        let mut dp = vec![0.0; 5];
        let mut dn = vec![0.0; 25];
        st.ln_fugacity_tp_into(p, v, &mut lnf).unwrap();
        st.fugacity_derivatives_into(v, &mut dt, &mut dp, &mut dn)
            .unwrap();
        for d in &mut dn {
            *d /= sumne;
        }
        app.tp_lnfug_apparent(nv, &nev, p, &lnf, Some(&dt), Some(&dp), Some(&dn))
            .unwrap()
    };
    let base = eval(&n_base);
    let dn_anal = base.dn.as_ref().unwrap();
    for j in 0..3 {
        let eps = 1.0e-6;
        let mut np = n_base;
        let mut nm = n_base;
        np[j] += eps;
        nm[j] -= eps;
        let rp = eval(&np);
        let rm = eval(&nm);
        for i in 0..3 {
            let fd = (rp.ln_fugacity[i] - rm.ln_fugacity[i]) / (2.0 * eps);
            assert_close(
                &format!("apparent dlnfug[{i},{j}]"),
                dn_anal[i * 3 + j],
                fd,
                1e-6,
                1e-6,
            );
        }
    }
}

#[test]
fn unifac_group_interaction_setter_and_critical_redefinition() {
    let mut m = Cubic::with_options(
        CubicKind::PengRobinson,
        &["C1", "C3"],
        CubicOptions {
            mixing_rule: MixingRule::Unifac,
            ..Default::default()
        },
    )
    .unwrap();
    let before = m.residual_helmholtz_tv(310.0, 2.0e-4, &[0.4, 0.6]).unwrap();
    let (a0, b0, c0) = m.unifac_group_interaction(1, 57).unwrap();
    m.set_unifac_group_interaction(1, 57, a0 + 25.0, b0 + 0.05, c0)
        .unwrap();
    let (a1, b1, c1) = m.unifac_group_interaction(1, 57).unwrap();
    assert_close("unifac a1", a1, a0 + 25.0, 1e-14, 1e-14);
    assert_close("unifac b1", b1, b0 + 0.05, 1e-14, 1e-14);
    assert_close("unifac c1", c1, c0, 1e-14, 1e-14);
    let after = m.residual_helmholtz_tv(310.0, 2.0e-4, &[0.4, 0.6]).unwrap();
    assert!((after.f - before.f).abs() > 1e-6);

    // Critical parameter redefinition
    let mut c = Cubic::new(CubicKind::PengRobinson, &["C1", "C3"]).unwrap();
    let tc_new = 195.0;
    let pc_new = 4.7e6;
    let acf_new = 0.02;
    c.redefine_critical_parameters(0, tc_new, pc_new, acf_new)
        .unwrap();
    assert_eq!(c.components()[0].critical_temperature, tc_new);
    assert_eq!(c.components()[0].critical_pressure, pc_new);
    assert_eq!(c.components()[0].acentric_factor, acf_new);
}
