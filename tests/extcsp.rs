use thermopack_rs::{
    AlphaCorrelation, Cubic, CubicKind, Error, ExtCsp, ExtCspOptions, ExtCspReferenceKind, Mbwr,
    MbwrKind, MixingRule, Multiparameter, Phase, R, Result,
};

fn within_physical_tol(
    expected: f64,
    actual: f64,
    phys_scale: f64,
    rtol: f64,
    atol_scale: f64,
) -> bool {
    if !expected.is_finite() || !actual.is_finite() {
        return false;
    }
    let tol = atol_scale * phys_scale + rtol * expected.abs().max(actual.abs());
    (expected - actual).abs() <= tol
}

fn assert_phys_close(
    label: &str,
    expected: f64,
    actual: f64,
    phys_scale: f64,
    rtol: f64,
    atol_scale: f64,
) {
    let err = (expected - actual).abs();
    let tol = atol_scale * phys_scale + rtol * expected.abs().max(actual.abs());
    assert!(
        within_physical_tol(expected, actual, phys_scale, rtol, atol_scale),
        "{label}: expected={expected:.15e}, actual={actual:.15e}, err={err:.3e}, tol={tol:.3e} (phys_scale={phys_scale:.3e})"
    );
}

fn fixture_physical_scales(nc: usize, sumn: f64, t: f64, p: f64, v_tv: f64) -> Vec<f64> {
    let v_tp = v_tv / sumn;
    let rt = R * t;
    let p_scale = p.abs().max(sumn * rt / v_tv);
    let mut scales = Vec::new();
    // F, F_T, F_V, F_n
    scales.extend([sumn, sumn / t, sumn / v_tv]);
    scales.extend(std::iter::repeat_n(1.0, nc));
    // F_TT, F_VV, F_TV, F_Tn, F_Vn, F_nn
    scales.extend([sumn / (t * t), sumn / (v_tv * v_tv), sumn / (t * v_tv)]);
    scales.extend(std::iter::repeat_n(1.0 / t, nc));
    scales.extend(std::iter::repeat_n(1.0 / v_tv, nc));
    scales.extend(std::iter::repeat_n(1.0 / sumn, nc * nc));
    // P_tv, dP/dT, dP/dV, dP/dn
    scales.extend([p_scale, p_scale / t, p_scale / v_tv]);
    scales.extend(std::iter::repeat_n(p_scale / sumn, nc));
    // U_tv, S_tv, H_tv, sos_tv
    scales.extend([sumn * rt, sumn * R, sumn * rt, 100.0]);
    // Z, dZ/dT, dZ/dP, dZ/dn
    scales.extend([1.0, 1.0 / t, 1.0 / p]);
    scales.extend(std::iter::repeat_n(1.0, nc));
    // v_tp, dv/dT, dv/dP, dV/dn
    scales.extend([v_tp, v_tp / t, v_tp / p]);
    scales.extend(std::iter::repeat_n(v_tp, nc));
    // h_res, dh_res/dT, dh_res/dP, dH_res/dn
    scales.extend([rt, R, v_tp]);
    scales.extend(std::iter::repeat_n(rt, nc));
    // h, dh/dT, dh/dP, dH/dn
    scales.extend([rt, R, v_tp]);
    scales.extend(std::iter::repeat_n(rt, nc));
    // s_res, ds_res/dT, ds_res/dP, dS_res/dn
    scales.extend([R, R / t, v_tp / t]);
    scales.extend(std::iter::repeat_n(R, nc));
    // s, ds/dT, ds/dP, dS/dn
    scales.extend([R, R / t, v_tp / t]);
    scales.extend(std::iter::repeat_n(R, nc));
    // ln_phi, dlnphi/dT, dlnphi/dP, dlnphi/dn
    scales.extend(std::iter::repeat_n(1.0, nc));
    scales.extend(std::iter::repeat_n(1.0 / t, nc));
    scales.extend(std::iter::repeat_n(1.0 / p, nc));
    scales.extend(std::iter::repeat_n(1.0, nc * nc));
    // sos_tp
    scales.push(100.0);
    scales
}

#[test]
fn test_extcsp_unit_aware_tolerance_catches_zeroed_dv_dp() -> Result<()> {
    let model = ExtCsp::new(&["C3"])?;
    let t = 280.0;
    let p = 2.0e6;
    let vd = model.volume_tp_derivatives(t, p, &[1.0], Phase::Liquid)?;
    // Liquid dv/dP is ~ -1e-12 m^3/(mol Pa); unit-aware tolerance must reject 0.0!
    assert!(vd.dp < 0.0 && vd.dp.abs() < 1e-10);
    let phys_dv_dp = vd.value / p;
    assert!(
        !within_physical_tol(vd.dp, 0.0, phys_dv_dp, 1e-7, 1e-11),
        "unit-aware comparator must reject zeroed liquid dv/dP"
    );
    Ok(())
}

#[test]
fn test_extcsp_fixture_parity() -> Result<()> {
    let fixture = include_str!("fixtures/extcsp.csv");
    let mut rows = 0usize;
    for line in fixture.lines().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        let mut cols = line.split(';');
        let label = cols.next().unwrap();
        let expected: Vec<f64> = cols.map(|s| s.parse::<f64>().unwrap()).collect();

        let fields: Vec<_> = label.split('|').collect();
        assert_eq!(fields[0], "EXTCSP");
        let sh_eos = match fields[1] {
            "SRK" => CubicKind::SoaveRedlichKwong,
            _ => CubicKind::PengRobinson,
        };
        let sh_mix = match fields[2] {
            "HV1" => MixingRule::HuronVidal1,
            "HV2" => MixingRule::HuronVidal2,
            "NRTL" => MixingRule::Nrtl,
            "UNIFAC" => MixingRule::Unifac,
            "WS" => MixingRule::WongSandler,
            _ => MixingRule::VdW,
        };
        let sh_alpha = AlphaCorrelation::from_name(fields[3])?;
        let ref_eos = match fields[4] {
            "NIST_MEOS" => ExtCspReferenceKind::NistMeos,
            "MBWR19" => ExtCspReferenceKind::Mbwr(MbwrKind::Bender19),
            "MBWR32" => ExtCspReferenceKind::Mbwr(MbwrKind::Mbwr32),
            "SRK" => ExtCspReferenceKind::Cubic {
                kind: CubicKind::SoaveRedlichKwong,
                alpha: AlphaCorrelation::Classic,
            },
            _ => ExtCspReferenceKind::Cubic {
                kind: CubicKind::PengRobinson,
                alpha: AlphaCorrelation::Classic,
            },
        };
        let ref_comp = fields[5];
        let param_ref = fields[6];
        let names: Vec<_> = fields[7].split(',').collect();
        let t: f64 = fields[8].parse().unwrap();
        let p: f64 = fields[9].parse().unwrap();
        let v_tv: f64 = fields[10].parse().unwrap();
        let moles: Vec<f64> = fields[11].split(',').map(|v| v.parse().unwrap()).collect();
        let sumn: f64 = moles.iter().sum();
        let x: Vec<f64> = moles.iter().map(|&v| v / sumn).collect();
        let phase = if fields[12] == "L" {
            Phase::Liquid
        } else {
            Phase::Vapor
        };

        let model = ExtCsp::with_options(
            &names,
            ExtCspOptions {
                shape_eos: sh_eos,
                shape_mixing: sh_mix,
                shape_alpha: sh_alpha,
                reference_component: ref_comp,
                reference_eos: ref_eos,
                parameter_reference: param_ref,
                wong_sandler: None,
            },
        )?;

        let mut actual = Vec::with_capacity(expected.len());
        let f = model.residual_helmholtz_tv(t, v_tv, &moles)?;
        actual.extend([f.f, f.f_t, f.f_v]);
        actual.extend(&f.f_n);
        actual.extend([f.f_tt, f.f_vv, f.f_tv]);
        actual.extend(&f.f_tn);
        actual.extend(&f.f_vn);
        actual.extend(&f.f_nn);
        let ptv = model.pressure_tv(t, v_tv, &moles)?;
        actual.extend([ptv.value, ptv.dt, ptv.dv]);
        actual.extend(&ptv.dn);
        let cal = model.caloric_tv(t, v_tv, &moles)?;
        let sostv = model.speed_of_sound_tv(t, v_tv, &moles)?;
        actual.extend([cal.internal_energy, cal.entropy, cal.enthalpy, sostv]);
        let z = model.compressibility_factor_derivatives(t, p, &x, phase)?;
        actual.extend([z.value, z.dt, z.dp]);
        actual.extend(&z.dn);
        let v = model.volume_tp_derivatives(t, p, &x, phase)?;
        actual.extend([v.value, v.dt, v.dp]);
        actual.extend(&v.dn);
        let hr = model.residual_enthalpy_derivatives(t, p, &x, phase)?;
        actual.extend([hr.value, hr.dt, hr.dp]);
        actual.extend(&hr.dn);
        let ht = model.enthalpy_derivatives(t, p, &x, phase)?;
        actual.extend([ht.value, ht.dt, ht.dp]);
        actual.extend(&ht.dn);
        let sr = model.residual_entropy_derivatives(t, p, &x, phase)?;
        actual.extend([sr.value, sr.dt, sr.dp]);
        actual.extend(&sr.dn);
        let st = model.entropy_derivatives(t, p, &x, phase)?;
        actual.extend([st.value, st.dt, st.dp]);
        actual.extend(&st.dn);
        let fug = model.ln_fugacity_derivatives_tp(t, p, &x, phase)?;
        actual.extend(&fug.ln_phi);
        actual.extend(&fug.dt);
        actual.extend(&fug.dp);
        actual.extend(&fug.dn);
        actual.push(model.speed_of_sound_tp(t, p, &x, phase)?);

        let scales = fixture_physical_scales(moles.len(), sumn, t, p, v_tv);
        assert_eq!(actual.len(), expected.len(), "length mismatch for {label}");
        assert_eq!(
            scales.len(),
            expected.len(),
            "scale length mismatch for {label}"
        );
        for (idx, ((&a, &b), &phys)) in expected.iter().zip(actual.iter()).zip(&scales).enumerate()
        {
            assert_phys_close(&format!("{label}[{idx}]"), a, b, phys, 1e-7, 1e-11);
        }
        rows += 1;
    }
    assert_eq!(rows, 24, "expected exact 24 ExtCSP fixture rows");
    Ok(())
}

fn verify_helmholtz_and_tp_fd(
    model: &ExtCsp,
    t: f64,
    p: f64,
    n: &[f64],
    phase: Phase,
) -> Result<()> {
    let nc = n.len();
    let sumn: f64 = n.iter().sum();
    let v = model.volume_tp(t, p, n, phase)?;
    let ptv = model.pressure_tv(t, v, n)?;
    assert_phys_close("TP/TV pressure closure", p, ptv.value, p, 1e-9, 1e-12);

    let f = model.residual_helmholtz_tv(t, v, n)?;

    // 1. Central finite differences of all 10 Helmholtz blocks with unit-aware scales
    let dt = 1e-5 * t;
    let dv = 1e-5 * v;
    let fp_t = model.residual_helmholtz_tv(t + dt, v, n)?;
    let fm_t = model.residual_helmholtz_tv(t - dt, v, n)?;
    let fp_v = model.residual_helmholtz_tv(t, v + dv, n)?;
    let fm_v = model.residual_helmholtz_tv(t, v - dv, n)?;

    assert_phys_close(
        "F_T",
        f.f_t,
        (fp_t.f - fm_t.f) / (2.0 * dt),
        sumn / t,
        1e-5,
        1e-7,
    );
    assert_phys_close(
        "F_V",
        f.f_v,
        (fp_v.f - fm_v.f) / (2.0 * dv),
        sumn / v,
        1e-5,
        1e-7,
    );
    assert_phys_close(
        "F_TT",
        f.f_tt,
        (fp_t.f_t - fm_t.f_t) / (2.0 * dt),
        sumn / (t * t),
        1e-5,
        1e-7,
    );
    assert_phys_close(
        "F_VV",
        f.f_vv,
        (fp_v.f_v - fm_v.f_v) / (2.0 * dv),
        sumn / (v * v),
        1e-5,
        1e-7,
    );
    assert_phys_close(
        "F_TV",
        f.f_tv,
        (fp_t.f_v - fm_t.f_v) / (2.0 * dt),
        sumn / (t * v),
        1e-5,
        1e-7,
    );

    let dn = 1e-6 * sumn;
    for j in 0..nc {
        let mut np = n.to_vec();
        let mut nm = n.to_vec();
        np[j] += dn;
        nm[j] -= dn;
        let fp = model.residual_helmholtz_tv(t, v, &np)?;
        let fm = model.residual_helmholtz_tv(t, v, &nm)?;
        assert_phys_close("F_n", f.f_n[j], (fp.f - fm.f) / (2.0 * dn), 1.0, 1e-5, 1e-7);
        assert_phys_close(
            "F_Tn",
            f.f_tn[j],
            (fp.f_t - fm.f_t) / (2.0 * dn),
            1.0 / t,
            1e-5,
            1e-7,
        );
        assert_phys_close(
            "F_Vn",
            f.f_vn[j],
            (fp.f_v - fm.f_v) / (2.0 * dn),
            1.0 / v,
            1e-5,
            1e-7,
        );
        for i in 0..nc {
            assert_phys_close(
                "F_nn",
                f.f_nn[i * nc + j],
                (fp.f_n[i] - fm.f_n[i]) / (2.0 * dn),
                1.0 / sumn,
                1e-5,
                1e-7,
            );
            assert_phys_close(
                "F_nn symmetry",
                f.f_nn[i * nc + j],
                f.f_nn[j * nc + i],
                1.0 / sumn,
                1e-12,
                1e-12,
            );
        }
    }

    // 2. Euler homogeneity of F(T, V, n): sum(n_i F_n[i]) + V F_V = F
    let euler_f: f64 = n.iter().zip(&f.f_n).map(|(ni, fni)| ni * fni).sum::<f64>() + v * f.f_v;
    assert_phys_close("Euler F", euler_f, f.f, sumn, 1e-10, 1e-10);
    for i in 0..nc {
        let euler_fn: f64 = (0..nc).map(|j| n[j] * f.f_nn[i * nc + j]).sum::<f64>() + v * f.f_vn[i];
        assert_phys_close("Euler F_n", euler_fn, 0.0, 1.0, 1e-10, 1e-10);
    }

    // 3. TP derivatives and thermodynamic identities
    let zd = model.compressibility_factor_derivatives(t, p, n, phase)?;
    let vd = model.volume_tp_derivatives(t, p, n, phase)?;
    let svd = model.specific_volume_derivatives(t, p, n, phase)?;
    let hrd = model.residual_enthalpy_derivatives(t, p, n, phase)?;
    let hd = model.enthalpy_derivatives(t, p, n, phase)?;
    let srd = model.residual_entropy_derivatives(t, p, n, phase)?;
    let sd = model.entropy_derivatives(t, p, n, phase)?;
    let grd = model.residual_gibbs_derivatives(t, p, n, phase)?;
    let gd = model.gibbs_derivatives(t, p, n, phase)?;
    let fug = model.ln_fugacity_derivatives_tp(t, p, n, phase)?;
    let (cp, cv) = model.heat_capacities_tp(t, p, n, phase)?;
    let sos_tp = model.speed_of_sound_tp(t, p, n, phase)?;
    let sos_tv = model.speed_of_sound_tv(t, v, n)?;

    assert_phys_close(
        "specific_volume == svd.value",
        model.specific_volume(t, p, n, phase)?,
        svd.value,
        v / sumn,
        1e-13,
        1e-13,
    );
    assert_phys_close("volume_tp == vd.value", v, vd.value, v, 1e-13, 1e-13);
    assert_phys_close(
        "gibbs == gd.value",
        model.gibbs(t, p, n, phase)?,
        gd.value,
        sumn * R * t,
        1e-12,
        1e-12,
    );
    assert_phys_close(
        "residual_gibbs == grd.value",
        model.residual_gibbs(t, p, n, phase)?,
        grd.value,
        sumn * R * t,
        1e-12,
        1e-12,
    );

    assert_phys_close("dS/dP = -dV/dT", sd.dp, -vd.dt, v / t, 1e-10, 1e-10);
    assert_phys_close(
        "dH/dP = V - T*dV/dT",
        hd.dp,
        vd.value - t * vd.dt,
        v,
        1e-10,
        1e-10,
    );
    assert_phys_close("dH/dT = T*dS/dT", hd.dt, t * sd.dt, sumn * R, 1e-10, 1e-10);
    assert_phys_close("Cp = dH/dT", cp, hd.dt, sumn * R, 1e-10, 1e-10);
    assert!(cp > cv && cv > 0.0);
    assert_phys_close("sos_tp == sos_tv", sos_tp, sos_tv, 100.0, 1e-11, 1e-11);
    assert_phys_close(
        "G = H - T*S",
        gd.value,
        hd.value - t * sd.value,
        sumn * R * t,
        1e-10,
        1e-10,
    );
    assert_phys_close("dG/dT = -S", gd.dt, -sd.value, sumn * R, 1e-10, 1e-10);
    assert_phys_close("dG/dP = V", gd.dp, vd.value, v, 1e-10, 1e-10);
    assert_phys_close(
        "dG_res/dT = -S_res",
        grd.dt,
        -srd.value,
        sumn * R,
        1e-10,
        1e-10,
    );
    for i in 0..nc {
        assert_phys_close(
            "dG_res/dn = R*T*lnphi",
            grd.dn[i],
            R * t * fug.ln_phi[i],
            R * t,
            1e-10,
            1e-10,
        );
        assert_phys_close(
            "dlnphi/dT = -H_bar_res/(RT^2)",
            fug.dt[i],
            -hrd.dn[i] / (R * t * t),
            1.0 / t,
            1e-10,
            1e-10,
        );
        assert_phys_close(
            "dlnphi/dP = V_bar/(RT) - 1/P",
            fug.dp[i],
            vd.dn[i] / (R * t) - 1.0 / p,
            1.0 / p,
            1e-10,
            1e-10,
        );
        let mut gd_sum = 0.0;
        for j in 0..nc {
            gd_sum += n[j] * fug.dn[i * nc + j];
            assert_phys_close(
                "dlnphi/dn symmetry",
                fug.dn[i * nc + j],
                fug.dn[j * nc + i],
                1.0 / sumn,
                1e-11,
                1e-11,
            );
        }
        assert_phys_close("Gibbs-Duhem", gd_sum, 0.0, 1.0, 1e-10, 1e-10);
    }

    // 4. Unit-aware TP FD for Z, V, v, H, S, lnphi w.r.t. T, P, n
    let dp = 1e-5 * p;
    assert_phys_close(
        "dZ/dT FD",
        zd.dt,
        (model.compressibility_factor(t + dt, p, n, phase)?
            - model.compressibility_factor(t - dt, p, n, phase)?)
            / (2.0 * dt),
        1.0 / t,
        1e-5,
        1e-7,
    );
    assert_phys_close(
        "dZ/dP FD",
        zd.dp,
        (model.compressibility_factor(t, p + dp, n, phase)?
            - model.compressibility_factor(t, p - dp, n, phase)?)
            / (2.0 * dp),
        1.0 / p,
        1e-5,
        1e-7,
    );
    assert_phys_close(
        "dV/dP FD",
        vd.dp,
        (model.volume_tp(t, p + dp, n, phase)? - model.volume_tp(t, p - dp, n, phase)?)
            / (2.0 * dp),
        v / p,
        1e-5,
        1e-7,
    );
    assert_phys_close(
        "dv/dP FD",
        svd.dp,
        (model.specific_volume(t, p + dp, n, phase)?
            - model.specific_volume(t, p - dp, n, phase)?)
            / (2.0 * dp),
        (v / sumn) / p,
        1e-5,
        1e-7,
    );
    assert_phys_close(
        "dH/dT FD",
        hd.dt,
        (model.enthalpy(t + dt, p, n, phase)? - model.enthalpy(t - dt, p, n, phase)?) / (2.0 * dt),
        sumn * R,
        1e-5,
        1e-7,
    );
    assert_phys_close(
        "dS/dP FD",
        sd.dp,
        (model.entropy(t, p + dp, n, phase)? - model.entropy(t, p - dp, n, phase)?) / (2.0 * dp),
        v / t,
        1e-5,
        1e-7,
    );
    for j in 0..nc {
        let mut np = n.to_vec();
        let mut nm = n.to_vec();
        np[j] += dn;
        nm[j] -= dn;
        let lnp_p = model.ln_fugacity_tp(t, p, &np, phase)?;
        let lnp_m = model.ln_fugacity_tp(t, p, &nm, phase)?;
        for i in 0..nc {
            assert_phys_close(
                "dlnphi/dn FD",
                fug.dn[i * nc + j],
                (lnp_p[i] - lnp_m[i]) / (2.0 * dn),
                1.0 / sumn,
                1e-5,
                1e-7,
            );
        }
    }
    Ok(())
}

#[test]
fn test_extcsp_tcpr_reference_untranslated_consistency() -> Result<()> {
    // Regression test: "tcPR" auto-enables volume translation inside Cubic::with_options,
    // so ExtCSP must explicitly disable translation on both shape and cubic reference EOS
    // so volume_tp and pure_residual_helmholtz_tv stay 100% consistent!
    let model = ExtCsp::with_options(
        &["C3"],
        ExtCspOptions {
            shape_eos: CubicKind::PengRobinson,
            shape_alpha: AlphaCorrelation::Twu,
            reference_component: "C3",
            reference_eos: ExtCspReferenceKind::Cubic {
                kind: CubicKind::PengRobinson,
                alpha: AlphaCorrelation::Twu,
            },
            parameter_reference: "tcPR",
            ..Default::default()
        },
    )?;
    assert!(!model.shape_eos().volume_translation_enabled());
    verify_helmholtz_and_tp_fd(&model, 300.0, 2.5e6, &[1.0], Phase::Liquid)?;
    verify_helmholtz_and_tp_fd(&model, 350.0, 2.0e6, &[2.5], Phase::Vapor)?;
    Ok(())
}

#[test]
fn test_extcsp_active_twu_alpha_and_amount_invariance() -> Result<()> {
    // Verify genuine Twu alpha (no Classic fallback!) using "tcPR" database parameters.
    let model = ExtCsp::with_options(
        &["C1", "C3"],
        ExtCspOptions {
            shape_eos: CubicKind::PengRobinson,
            shape_alpha: AlphaCorrelation::Twu,
            parameter_reference: "tcPR",
            ..Default::default()
        },
    )?;
    for idx in 0..2 {
        let a = model.shape_eos().alpha(idx)?;
        assert_eq!(a.correlation, AlphaCorrelation::Twu);
        assert!(
            !a.used_fallback,
            "component {idx} unexpectedly fell back from Twu"
        );
    }
    let ref_alpha = model.reference_alpha();
    assert_eq!(ref_alpha.correlation, AlphaCorrelation::Twu);
    assert!(
        !ref_alpha.used_fallback,
        "reference C3 unexpectedly fell back from Twu"
    );

    // Independently check the implicit Twu equation: F * alpha_Twu(N*T/F; L,M,N_exp) == D / (H * ac0)
    let t = 300.0;
    let p = 2.5e6;
    let base_n = [0.4, 0.6];
    let (h_base, f_base) = model.shape_factors(t, &base_n)?;
    let [l_param, m_param, n_exp] = ref_alpha.parameters;
    let tc0 = model.reference_critical_temperature();
    let tr0 = (1.0 * t / f_base) / tc0;
    let alpha_twu_indep =
        tr0.powf(n_exp * (m_param - 1.0)) * (l_param * (1.0 - tr0.powf(n_exp * m_param))).exp();
    let prep = model.shape_eos().prepare(t, &base_n)?;
    let b_molar = h_base * model.reference_bc0();
    let v_probe = 1.0e-2;
    let p_pr = prep.pressure(v_probe)?.value;
    let a_molar = (R * t / (v_probe - b_molar) - p_pr)
        * (v_probe * v_probe + 2.0 * b_molar * v_probe - b_molar * b_molar);
    let rhs = a_molar / (h_base * model.reference_ac0());
    assert_phys_close(
        "independent Twu implicit equation",
        f_base * alpha_twu_indep,
        rhs,
        rhs,
        1e-13,
        1e-14,
    );

    // Exact scale invariance of intensive shape factors (H/N, F/N) across 18 orders of magnitude
    for scale in [1.0e-12_f64, 0.37, 2.5, 1.0e6] {
        let scaled = [base_n[0] * scale, base_n[1] * scale];
        let (h_s, f_s) = model.shape_factors(t, &scaled)?;
        assert_phys_close(
            &format!("H/N invariance at scale {scale:e}"),
            h_s / scale,
            h_base,
            h_base,
            1e-14,
            1e-15,
        );
        assert_phys_close(
            &format!("F/N invariance at scale {scale:e}"),
            f_s / scale,
            f_base,
            f_base,
            1e-14,
            1e-15,
        );
    }

    // Full Helmholtz + TP FD at non-unit moles
    verify_helmholtz_and_tp_fd(&model, t, p, &[1.2, 1.8], Phase::Vapor)?;
    verify_helmholtz_and_tp_fd(&model, 240.0, 5.0e6, &[0.8, 1.2], Phase::Liquid)?;
    Ok(())
}

#[test]
fn test_extcsp_volume_vs_specific_volume_derivatives() -> Result<()> {
    let model = ExtCsp::new(&["C1", "C3"])?;
    let t = 320.0;
    let p = 2.0e6;
    let n = [1.2, 0.8];
    let sumn = 2.0;
    let phase = Phase::Vapor;

    let v_ext = model.volume_tp(t, p, &n, phase)?;
    let v_mol = model.specific_volume(t, p, &n, phase)?;
    let vd = model.volume_tp_derivatives(t, p, &n, phase)?;
    let svd = model.specific_volume_derivatives(t, p, &n, phase)?;

    assert_phys_close("v_ext == vd.value", v_ext, vd.value, v_ext, 1e-14, 1e-14);
    assert_phys_close("v_mol == svd.value", v_mol, svd.value, v_mol, 1e-14, 1e-14);
    assert_phys_close(
        "vd.value == N * svd.value",
        vd.value,
        sumn * svd.value,
        v_ext,
        1e-14,
        1e-14,
    );

    // Euler homogeneity: sum(n_i * dV/dn_i) = V (degree 1), sum(n_i * dv/dn_i) = 0 (degree 0)
    let euler_ext = n[0] * vd.dn[0] + n[1] * vd.dn[1];
    let euler_mol = n[0] * svd.dn[0] + n[1] * svd.dn[1];
    assert_phys_close("Euler V", euler_ext, v_ext, v_ext, 1e-12, 1e-12);
    assert_phys_close("Euler v", euler_mol, 0.0, v_mol, 1e-12, 1e-12);

    // Central FD of both extensive volume_tp and molar specific_volume w.r.t. T, P, n
    let dt = 1e-5 * t;
    let dp = 1e-5 * p;
    let dn = 1e-6 * sumn;
    assert_phys_close(
        "dV/dT FD",
        vd.dt,
        (model.volume_tp(t + dt, p, &n, phase)? - model.volume_tp(t - dt, p, &n, phase)?)
            / (2.0 * dt),
        v_ext / t,
        1e-6,
        1e-8,
    );
    assert_phys_close(
        "dv/dT FD",
        svd.dt,
        (model.specific_volume(t + dt, p, &n, phase)?
            - model.specific_volume(t - dt, p, &n, phase)?)
            / (2.0 * dt),
        v_mol / t,
        1e-6,
        1e-8,
    );
    assert_phys_close(
        "dV/dP FD",
        vd.dp,
        (model.volume_tp(t, p + dp, &n, phase)? - model.volume_tp(t, p - dp, &n, phase)?)
            / (2.0 * dp),
        v_ext / p,
        1e-6,
        1e-8,
    );
    assert_phys_close(
        "dv/dP FD",
        svd.dp,
        (model.specific_volume(t, p + dp, &n, phase)?
            - model.specific_volume(t, p - dp, &n, phase)?)
            / (2.0 * dp),
        v_mol / p,
        1e-6,
        1e-8,
    );
    for j in 0..2 {
        let mut np = n;
        let mut nm = n;
        np[j] += dn;
        nm[j] -= dn;
        assert_phys_close(
            "dV/dn FD",
            vd.dn[j],
            (model.volume_tp(t, p, &np, phase)? - model.volume_tp(t, p, &nm, phase)?) / (2.0 * dn),
            v_mol,
            1e-6,
            1e-8,
        );
        assert_phys_close(
            "dv/dn FD",
            svd.dn[j],
            (model.specific_volume(t, p, &np, phase)? - model.specific_volume(t, p, &nm, phase)?)
                / (2.0 * dn),
            v_mol / sumn,
            1e-6,
            1e-8,
        );
    }
    Ok(())
}

/// Independently recover untranslated two-parameter cubic `(a_mix, b_mix)` from
/// standalone `Cubic::pressure(v)` at two moderate molar volumes without calling
/// any `ExtCsp` method.
fn recover_cubic_ab_from_pressure(
    cubic: &Cubic,
    kind: CubicKind,
    t: f64,
    x: &[f64],
) -> Result<(f64, f64)> {
    let (d1, d2) = match kind {
        CubicKind::SoaveRedlichKwong => (1.0, 0.0),
        CubicKind::PengRobinson => (1.0 + 2.0_f64.sqrt(), 1.0 - 2.0_f64.sqrt()),
        _ => unreachable!(),
    };
    let prep = cubic.prepare(t, x)?;
    let v1 = 5.0e-3;
    let v2 = 2.0e-2;
    let p1 = prep.pressure(v1)?.value;
    let p2 = prep.pressure(v2)?.value;
    let rt = R * t;
    let a_from_v = |v: f64, p: f64, b: f64| (rt / (v - b) - p) * (v + d1 * b) * (v + d2 * b);
    let residual = |b: f64| a_from_v(v1, p1, b) - a_from_v(v2, p2, b);

    // Bracket b inside (0, v1 * 0.5)
    let mut lo = 1.0e-7;
    let mut hi = 1.0e-3;
    let r_lo = residual(lo);
    while residual(hi) * r_lo > 0.0 && hi < v1 * 0.8 {
        hi *= 1.5;
    }
    assert!(
        residual(lo) * residual(hi) <= 0.0,
        "failed to bracket cubic covolume b from two-point pressure"
    );
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        let rm = residual(mid);
        if rm == 0.0 || (hi - lo) <= 1.0e-15 * mid {
            lo = mid;
            hi = mid;
            break;
        }
        if rm * r_lo <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let b_mix = 0.5 * (lo + hi);
    let a_mix = a_from_v(v1, p1, b_mix);
    Ok((a_mix, b_mix))
}

/// Independently compute ExtCSP shape factors `(H_expected, F_expected)` directly
/// from standalone untranslated `Cubic` parameters and reference `(tc0, pc0)`,
/// without calling any method on `ExtCsp`.
#[allow(clippy::too_many_arguments)]
fn independent_shape_factors(
    shape_kind: CubicKind,
    shape_mixing: MixingRule,
    shape_alpha: AlphaCorrelation,
    param_ref: &str,
    ws_data: Option<&thermopack_rs::wong_sandler::WongSandlerData>,
    comps: &[&str],
    ref_comp: &str,
    tc0: f64,
    pc0: f64,
    t: f64,
    n: &[f64],
) -> Result<(f64, f64)> {
    let sumn: f64 = n.iter().sum();
    let x: Vec<f64> = n.iter().map(|&ni| ni / sumn).collect();

    let mut mix_cubic = Cubic::with_options(
        shape_kind,
        comps,
        thermopack_rs::CubicOptions {
            alpha: Some(shape_alpha),
            mixing_rule: shape_mixing,
            parameter_reference: param_ref,
            volume_translation: false,
            wong_sandler: ws_data,
        },
    )?;
    mix_cubic.disable_volume_translation();
    let (a_mix, b_mix) = recover_cubic_ab_from_pressure(&mix_cubic, shape_kind, t, &x)?;

    let mut ref_cubic = Cubic::with_options(
        shape_kind,
        &[ref_comp],
        thermopack_rs::CubicOptions {
            alpha: Some(shape_alpha),
            mixing_rule: MixingRule::VdW,
            parameter_reference: param_ref,
            volume_translation: false,
            ..Default::default()
        },
    )?;
    ref_cubic.disable_volume_translation();
    let ref_c = ref_cubic.components()[0];
    let (a_ref_tc, b_ref_tc) =
        recover_cubic_ab_from_pressure(&ref_cubic, shape_kind, ref_c.critical_temperature, &[1.0])?;
    let omega_b = b_ref_tc * ref_c.critical_pressure / (R * ref_c.critical_temperature);
    let omega_a = a_ref_tc * ref_c.critical_pressure
        / (R * ref_c.critical_temperature * R * ref_c.critical_temperature);

    let bc0 = omega_b * R * tc0 / pc0;
    let ac0 = omega_a * (R * tc0) * (R * tc0) / pc0;
    let h_molar = b_mix / bc0;
    let h_expected = sumn * h_molar;

    // Solve a_mix / b_mix - a0(T0) / bc0 = 0, i.e., f * alpha0(T / f) = a_mix / (h_molar * ac0)
    let ref_alpha = *ref_cubic.alpha(0)?;
    let target = a_mix / (h_molar * ac0);
    let g = |f_val: f64| -> Result<f64> {
        let av = ref_alpha.evaluate(t / f_val, tc0)?;
        Ok(f_val * av.value - target)
    };
    let mut lo = 1.0e-4;
    let mut hi = 1.0e2;
    while g(lo)? > 0.0 {
        lo *= 0.5;
    }
    while g(hi)? < 0.0 {
        hi *= 2.0;
    }
    for _ in 0..120 {
        let mid = 0.5 * (lo + hi);
        let gm = g(mid)?;
        if gm.abs() <= 1.0e-15 * target || (hi - lo) <= 1.0e-15 * mid {
            lo = mid;
            hi = mid;
            break;
        }
        if gm < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let f_expected = sumn * 0.5 * (lo + hi);
    Ok((h_expected, f_expected))
}

#[test]
fn test_extcsp_all_backends_wong_sandler_zero_components_and_minimum_gibbs() -> Result<()> {
    // 1. NIST MEOS + VdW (classic)
    let m_nist = ExtCsp::new(&["C1", "C2", "C3"])?;
    verify_helmholtz_and_tp_fd(&m_nist, 310.0, 4.0e6, &[1.0, 0.6, 0.4], Phase::Vapor)?;
    verify_helmholtz_and_tp_fd(&m_nist, 230.0, 6.0e6, &[0.5, 0.3, 0.2], Phase::Liquid)?;

    // 2. Cubic SRK & PR reference backends + truly independent (non-circular) mapped checks
    let m_cb_srk = ExtCsp::with_options(
        &["C1", "C3"],
        ExtCspOptions {
            shape_eos: CubicKind::SoaveRedlichKwong,
            reference_eos: ExtCspReferenceKind::Cubic {
                kind: CubicKind::SoaveRedlichKwong,
                alpha: AlphaCorrelation::Classic,
            },
            ..Default::default()
        },
    )?;
    verify_helmholtz_and_tp_fd(&m_cb_srk, 320.0, 3.0e6, &[0.6, 0.4], Phase::Vapor)?;

    let m_cb_pr = ExtCsp::with_options(
        &["CO2", "C1"],
        ExtCspOptions {
            shape_eos: CubicKind::PengRobinson,
            reference_eos: ExtCspReferenceKind::Cubic {
                kind: CubicKind::PengRobinson,
                alpha: AlphaCorrelation::Classic,
            },
            ..Default::default()
        },
    )?;
    verify_helmholtz_and_tp_fd(&m_cb_pr, 250.0, 6.0e6, &[0.4, 0.6], Phase::Liquid)?;

    // Truly independent mapped reference-state checks for both Cubic SRK and PR backends
    // (H, F recovered from standalone Cubic pressure without calling ExtCsp::shape_factors):
    for (label, kind, comps, model, t, p, n, ph) in [
        (
            "Cubic SRK",
            CubicKind::SoaveRedlichKwong,
            &["C1", "C3"][..],
            &m_cb_srk,
            320.0,
            3.0e6,
            [1.2, 0.8],
            Phase::Vapor,
        ),
        (
            "Cubic PR",
            CubicKind::PengRobinson,
            &["CO2", "C1"][..],
            &m_cb_pr,
            250.0,
            6.0e6,
            [0.8, 1.2],
            Phase::Liquid,
        ),
    ] {
        let sumn = n[0] + n[1];
        let ref_cubic = Cubic::new(kind, &["C3"])?;
        let c3 = ref_cubic.components()[0];
        let (h_exp, f_exp) = independent_shape_factors(
            kind,
            MixingRule::VdW,
            AlphaCorrelation::Classic,
            "DEFAULT",
            None,
            comps,
            "C3",
            c3.critical_temperature,
            c3.critical_pressure,
            t,
            &n,
        )?;
        let (h_act, f_act) = model.shape_factors(t, &n)?;
        assert_phys_close(
            &format!("{label} independent H"),
            h_exp,
            h_act,
            h_exp,
            1e-11,
            1e-12,
        );
        assert_phys_close(
            &format!("{label} independent F"),
            f_exp,
            f_act,
            f_exp,
            1e-11,
            1e-12,
        );
        let t0 = sumn * t / f_exp;
        let p0 = h_exp * p / f_exp;
        let ref_state = ref_cubic.prepare(t0, &[1.0])?;
        let v0_ref = ref_state.volume(p0, ph)?;
        let v_exp = h_exp * v0_ref;
        let v_csp = model.volume_tp(t, p, &n, ph)?;
        assert_phys_close(
            &format!("{label} independent mapped V == H*v0"),
            v_csp,
            v_exp,
            v_exp,
            1e-11,
            1e-12,
        );
        let f_ref = ref_state.residual(v0_ref)?.helmholtz_tv / (R * t0);
        let f_csp = model.residual_helmholtz_tv(t, v_exp, &n)?;
        assert_phys_close(
            &format!("{label} independent mapped F_res"),
            f_csp.f,
            sumn * f_ref,
            sumn,
            1e-11,
            1e-12,
        );
        let p_ref = ref_state.pressure(v0_ref)?.value;
        let p_csp = model.pressure_tv(t, v_exp, &n)?.value;
        assert_phys_close(
            &format!("{label} independent mapped P == (F/H)*P0"),
            p_csp,
            (f_exp / h_exp) * p_ref,
            p,
            1e-11,
            1e-12,
        );
    }

    // 3. MBWR19 (Bender) & MBWR32 reference backends + truly independent mapped checks
    let m_mbwr19 = ExtCsp::with_options(
        &["C1", "C3"],
        ExtCspOptions {
            reference_eos: ExtCspReferenceKind::Mbwr(MbwrKind::Bender19),
            ..Default::default()
        },
    )?;
    verify_helmholtz_and_tp_fd(&m_mbwr19, 320.0, 2.5e6, &[0.5, 0.5], Phase::Vapor)?;

    let m_mbwr32 = ExtCsp::with_options(
        &["C1", "C3"],
        ExtCspOptions {
            reference_eos: ExtCspReferenceKind::Mbwr(MbwrKind::Mbwr32),
            ..Default::default()
        },
    )?;
    verify_helmholtz_and_tp_fd(&m_mbwr32, 240.0, 5.0e6, &[0.5, 0.5], Phase::Liquid)?;

    for (kind, model, t, p, ph) in [
        (MbwrKind::Bender19, &m_mbwr19, 320.0, 2.5e6, Phase::Vapor),
        (MbwrKind::Mbwr32, &m_mbwr32, 240.0, 5.0e6, Phase::Liquid),
    ] {
        let n = [0.75, 0.75];
        let sumn = 1.5;
        let mbwr = Mbwr::new(kind, "C3")?;
        let dp_ref = mbwr.source_density_parameters()?;
        let (h_exp, f_exp) = independent_shape_factors(
            CubicKind::SoaveRedlichKwong,
            MixingRule::VdW,
            AlphaCorrelation::Classic,
            "DEFAULT",
            None,
            &["C1", "C3"],
            "C3",
            dp_ref.critical_temperature,
            dp_ref.critical_pressure,
            t,
            &n,
        )?;
        let (h_act, f_act) = model.shape_factors(t, &n)?;
        assert_phys_close("MBWR independent H", h_exp, h_act, h_exp, 1e-11, 1e-12);
        assert_phys_close("MBWR independent F", f_exp, f_act, f_exp, 1e-11, 1e-12);
        let t0 = sumn * t / f_exp;
        let p0 = h_exp * p / f_exp;
        let iter = mbwr.prepare_baseline_density(t0, dp_ref)?.density(p0, ph)?;
        let v_exp = h_exp / iter.density;
        let v_csp = model.volume_tp(t, p, &n, ph)?;
        assert_phys_close(
            "MBWR independent mapped V",
            v_csp,
            v_exp,
            v_exp,
            1e-10,
            1e-11,
        );
        let f_csp = model.residual_helmholtz_tv(t, v_exp, &n)?;
        let f_ref = mbwr.prepare(t0)?.residual_helmholtz(v_exp / h_exp, 1.0)?;
        assert_phys_close(
            "MBWR independent mapped F_res",
            f_csp.f,
            sumn * f_ref.value,
            sumn,
            1e-11,
            1e-12,
        );
        let p_csp = model.pressure_tv(t, v_exp, &n)?.value;
        let p_ref = mbwr.prepare(t0)?.pressure_tv(v_exp / h_exp, 1.0)?.value;
        assert_phys_close(
            "MBWR independent mapped P",
            p_csp,
            (f_exp / h_exp) * p_ref,
            p,
            1e-10,
            1e-11,
        );
    }

    // 4. Wong-Sandler mixing rule: verify non-zero H_T and independent mapped state!
    let mut ws = thermopack_rs::wong_sandler::WongSandlerData::new(2);
    ws.set_kij(0, 1, thermopack_rs::wong_sandler::Fraction::constant(0.1))?;
    ws.set_kij(1, 0, thermopack_rs::wong_sandler::Fraction::constant(0.1))?;
    ws.set_tauij(0, 1, thermopack_rs::wong_sandler::Fraction::constant(0.25))?;
    ws.set_tauij(1, 0, thermopack_rs::wong_sandler::Fraction::constant(-0.15))?;
    ws.set_alphaij(0, 1, 0.3)?;
    ws.set_alphaij(1, 0, 0.3)?;
    let m_ws = ExtCsp::with_options(
        &["CO2", "C1"],
        ExtCspOptions {
            shape_eos: CubicKind::PengRobinson,
            shape_mixing: MixingRule::WongSandler,
            wong_sandler: Some(&ws),
            ..Default::default()
        },
    )?;
    let ws_t = 320.0;
    let ws_p = 3.0e6;
    let ws_n = [0.5, 0.5];
    let meos_c3 = Multiparameter::new("C3")?;
    let (h_ws_exp, f_ws_exp) = independent_shape_factors(
        CubicKind::PengRobinson,
        MixingRule::WongSandler,
        AlphaCorrelation::Classic,
        "DEFAULT",
        Some(&ws),
        &["CO2", "C1"],
        "C3",
        meos_c3.tc(),
        meos_c3.pc(),
        ws_t,
        &ws_n,
    )?;
    let (h_ws_act, f_ws_act) = m_ws.shape_factors(ws_t, &ws_n)?;
    assert_phys_close(
        "WS independent H",
        h_ws_exp,
        h_ws_act,
        h_ws_exp,
        1e-11,
        1e-12,
    );
    assert_phys_close(
        "WS independent F",
        f_ws_exp,
        f_ws_act,
        f_ws_exp,
        1e-11,
        1e-12,
    );
    let ws_v = m_ws.volume_tp(ws_t, ws_p, &ws_n, Phase::Vapor)?;
    let f_ws_csp = m_ws.residual_helmholtz_tv(ws_t, ws_v, &ws_n)?;
    let f_ws_ref = meos_c3
        .prepare(ws_t / f_ws_exp)?
        .residual_helmholtz(ws_v / h_ws_exp, 1.0)?;
    assert_phys_close(
        "WS independent mapped F_res",
        f_ws_csp.f,
        f_ws_ref.value,
        1.0,
        1e-11,
        1e-12,
    );

    let ws_coeffs = m_ws.shape_eos().wong_sandler_coefficients(ws_t, &ws_n)?;
    let h_t_expected = ws_coeffs.b_t / m_ws.reference_bc0();
    assert!(
        h_t_expected.abs() > 1e-7,
        "Wong-Sandler H_T must be genuinely non-zero (got {h_t_expected:.3e})"
    );
    let dt = 1e-4 * ws_t;
    let (h_p, _) = m_ws.shape_factors(ws_t + dt, &ws_n)?;
    let (h_m, _) = m_ws.shape_factors(ws_t - dt, &ws_n)?;
    assert_phys_close(
        "Wong-Sandler H_T vs shape_factors FD",
        h_t_expected,
        (h_p - h_m) / (2.0 * dt),
        h_t_expected.abs(),
        1e-6,
        1e-8,
    );
    verify_helmholtz_and_tp_fd(&m_ws, ws_t, ws_p, &ws_n, Phase::Vapor)?;

    // 4b. Independent fixed-V finite-difference Cv check across all backends & mixing rules:
    //     Cv = (dU/dT)_{V,n} where U(T,V,n) = sum n_i (h_id_i(T) - R T) - R T^2 * F_T(T,V,n).
    let internal_energy_from_helmholtz =
        |model: &ExtCsp, t_eval: f64, v_fixed: f64, n_vec: &[f64]| -> Result<f64> {
            let f = model.residual_helmholtz_tv(t_eval, v_fixed, n_vec)?;
            let mut u_id = 0.0;
            for (i, &ni) in n_vec.iter().enumerate() {
                if ni > 0.0 {
                    let h_id = model.ideal_gas(i)?.properties(t_eval, 1.0e5)?.enthalpy;
                    u_id += ni * (h_id - R * t_eval);
                }
            }
            Ok(u_id - R * t_eval * t_eval * f.f_t)
        };
    for (label, model, t_c, p_c, n_c, ph_c) in [
        (
            "NIST ternary",
            &m_nist,
            310.0,
            4.0e6,
            &[1.0, 0.6, 0.4][..],
            Phase::Vapor,
        ),
        (
            "Cubic PR liquid",
            &m_cb_pr,
            250.0,
            6.0e6,
            &[0.4, 0.6][..],
            Phase::Liquid,
        ),
        (
            "MBWR32 liquid",
            &m_mbwr32,
            240.0,
            5.0e6,
            &[0.5, 0.5][..],
            Phase::Liquid,
        ),
        (
            "Wong-Sandler vapor",
            &m_ws,
            320.0,
            3.0e6,
            &[0.5, 0.5][..],
            Phase::Vapor,
        ),
    ] {
        let sumn: f64 = n_c.iter().sum();
        let v_fixed = model.volume_tp(t_c, p_c, n_c, ph_c)?;
        let (_cp, cv_analytic) = model.heat_capacities_tp(t_c, p_c, n_c, ph_c)?;
        let dt_cv = 1.0e-5 * t_c;
        let u_plus = internal_energy_from_helmholtz(model, t_c + dt_cv, v_fixed, n_c)?;
        let u_minus = internal_energy_from_helmholtz(model, t_c - dt_cv, v_fixed, n_c)?;
        let cv_fd = (u_plus - u_minus) / (2.0 * dt_cv);
        assert_phys_close(
            &format!("{label} fixed-V FD Cv"),
            cv_analytic,
            cv_fd,
            sumn * R,
            1e-6,
            1e-8,
        );
    }

    // 5. Phase::MinimumGibbs selects the lower-Gibbs phase across all three reference backends
    for model in [&m_nist, &m_cb_srk, &m_mbwr32] {
        let n_len = model.components().len();
        let n_vec: Vec<f64> = if n_len == 3 {
            vec![0.5, 0.3, 0.2]
        } else {
            vec![0.5, 0.5]
        };
        let t = 260.0;
        let p = 3.0e6;
        let g_min = model.gibbs(t, p, &n_vec, Phase::MinimumGibbs)?;
        let g_liq = model.gibbs(t, p, &n_vec, Phase::Liquid);
        let g_vap = model.gibbs(t, p, &n_vec, Phase::Vapor);
        match (g_liq, g_vap) {
            (Ok(gl), Ok(gv)) => {
                assert_phys_close(
                    "MinimumGibbs selects lower G",
                    g_min,
                    gl.min(gv),
                    R * t,
                    1e-11,
                    1e-11,
                );
            }
            (Ok(gl), Err(_)) => assert_phys_close(
                "MinimumGibbs liquid fallback",
                g_min,
                gl,
                R * t,
                1e-11,
                1e-11,
            ),
            (Err(_), Ok(gv)) => assert_phys_close(
                "MinimumGibbs vapor fallback",
                g_min,
                gv,
                R * t,
                1e-11,
                1e-11,
            ),
            (Err(e), Err(_)) => return Err(e),
        }
    }

    // 6. Zero-component and trace-component 1e-100 log cutoff numerical assertions across 3 regimes:
    //    (a) absent (x_1 = 0.0), (b) below cutoff (x_1 = 1e-120), (c) above cutoff (x_1 = 1e-50)
    {
        let t = 320.0;
        let p = 3.0e6;
        let m_bin = ExtCsp::new(&["C1", "C3"])?;
        let n_tern = [0.6, 0.0, 0.4];
        let n_bin = [0.6, 0.4];

        let v_tern = m_nist.volume_tp(t, p, &n_tern, Phase::Vapor)?;
        let v_bin = m_bin.volume_tp(t, p, &n_bin, Phase::Vapor)?;
        assert_phys_close("zero-comp volume", v_tern, v_bin, v_bin, 1e-12, 1e-12);

        let fug_tern = m_nist.ln_fugacity_derivatives_tp(t, p, &n_tern, Phase::Vapor)?;
        let fug_bin = m_bin.ln_fugacity_derivatives_tp(t, p, &n_bin, Phase::Vapor)?;
        assert!(fug_tern.ln_phi.iter().all(|v| v.is_finite()));
        assert!(fug_tern.dn.iter().all(|v| v.is_finite()));
        assert_phys_close(
            "zero-comp ln_phi[0]",
            fug_tern.ln_phi[0],
            fug_bin.ln_phi[0],
            1.0,
            1e-12,
            1e-12,
        );
        assert_phys_close(
            "zero-comp ln_phi[2]",
            fug_tern.ln_phi[2],
            fug_bin.ln_phi[1],
            1.0,
            1e-12,
            1e-12,
        );

        let ideal1 = m_nist.ideal_gas(1)?.properties(t, p)?;
        let exp_s_cap = ideal1.entropy - R * 1.0e-100_f64.ln();
        let exp_g_cap = ideal1.enthalpy - t * exp_s_cap;
        let exp_s_above = ideal1.entropy - R * 1.0e-50_f64.ln();
        let exp_g_above = ideal1.enthalpy - t * exp_s_above;

        let s_bin_tot = m_bin.entropy(t, p, &n_bin, Phase::Vapor)?;
        let g_bin_tot = m_bin.gibbs(t, p, &n_bin, Phase::Vapor)?;
        let s_abs = m_nist.entropy_derivatives(t, p, &n_tern, Phase::Vapor)?;
        let g_abs = m_nist.gibbs_derivatives(t, p, &n_tern, Phase::Vapor)?;
        assert_phys_close(
            "zero-comp exact 0*ln(0) S",
            s_abs.value,
            s_bin_tot,
            R,
            1e-12,
            1e-12,
        );
        assert_phys_close(
            "zero-comp exact 0*ln(0) G",
            g_abs.value,
            g_bin_tot,
            R * t,
            1e-11,
            1e-11,
        );

        let mut s_dn1_below = 0.0;
        let mut g_dn1_below = 0.0;
        let mut s_dn1_above = 0.0;
        let mut g_dn1_above = 0.0;
        for (label, n_vec, exp_s_id, exp_g_id) in [
            ("absent x=0", [0.6, 0.0, 0.4], exp_s_cap, exp_g_cap),
            (
                "below cutoff x=1e-120",
                [0.6, 1.0e-120, 0.4],
                exp_s_cap,
                exp_g_cap,
            ),
            (
                "above cutoff x=1e-50",
                [0.6, 1.0e-50, 0.4],
                exp_s_above,
                exp_g_above,
            ),
        ] {
            let s_res = m_nist.residual_entropy_derivatives(t, p, &n_vec, Phase::Vapor)?;
            let s_tot = m_nist.entropy_derivatives(t, p, &n_vec, Phase::Vapor)?;
            let g_res = m_nist.residual_gibbs_derivatives(t, p, &n_vec, Phase::Vapor)?;
            let g_tot = m_nist.gibbs_derivatives(t, p, &n_vec, Phase::Vapor)?;
            assert!(s_res.dn.iter().all(|v| v.is_finite()));
            assert!(s_tot.dn.iter().all(|v| v.is_finite()));
            assert!(g_res.dn.iter().all(|v| v.is_finite()));
            assert!(g_tot.dn.iter().all(|v| v.is_finite()));
            assert_phys_close(
                &format!("{label} entropy.dn[1] ideal mixing"),
                s_tot.dn[1] - s_res.dn[1],
                exp_s_id,
                R * 250.0,
                1e-12,
                1e-12,
            );
            assert_phys_close(
                &format!("{label} gibbs.dn[1] ideal mixing"),
                g_tot.dn[1] - g_res.dn[1],
                exp_g_id,
                R * t * 250.0,
                1e-12,
                1e-12,
            );
            if label.contains("below") {
                s_dn1_below = s_tot.dn[1];
                g_dn1_below = g_tot.dn[1];
            } else if label.contains("above") {
                s_dn1_above = s_tot.dn[1];
                g_dn1_above = g_tot.dn[1];
            }
        }
        assert!((s_dn1_above - s_dn1_below).abs() > 100.0 * R);
        assert!((g_dn1_above - g_dn1_below).abs() > 100.0 * R * t);

        assert_phys_close(
            "zero-comp total entropy == scalar",
            s_abs.value,
            m_nist.entropy(t, p, &n_tern, Phase::Vapor)?,
            R,
            1e-12,
            1e-12,
        );

        // Pure component embedded in binary [0.0, 1.0]
        let m_pure = ExtCsp::new(&["C3"])?;
        assert_phys_close(
            "embedded pure C3 volume",
            m_bin.volume_tp(t, p, &[0.0, 1.0], Phase::Vapor)?,
            m_pure.volume_tp(t, p, &[1.0], Phase::Vapor)?,
            1e-3,
            1e-12,
            1e-12,
        );
    }

    // 7. Coherent shape alpha setter (replaces unsafe &mut Cubic exposure)
    {
        let mut custom = ExtCsp::new(&["C1", "C3"])?;
        custom.set_shape_alpha(0, AlphaCorrelation::Classic, &[0.52])?;
        assert_phys_close(
            "custom shape alpha",
            custom.shape_eos().alpha(0)?.parameters[0],
            0.52,
            1.0,
            1e-15,
            1e-15,
        );
        verify_helmholtz_and_tp_fd(&custom, 320.0, 2.0e6, &[0.5, 0.5], Phase::Vapor)?;
    }

    Ok(())
}

#[test]
fn test_extcsp_nist_c3_only_constraint() -> Result<()> {
    // Propane aliases accepted for NIST_MEOS
    assert!(
        ExtCsp::with_options(
            &["C1", "C3"],
            ExtCspOptions {
                reference_component: "C3",
                reference_eos: ExtCspReferenceKind::NistMeos,
                ..Default::default()
            },
        )
        .is_ok()
    );
    assert!(
        ExtCsp::with_options(
            &["C1", "C3"],
            ExtCspOptions {
                reference_component: "Propane",
                reference_eos: ExtCspReferenceKind::NistMeos,
                ..Default::default()
            },
        )
        .is_ok()
    );

    // Non-propane NIST fluids rejected
    for bad_ref in ["R134A", "N-H2", "O-H2", "P-H2"] {
        let err = ExtCsp::with_options(
            &["C1", "C3"],
            ExtCspOptions {
                reference_component: bad_ref,
                reference_eos: ExtCspReferenceKind::NistMeos,
                ..Default::default()
            },
        )
        .unwrap_err();
        assert!(
            matches!(err, Error::InvalidInput(_)),
            "expected InvalidInput for NIST_MEOS reference {bad_ref}, got {err:?}"
        );
    }
    Ok(())
}

#[test]
fn test_extcsp_speed_of_sound_rejects_mechanically_unstable_state() -> Result<()> {
    let model = ExtCsp::new(&["C3"])?;
    // Inside the two-phase spinodal for C3: T = 330 K, V = 3.8e-4 m^3/mol has dP/dV > 0.
    let p_tv = model.pressure_tv(330.0, 3.8e-4, &[1.0])?;
    assert!(
        p_tv.dv > 0.0,
        "expected mechanically unstable state (dP/dV > 0), got {}",
        p_tv.dv
    );
    let err = model.speed_of_sound_tv(330.0, 3.8e-4, &[1.0]).unwrap_err();
    assert!(matches!(err, Error::NoPhysicalRoot));
    Ok(())
}

#[test]
fn test_extcsp_pure_reference_fluid_identities() -> Result<()> {
    // When the mixture is the pure reference fluid itself (same shape EOS & alpha),
    // H == N, F == N, and ExtCSP reduces identically to the reference EOS!
    let t = 350.0;
    let p = 2.0e6;
    let n = 2.5;

    // 1. NIST MEOS C3
    let csp_nist = ExtCsp::new(&["C3"])?;
    let meos = Multiparameter::new("C3")?;
    let (h, f) = csp_nist.shape_factors(t, &[n])?;
    let t0 = n * t / f;
    let v = csp_nist.volume_tp(t, p, &[n], Phase::Vapor)?;
    let f_csp = csp_nist.residual_helmholtz_tv(t, v, &[n])?;
    let f_meos = meos.prepare(t0)?.residual_helmholtz(v / h, 1.0)?;
    assert_phys_close(
        "NIST C3 mapped F_res",
        f_csp.f,
        n * f_meos.value,
        n,
        1e-12,
        1e-12,
    );

    // 2. Cubic SRK C3: tc0, pc0 come directly from compdb, so H == N and F == N IDENTICALLY!
    let csp_srk = ExtCsp::with_options(
        &["C3"],
        ExtCspOptions {
            shape_eos: CubicKind::SoaveRedlichKwong,
            reference_eos: ExtCspReferenceKind::Cubic {
                kind: CubicKind::SoaveRedlichKwong,
                alpha: AlphaCorrelation::Classic,
            },
            ..Default::default()
        },
    )?;
    let srk = Cubic::new(CubicKind::SoaveRedlichKwong, &["C3"])?;
    let (h, f) = csp_srk.shape_factors(t, &[n])?;
    assert_phys_close("Pure SRK H == N", h, n, n, 1e-13, 1e-13);
    assert_phys_close("Pure SRK F == N", f, n, n, 1e-13, 1e-13);
    let v_csp = csp_srk.specific_volume(t, p, &[n], Phase::Vapor)?;
    let srk_state = srk.prepare(t, &[1.0])?;
    let v_srk = srk_state.volume(p, Phase::Vapor)?;
    assert_phys_close(
        "Pure SRK volume identity",
        v_csp,
        v_srk,
        v_srk,
        1e-12,
        1e-12,
    );
    assert_phys_close(
        "Pure SRK enthalpy identity",
        csp_srk.enthalpy(t, p, &[1.0], Phase::Vapor)?,
        srk_state.enthalpy(v_srk)?,
        R * t,
        1e-11,
        1e-11,
    );
    assert_phys_close(
        "Pure SRK entropy identity",
        csp_srk.entropy(t, p, &[1.0], Phase::Vapor)?,
        srk_state.entropy(v_srk)?,
        R,
        1e-11,
        1e-11,
    );

    // 3. Cubic PR C3 with genuine Twu alpha ("tcPR", untranslated): H == N and F == N IDENTICALLY!
    let csp_pr_twu = ExtCsp::with_options(
        &["C3"],
        ExtCspOptions {
            shape_eos: CubicKind::PengRobinson,
            shape_alpha: AlphaCorrelation::Twu,
            parameter_reference: "tcPR",
            reference_eos: ExtCspReferenceKind::Cubic {
                kind: CubicKind::PengRobinson,
                alpha: AlphaCorrelation::Twu,
            },
            ..Default::default()
        },
    )?;
    assert_eq!(
        csp_pr_twu.reference_alpha().correlation,
        AlphaCorrelation::Twu
    );
    assert!(!csp_pr_twu.reference_alpha().used_fallback);
    let mut pr_twu = Cubic::with_alpha(
        CubicKind::PengRobinson,
        &["C3"],
        AlphaCorrelation::Twu,
        "tcPR",
    )?;
    pr_twu.disable_volume_translation();
    let (h, f) = csp_pr_twu.shape_factors(t, &[n])?;
    assert_phys_close("Pure PR+Twu H == N", h, n, n, 1e-13, 1e-13);
    assert_phys_close("Pure PR+Twu F == N", f, n, n, 1e-13, 1e-13);
    let v_csp = csp_pr_twu.specific_volume(t, p, &[n], Phase::Vapor)?;
    let v_pr = pr_twu.prepare(t, &[1.0])?.volume(p, Phase::Vapor)?;
    assert_phys_close(
        "Pure PR+Twu volume identity",
        v_csp,
        v_pr,
        v_pr,
        1e-12,
        1e-12,
    );

    // 4. MBWR19 & MBWR32 C3 mapped residual Helmholtz identity
    for kind in [MbwrKind::Bender19, MbwrKind::Mbwr32] {
        let csp_mbwr = ExtCsp::with_options(
            &["C3"],
            ExtCspOptions {
                reference_eos: ExtCspReferenceKind::Mbwr(kind),
                ..Default::default()
            },
        )?;
        let mbwr = Mbwr::new(kind, "C3")?;
        let (h, f) = csp_mbwr.shape_factors(t, &[n])?;
        let t0 = n * t / f;
        let v = csp_mbwr.volume_tp(t, p, &[n], Phase::Vapor)?;
        let f_csp = csp_mbwr.residual_helmholtz_tv(t, v, &[n])?;
        let f_mbwr = mbwr.prepare(t0)?.residual_helmholtz(v / h, 1.0)?;
        assert_phys_close(
            "MBWR C3 mapped F_res",
            f_csp.f,
            n * f_mbwr.value,
            n,
            1e-12,
            1e-12,
        );
    }
    Ok(())
}

#[test]
fn test_extcsp_non_finite_wong_sandler_activity_sum_propagates_error() -> Result<()> {
    let mut ws = thermopack_rs::wong_sandler::WongSandlerData::new(2);
    ws.set_tauij(
        0,
        1,
        thermopack_rs::wong_sandler::Fraction::constant(-1000.0),
    )?;
    ws.set_alphaij(0, 1, 1.0)?;

    let model = ExtCsp::with_options(
        &["CO2", "C1"],
        ExtCspOptions {
            shape_eos: CubicKind::PengRobinson,
            shape_mixing: MixingRule::WongSandler,
            wong_sandler: Some(&ws),
            ..Default::default()
        },
    )?;

    let t = 320.0;
    let p = 3.0e6;
    let n = [1.0, 0.0];
    let v = 7.75e-4;

    assert!(
        model.shape_factors(t, &n).is_err(),
        "shape_factors must propagate error on non-finite activity sum"
    );
    assert!(
        model.volume_tp(t, p, &n, Phase::Vapor).is_err(),
        "volume_tp must propagate error on non-finite activity sum"
    );
    assert!(
        model.volume_tp_derivatives(t, p, &n, Phase::Vapor).is_err(),
        "volume_tp_derivatives must propagate error on non-finite activity sum"
    );
    assert!(
        model.pressure_tv(t, v, &n).is_err(),
        "pressure_tv must propagate error on non-finite activity sum"
    );
    Ok(())
}
