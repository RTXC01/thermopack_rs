use thermopack_rs::{
    AlphaCorrelation, Cubic, CubicKind, ExtCsp, ExtCspOptions, ExtCspReferenceKind, Mbwr, MbwrKind,
    MixingRule, Multiparameter, Phase, R, Result,
};

fn assert_close(label: &str, a: f64, b: f64, rtol: f64, atol: f64) {
    let err = (a - b).abs();
    let tol = atol + rtol * a.abs().max(b.abs());
    assert!(
        err <= tol,
        "{label}: expected={a:.15e}, actual={b:.15e}, err={err:.3e}, tol={tol:.3e}"
    );
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
        let v = model.specific_volume_derivatives(t, p, &x, phase)?;
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

        assert_eq!(actual.len(), expected.len(), "length mismatch for {label}");
        for (idx, (&a, &b)) in expected.iter().zip(actual.iter()).enumerate() {
            assert_close(&format!("{label}[{idx}]"), a, b, 1e-7, 1e-9);
        }
        rows += 1;
    }
    assert!(rows >= 20, "expected at least 20 ExtCSP fixture rows");
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
    let v = model.volume_tp(t, p, n, phase)?;
    let f = model.residual_helmholtz_tv(t, v, n)?;

    // 1. Central finite differences of all 10 Helmholtz blocks
    let dt = 1e-5 * t;
    let dv = 1e-5 * v;
    let fp_t = model.residual_helmholtz_tv(t + dt, v, n)?;
    let fm_t = model.residual_helmholtz_tv(t - dt, v, n)?;
    let fp_v = model.residual_helmholtz_tv(t, v + dv, n)?;
    let fm_v = model.residual_helmholtz_tv(t, v - dv, n)?;

    assert_close("F_T", f.f_t, (fp_t.f - fm_t.f) / (2.0 * dt), 1e-5, 1e-7);
    assert_close("F_V", f.f_v, (fp_v.f - fm_v.f) / (2.0 * dv), 1e-5, 1e-7);
    assert_close(
        "F_TT",
        f.f_tt,
        (fp_t.f_t - fm_t.f_t) / (2.0 * dt),
        1e-5,
        1e-7,
    );
    assert_close(
        "F_VV",
        f.f_vv,
        (fp_v.f_v - fm_v.f_v) / (2.0 * dv),
        1e-5,
        1e-7,
    );
    assert_close(
        "F_TV",
        f.f_tv,
        (fp_t.f_v - fm_t.f_v) / (2.0 * dt),
        1e-5,
        1e-7,
    );

    let sumn: f64 = n.iter().sum();
    let dn = 1e-6 * sumn;
    for j in 0..nc {
        let mut np = n.to_vec();
        let mut nm = n.to_vec();
        np[j] += dn;
        nm[j] -= dn;
        let fp = model.residual_helmholtz_tv(t, v, &np)?;
        let fm = model.residual_helmholtz_tv(t, v, &nm)?;
        assert_close("F_n", f.f_n[j], (fp.f - fm.f) / (2.0 * dn), 1e-5, 1e-7);
        assert_close(
            "F_Tn",
            f.f_tn[j],
            (fp.f_t - fm.f_t) / (2.0 * dn),
            1e-5,
            1e-7,
        );
        assert_close(
            "F_Vn",
            f.f_vn[j],
            (fp.f_v - fm.f_v) / (2.0 * dn),
            1e-5,
            1e-7,
        );
        for i in 0..nc {
            assert_close(
                "F_nn",
                f.f_nn[i * nc + j],
                (fp.f_n[i] - fm.f_n[i]) / (2.0 * dn),
                1e-5,
                1e-7,
            );
            assert_close(
                "F_nn symmetry",
                f.f_nn[i * nc + j],
                f.f_nn[j * nc + i],
                1e-12,
                1e-12,
            );
        }
    }

    // 2. Euler homogeneity of F(T, V, n): sum(n_i F_n[i]) + V F_V = F
    let euler_f: f64 = n.iter().zip(&f.f_n).map(|(ni, fni)| ni * fni).sum::<f64>() + v * f.f_v;
    assert_close("Euler F", euler_f, f.f, 1e-10, 1e-10);
    for i in 0..nc {
        let euler_fn: f64 = (0..nc).map(|j| n[j] * f.f_nn[i * nc + j]).sum::<f64>() + v * f.f_vn[i];
        assert_close("Euler F_n", euler_fn, 0.0, 1e-10, 1e-10);
    }

    // 3. TP derivatives and thermodynamic identities
    let zd = model.compressibility_factor_derivatives(t, p, n, phase)?;
    let vd = model.specific_volume_derivatives(t, p, n, phase)?;
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

    assert_close("dS/dP = -dV/dT", sd.dp, -vd.dt, 1e-10, 1e-10);
    assert_close(
        "dH/dP = V - T*dV/dT",
        hd.dp,
        vd.value - t * vd.dt,
        1e-10,
        1e-10,
    );
    assert_close("dH/dT = T*dS/dT", hd.dt, t * sd.dt, 1e-10, 1e-10);
    assert_close("Cp = dH/dT", cp, hd.dt, 1e-10, 1e-10);
    assert!(cp > cv && cv > 0.0);
    assert_close("sos_tp == sos_tv", sos_tp, sos_tv, 1e-11, 1e-11);
    assert_close(
        "G = H - T*S",
        gd.value,
        hd.value - t * sd.value,
        1e-10,
        1e-10,
    );
    assert_close("dG/dT = -S", gd.dt, -sd.value, 1e-10, 1e-10);
    assert_close("dG/dP = V", gd.dp, vd.value, 1e-10, 1e-10);
    assert_close("dG_res/dT = -S_res", grd.dt, -srd.value, 1e-10, 1e-10);
    for i in 0..nc {
        assert_close(
            "dG_res/dn = R*T*lnphi",
            grd.dn[i],
            R * t * fug.ln_phi[i],
            1e-10,
            1e-10,
        );
        assert_close(
            "dlnphi/dT = -H_bar_res/(RT^2)",
            fug.dt[i],
            -hrd.dn[i] / (R * t * t),
            1e-10,
            1e-10,
        );
        assert_close(
            "dlnphi/dP = V_bar/(RT) - 1/P",
            fug.dp[i],
            vd.dn[i] / (R * t) - 1.0 / p,
            1e-10,
            1e-10,
        );
        let mut gd_sum = 0.0;
        for j in 0..nc {
            gd_sum += n[j] * fug.dn[i * nc + j];
            assert_close(
                "dlnphi/dn symmetry",
                fug.dn[i * nc + j],
                fug.dn[j * nc + i],
                1e-11,
                1e-11,
            );
        }
        assert_close("Gibbs-Duhem", gd_sum, 0.0, 1e-10, 1e-10);
    }

    // Spot-check TP FD for Z, H, S, lnphi w.r.t. T, P, n
    let dp = 1e-5 * p;
    assert_close(
        "dZ/dT FD",
        zd.dt,
        (model.compressibility_factor(t + dt, p, n, phase)?
            - model.compressibility_factor(t - dt, p, n, phase)?)
            / (2.0 * dt),
        1e-5,
        1e-7,
    );
    assert_close(
        "dZ/dP FD",
        zd.dp,
        (model.compressibility_factor(t, p + dp, n, phase)?
            - model.compressibility_factor(t, p - dp, n, phase)?)
            / (2.0 * dp),
        1e-5,
        1e-7,
    );
    assert_close(
        "dH/dT FD",
        hd.dt,
        (model.enthalpy(t + dt, p, n, phase)? - model.enthalpy(t - dt, p, n, phase)?) / (2.0 * dt),
        1e-5,
        1e-6,
    );
    assert_close(
        "dS/dP FD",
        sd.dp,
        (model.entropy(t, p + dp, n, phase)? - model.entropy(t, p - dp, n, phase)?) / (2.0 * dp),
        1e-5,
        1e-6,
    );
    for j in 0..nc {
        let mut np = n.to_vec();
        let mut nm = n.to_vec();
        np[j] += dn;
        nm[j] -= dn;
        let lnp_p = model.ln_fugacity_tp(t, p, &np, phase)?;
        let lnp_m = model.ln_fugacity_tp(t, p, &nm, phase)?;
        for i in 0..nc {
            assert_close(
                "dlnphi/dn FD",
                fug.dn[i * nc + j],
                (lnp_p[i] - lnp_m[i]) / (2.0 * dn),
                1e-5,
                1e-7,
            );
        }
    }
    Ok(())
}

#[test]
fn test_extcsp_all_backends_mixing_rules_and_twu_alpha() -> Result<()> {
    // 1. NIST MEOS + VdW (classic)
    let m_nist = ExtCsp::new(&["C1", "C2", "C3"])?;
    verify_helmholtz_and_tp_fd(&m_nist, 310.0, 4.0e6, &[1.0, 0.6, 0.4], Phase::Vapor)?;
    verify_helmholtz_and_tp_fd(&m_nist, 230.0, 6.0e6, &[0.5, 0.3, 0.2], Phase::Liquid)?;

    // 2. Cubic SRK & PR reference backends (upstream bug #1 fixed!)
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

    // 3. MBWR19 (Bender) & MBWR32 reference backends (upstream bugs #2 & #3 fixed!)
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

    // 4. Twu alpha correlation (upstream bug #4 fixed!)
    let m_twu = ExtCsp::with_options(
        &["C1", "C3"],
        ExtCspOptions {
            shape_eos: CubicKind::PengRobinson,
            shape_alpha: AlphaCorrelation::Twu,
            ..Default::default()
        },
    )?;
    verify_helmholtz_and_tp_fd(&m_twu, 300.0, 2.5e6, &[0.4, 0.6], Phase::Vapor)?;

    // 5. Wong-Sandler mixing rule (exercises nonzero H_T, H_TT, H_iT!)
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
    verify_helmholtz_and_tp_fd(&m_ws, 320.0, 3.0e6, &[0.5, 0.5], Phase::Vapor)?;
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
    // Note: C3 tc/pc in compdb vs NIST_MEOS record differ by ~0.01 K, so shape factors
    // are ~N (not exact 1.0 unless tc/pc match). Check reference EOS direct identity for Cubic & MBWR,
    // and check NIST_MEOS consistency via exact t0, v0 mapping:
    let t0 = n * t / f;
    let v = csp_nist.volume_tp(t, p, &[n], Phase::Vapor)?;
    let f_csp = csp_nist.residual_helmholtz_tv(t, v, &[n])?;
    let f_meos = meos.prepare(t0)?.residual_helmholtz(v / h, 1.0)?;
    assert_close(
        "NIST C3 mapped F_res",
        f_csp.f,
        n * f_meos.value,
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
    assert_close("Pure SRK H == N", h, n, 1e-13, 1e-13);
    assert_close("Pure SRK F == N", f, n, 1e-13, 1e-13);
    let v_csp = csp_srk.specific_volume(t, p, &[n], Phase::Vapor)?;
    let srk_state = srk.prepare(t, &[1.0])?;
    let v_srk = srk_state.volume(p, Phase::Vapor)?;
    assert_close("Pure SRK volume identity", v_csp, v_srk, 1e-12, 1e-12);
    assert_close(
        "Pure SRK enthalpy identity",
        csp_srk.enthalpy(t, p, &[1.0], Phase::Vapor)?,
        srk_state.enthalpy(v_srk)?,
        1e-11,
        1e-11,
    );
    assert_close(
        "Pure SRK entropy identity",
        csp_srk.entropy(t, p, &[1.0], Phase::Vapor)?,
        srk_state.entropy(v_srk)?,
        1e-11,
        1e-11,
    );

    // 3. Cubic PR C3 with Twu alpha: H == N and F == N IDENTICALLY!
    let csp_pr_twu = ExtCsp::with_options(
        &["C3"],
        ExtCspOptions {
            shape_eos: CubicKind::PengRobinson,
            shape_alpha: AlphaCorrelation::Twu,
            reference_eos: ExtCspReferenceKind::Cubic {
                kind: CubicKind::PengRobinson,
                alpha: AlphaCorrelation::Twu,
            },
            ..Default::default()
        },
    )?;
    let pr_twu = Cubic::with_alpha(
        CubicKind::PengRobinson,
        &["C3"],
        AlphaCorrelation::Twu,
        "DEFAULT",
    )?;
    let (h, f) = csp_pr_twu.shape_factors(t, &[n])?;
    assert_close("Pure PR+Twu H == N", h, n, 1e-13, 1e-13);
    assert_close("Pure PR+Twu F == N", f, n, 1e-13, 1e-13);
    let v_csp = csp_pr_twu.specific_volume(t, p, &[n], Phase::Vapor)?;
    let v_pr = pr_twu.prepare(t, &[1.0])?.volume(p, Phase::Vapor)?;
    assert_close("Pure PR+Twu volume identity", v_csp, v_pr, 1e-12, 1e-12);

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
        assert_close(
            "MBWR C3 mapped F_res",
            f_csp.f,
            n * f_mbwr.value,
            1e-12,
            1e-12,
        );
    }
    Ok(())
}
