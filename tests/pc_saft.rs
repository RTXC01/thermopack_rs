use thermopack_rs::{PcSaft, PcSaftKind, Phase};

#[test]
fn test_pc_saft_pure_propane_states() {
    let eos = PcSaft::new(&["C3"]).expect("failed to load C3");
    assert_eq!(eos.num_components(), 1);

    // State 1: T = 250.0 K, rho = 50.0 mol/m^3 (vapor)
    let prep1 = eos.prepare(250.0, &[1.0]).unwrap();
    let st1 = prep1.state(50.0).unwrap();
    assert!((st1.pressure - 101295.24000317785).abs() < 1e-4);
    assert!((st1.dp_dt - 423.25629840601033).abs() < 1e-4);
    assert!((st1.dp_dv - (-4933674.68354463)).abs() < 1e-1);
    assert!((st1.ln_phi[0] - (-0.025097978831742467)).abs() < 1e-6);
    assert!((st1.h_total - (-108145.27187732508)).abs() < 1e-2);
    assert!((st1.s_total - 257.7543893472314).abs() < 1e-3);
    assert!((st1.sound_speed - 228.11416867657903).abs() < 1e-2);

    // State 2: T = 300.0 K, rho = 665.258 mol/m^3 (dense gas)
    let prep2 = eos.prepare(300.0, &[1.0]).unwrap();
    let st2 = prep2.state(665.258).unwrap();

    assert!((st2.pressure - 1277163.7668680672).abs() < 1e-4);
    assert!((st2.dp_dt - 6528.101604953861).abs() < 1e-4);
    assert!((st2.dp_dv - (-612845684.0165149)).abs() < 1e-1);
    assert!((st2.ln_phi[0] - (-0.20689553787503584)).abs() < 1e-6);
    assert!((st2.h_total - (-106200.71695850792)).abs() < 1e-2);
    assert!((st2.s_total - 245.82041757608926).abs() < 1e-3);
    assert!((st2.sound_speed - 202.85317377824083).abs() < 1e-2);

    // State 3: T = 300.0 K, rho = 11060.948308682653 mol/m^3 (liquid at 5 bar)
    let st3 = prep2.state(11060.948308682653).unwrap();
    assert!((st3.pressure - 500000.0).abs() < 1e-2);
    assert!((st3.dp_dt - 486868.6587188868).abs() < 1e-1);
    assert!((st3.h_total - (-120697.60633388176)).abs() < 1e-1);
    assert!((st3.s_total - 199.28704195871032).abs() < 1e-3);
    assert!((st3.sound_speed - 673.4579332739784).abs() < 1e-1);
}

#[test]
fn test_pc_saft_pure_propane_volumes() {
    let eos = PcSaft::new(&["C3"]).expect("failed to load C3");
    let prep = eos.prepare(300.0, &[1.0]).unwrap();

    let vl = prep.volume(500000.0, Phase::Liquid).unwrap();
    let vv = prep.volume(500000.0, Phase::Vapor).unwrap();

    assert!((vl - 9.040816140646977e-05).abs() < 1e-10);
    assert!((vv - 0.004595162429259163).abs() < 1e-8);

    // Verify consistency: pressure at recovered volumes equals target pressure
    let pl = prep.pressure(1.0 / vl);
    let pv = prep.pressure(1.0 / vv);
    assert!((pl - 500000.0).abs() < 1e-3);
    assert!((pv - 500000.0).abs() < 1e-3);
}

#[test]
fn test_pc_saft_pure_propane_saturation() {
    let eos = PcSaft::new(&["C3"]).expect("failed to load C3");
    let prep = eos.prepare(300.0, &[1.0]).unwrap();

    let sat = prep.saturation_pressure(300.0).unwrap();
    assert!((sat.pressure - 998660.895732584).abs() < 1e-1);
    assert!((sat.liquid_volume - 9.008805107115614e-05).abs() < 1e-9);
    assert!((sat.vapor_volume - 0.002072486771700257).abs() < 1e-7);
    assert!(sat.enthalpy_vaporization > 0.0);
    assert!(sat.entropy_vaporization > 0.0);
}

#[test]
fn test_pc_saft_vs_spc_saft_pure_fluid_identity() {
    let orig = PcSaft::new_with_options(PcSaftKind::PcSaft, &["C3"], "Default").unwrap();
    let simp = PcSaft::new_with_options(PcSaftKind::SPcSaft, &["C3"], "Default").unwrap();

    let prep_orig = orig.prepare(320.0, &[1.0]).unwrap();
    let prep_simp = simp.prepare(320.0, &[1.0]).unwrap();

    for &rho in &[10.0, 100.0, 1000.0, 8000.0, 11000.0] {
        let a_orig = prep_orig.residual_helmholtz_energy(rho);
        let a_simp = prep_simp.residual_helmholtz_energy(rho);
        let p_orig = prep_orig.pressure(rho);
        let p_simp = prep_simp.pressure(rho);

        assert!((a_orig - a_simp).abs() < 1e-14);
        assert!((p_orig - p_simp).abs() < 1e-10);
    }
}

#[test]
fn test_pc_saft_binary_co2_c1() {
    let eos = PcSaft::new(&["CO2", "C1"]).expect("failed to load CO2,C1");
    assert_eq!(eos.num_components(), 2);
    assert!((eos.get_kij(0, 1) - 0.065).abs() < 1e-6);

    let prep = eos.prepare(280.0, &[0.3, 0.7]).unwrap();
    let v = 0.002;
    let rho = 1.0 / v;
    let st = prep.state(rho).unwrap();

    assert!((st.pressure - 1122950.9879415003).abs() < 1.0);
    assert!((st.dp_dt - 4328.035670293376).abs() < 1e-1);
    println!("CO2,C1 dp_dv = {}, expected = -541467953.0951648", st.dp_dv);
    assert!((st.dp_dv - (-541467953.0951648)).abs() < 1e2);
    assert!((st.ln_phi[0] - (-0.059862315254152296)).abs() < 1e-4);
    assert!((st.ln_phi[1] - (-0.02449780863596171)).abs() < 1e-4);
}
