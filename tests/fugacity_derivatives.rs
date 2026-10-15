use thermopack_rs::{Cubic, CubicKind, Phase};

fn phi(model: &Cubic, t: f64, p: f64, x: &[f64], phase: Phase) -> Vec<f64> {
    let state = model.prepare(t, x).unwrap();
    let v = state.volume(p, phase).unwrap();
    let mut result = vec![0.0; x.len()];
    state.ln_fugacity_into(v, &mut result).unwrap();
    result
}

#[test]
fn analytic_fugacity_derivatives() {
    for kind in [
        CubicKind::PengRobinson,
        CubicKind::SoaveRedlichKwong,
        CubicKind::RedlichKwong,
        CubicKind::VanDerWaals,
    ] {
        let mut model = Cubic::new(kind, &["CO2", "C1", "N2"]).unwrap();
        for lij in [0.0, 0.07] {
            model.set_lij(0, 1, lij).unwrap();
            for (t, p) in [(180.0, 1e5), (300.0, 5e6), (500.0, 1e3)] {
                let x = [0.7, 0.2, 0.1];
                for phase in [Phase::Liquid, Phase::Vapor] {
                    let state = model.prepare(t, &x).unwrap();
                    let v = state.volume(p, phase).unwrap();
                    let (mut dt, mut dp, mut dn) = ([0.0; 3], [0.0; 3], [0.0; 9]);
                    state
                        .fugacity_derivatives_into(v, &mut dt, &mut dp, &mut dn)
                        .unwrap();
                    let ht = 1e-3;
                    let hp = p * 1e-5;
                    let pt = phi(&model, t + ht, p, &x, phase);
                    let mt = phi(&model, t - ht, p, &x, phase);
                    let pp = phi(&model, t, p + hp, &x, phase);
                    let mp = phi(&model, t, p - hp, &x, phase);
                    for i in 0..3 {
                        let tfd = (pt[i] - mt[i]) / (2.0 * ht);
                        let pfd = (pp[i] - mp[i]) / (2.0 * hp);
                        assert!(
                            (dt[i] - tfd).abs() < 1e-7 * tfd.abs().max(1.0),
                            "T {kind:?} {t} {p} {dt:?} {tfd}"
                        );
                        assert!(
                            (dp[i] - pfd).abs() < 1e-10 + 1e-5 * pfd.abs(),
                            "P {kind:?} {t} {p}"
                        );
                    }
                    for j in 0..3 {
                        let h = 1e-5;
                        let mut plus = x;
                        let mut minus = x;
                        plus[j] += h;
                        minus[j] -= h;
                        for xi in &mut plus {
                            *xi /= 1.0 + h;
                        }
                        for xi in &mut minus {
                            *xi /= 1.0 - h;
                        }
                        let pp = phi(&model, t, p, &plus, phase);
                        let mp = phi(&model, t, p, &minus, phase);
                        let mut gibbs_duhem = 0.0;
                        for i in 0..3 {
                            let fd = (pp[i] - mp[i]) / (2.0 * h);
                            assert!(
                                (dn[i * 3 + j] - fd).abs() < 2e-6 * fd.abs().max(1.0),
                                "n {kind:?} {t} {p} {i} {j}: {} {fd}",
                                dn[i * 3 + j]
                            );
                            assert!((dn[i * 3 + j] - dn[j * 3 + i]).abs() < 1e-10);
                            gibbs_duhem += x[i] * dn[i * 3 + j];
                        }
                        assert!(gibbs_duhem.abs() < 1e-9, "Gibbs-Duhem {gibbs_duhem}");
                    }
                }
            }
        }
    }
}
