use thermopack_rs::{PcSaft, PcSaftKind};

#[test]
fn mole_derivatives_match_extensive_energy_and_euler_relation() {
    for kind in [PcSaftKind::PcSaft, PcSaftKind::SPcSaft] {
        let model = PcSaft::new_with_options(kind, &["CO2", "C1", "C3"], "Default").unwrap();
        for n in [[0.2, 0.3, 0.5], [0.0, 0.4, 0.6], [0.97, 0.02, 0.01]] {
            for rho in [200.0, 5000.0, 12000.0] {
                let temperature = 320.0;
                let v = 1.0 / rho;
                let prep = model.prepare(temperature, &n).unwrap();
                let z = prep.pressure(rho) / (rho * thermopack_rs::pc_saft::RGAS * temperature);
                let phi = prep.fugacity_coefficients(rho, z);
                let energy = |amounts: [f64; 3]| {
                    let total: f64 = amounts.iter().sum();
                    total
                        * model
                            .prepare(temperature, &amounts)
                            .unwrap()
                            .residual_helmholtz_energy(total / v)
                };
                for i in 0..3 {
                    let h = 1e-5;
                    let mut plus = n;
                    plus[i] += h;
                    let numeric = if n[i] == 0.0 {
                        let mut plus2 = n;
                        plus2[i] += 2.0 * h;
                        (-3.0 * energy(n) + 4.0 * energy(plus) - energy(plus2)) / (2.0 * h)
                    } else {
                        let mut minus = n;
                        minus[i] -= h;
                        (energy(plus) - energy(minus)) / (2.0 * h)
                    };
                    let analytic = phi[i] + z.ln();
                    assert!(
                        (numeric - analytic).abs() < 3e-6 * (1.0 + numeric.abs()),
                        "{kind:?} {n:?} {rho} {i}: {numeric} {analytic}"
                    );
                }
                let weighted: f64 = n.iter().zip(&phi).map(|(x, p)| x * p).sum();
                let expected = prep.residual_helmholtz_energy(rho) + z - 1.0 - z.ln();
                assert!((weighted - expected).abs() < 1e-10 * (1.0 + expected.abs()));
                let scaled = model
                    .prepare(temperature, &[n[0] * 7.0, n[1] * 7.0, n[2] * 7.0])
                    .unwrap()
                    .fugacity_coefficients(rho, z);
                for i in 0..3 {
                    assert!((scaled[i] - phi[i]).abs() < 1e-11 * (1.0 + phi[i].abs()));
                }
            }
        }
    }
}

#[test]
fn mole_derivatives_match_exact_fortran_source() {
    let mut rows = include_str!("fixtures/pcsaft_mole_derivatives.csv").lines();
    assert_eq!(
        rows.next().unwrap(),
        "components;temperature;volume;moles;f_n;simplified;reference"
    );
    let mut count = 0;
    let mut unavailable = Vec::new();
    for row in rows {
        let fields: Vec<_> = row.split(';').collect();
        let components: Vec<_> = fields[0].split(',').collect();
        let temperature: f64 = fields[1].parse().unwrap();
        let volume: f64 = fields[2].parse().unwrap();
        let n: Vec<f64> = fields[3].split(',').map(|x| x.parse().unwrap()).collect();
        let source: Vec<f64> = fields[4].split(',').map(|x| x.parse().unwrap()).collect();
        let kind = if fields[5] == "true" {
            PcSaftKind::SPcSaft
        } else {
            PcSaftKind::PcSaft
        };
        let model = match PcSaft::new_nonassoc(kind, &components, fields[6]) {
            Ok(model) => model,
            Err(error) => {
                // The fit exists, but the public constructor also requires an
                // ideal heat-capacity correlation, which BUTANAL lacks.
                assert_eq!(
                    components[0], "BUTANAL",
                    "unexpected unavailable row {row}: {error:?}"
                );
                assert!(
                    matches!(
                        error,
                        thermopack_rs::Error::Unsupported(
                            "component has no default ideal heat-capacity correlation"
                        )
                    ),
                    "unexpected BUTANAL error for {row}: {error:?}"
                );
                unavailable.push(format!("{row}: {error:?}"));
                continue;
            }
        };
        let prep = model.prepare(temperature, &n).unwrap();
        let rho = n.iter().sum::<f64>() / volume;
        // For mixture states z=1 makes ln(phi_i) exactly F_n,i; this is a
        // test adapter for the residual kernel even when the physical Z <= 0.
        let gradients = prep.fugacity_coefficients(rho, 1.0);
        for i in 0..n.len() {
            let actual = gradients[i];
            assert!(
                actual.is_finite() && source[i].is_finite(),
                "nonfinite {row}: component {i}"
            );
            assert!(
                (actual - source[i]).abs() < 2e-11 * (1.0 + source[i].abs()),
                "{row}: component {i}: source {} actual {actual}",
                source[i]
            );
        }
        count += 1;
    }
    assert_eq!(unavailable.len(), 4, "{unavailable:?}");
    assert_eq!(count, 340);
}

#[test]
fn gibbs_duhem_holds_along_a_fixed_pressure_composition_change() {
    for kind in [PcSaftKind::PcSaft, PcSaftKind::SPcSaft] {
        let model = PcSaft::new_with_options(kind, &["CO2", "C1", "C3"], "Default").unwrap();
        let temperature = 350.0;
        let pressure = 1.0e6;
        let x = [0.2, 0.3, 0.5];
        let h = 1e-5;
        let evaluate = |x: [f64; 3]| {
            let prepared = model.prepare(temperature, &x).unwrap();
            let volume = prepared
                .volume(pressure, thermopack_rs::Phase::Vapor)
                .unwrap();
            let rho = 1.0 / volume;
            prepared.fugacity_coefficients(rho, prepared.compressibility_factor(rho))
        };
        let plus = evaluate([x[0] + h, x[1] - h, x[2]]);
        let minus = evaluate([x[0] - h, x[1] + h, x[2]]);
        let sum: f64 = (0..3)
            .map(|i| x[i] * (plus[i] - minus[i]) / (2.0 * h))
            .sum();
        assert!(sum.abs() < 2e-7, "{kind:?}: {sum}");
    }
}
