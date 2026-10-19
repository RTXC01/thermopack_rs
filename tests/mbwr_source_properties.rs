//! Source-rounding reproducibility runs in the oracle's process-main context.
//! Windows CRT exp differs by an ulp on new threads; huge extrapolated terms
//! amplify that into pure n-derivative residues. Physical tests run normally.
use thermopack_rs::Mbwr;
fn main() {
    let mut count = 0;
    for row in include_str!("fixtures/mbwr_tp_properties.csv")
        .lines()
        .skip(1)
    {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let model = Mbwr::from_record(fields[1]).unwrap();
        let t = fields[2].parse().unwrap();
        let v = fields[3].parse().unwrap();
        let n = fields[4].parse().unwrap();
        let p = fields[5].split('@').nth(1).unwrap().parse().unwrap();
        let state = model.baseline_residual_tp(t, p, v, n).unwrap();
        let mut actual = vec![state.pressure, state.gibbs_energy];
        for q in [
            state.enthalpy,
            state.entropy,
            state.compressibility,
            state.ln_fugacity_coefficient,
        ] {
            let d = q.derivatives.unwrap();
            actual.extend([q.value, d.dt, d.dp, d.dn]);
        }
        let atol = [
            0.001, 1e-6, 1e-6, 1e-5, 1e-9, 1e-6, 1e-8, 1e-7, 1e-10, 1e-8, 1e-10, 1e-9, 1e-14, 1e-9,
            1e-8, 1e-8, 1e-14, 1e-8,
        ];
        for (i, (&a, b)) in actual.iter().zip(&columns[1..]).enumerate() {
            let b: f64 = b.parse().unwrap();
            assert!(
                (a - b).abs() <= atol[i] + 2e-7 * a.abs().max(b.abs()),
                "{} field {i}: {a} != {b}",
                columns[0]
            );
        }
        count += 1;
    }
    assert_eq!(count, 5983);
    println!("5983 MBWR source property states verified in process-main context");
}
