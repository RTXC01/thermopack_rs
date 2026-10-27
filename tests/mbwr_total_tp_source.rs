use thermopack_rs::{
    EntropyReference, IdealGas, Mbwr, MbwrDensityParameters, MbwrFluid, Phase, R, component,
};

#[path = "../rust/support/mbwr_total_tp_coordinates.rs"]
mod protocol;

const RELATIVE: f64 = 2.0e-7;

fn parse_vector(cell: &str) -> Option<Vec<f64>> {
    if cell.is_empty() || cell.starts_with("ERROR") {
        None
    } else {
        Some(cell.split(',').map(|x| x.parse().unwrap()).collect())
    }
}

fn matches_error(cell: &str, result: Result<Vec<f64>, String>, case: &str, label: &str) {
    if cell.starts_with("ERROR") {
        assert!(
            result.is_err(),
            "{case} {label}: expected {cell}, got {result:?}"
        );
    } else if !cell.is_empty() {
        assert!(
            result.is_ok(),
            "{case} {label}: unexpected error {:?}",
            result.err()
        );
    }
}

fn floors(t: f64, v: f64, n: f64, p: f64, source: bool) -> Vec<f64> {
    let rt = R * t;
    let e = 2.0e-10;
    let volume = [v, 0.0, v / t, v / p, v / n];
    let energy = [n * rt, 0.0, n * R, n * rt / p, rt];
    let entropy = [n * R, 0.0, n * R / t, n * R / p, R];
    let intensive = [1.0, 0.0, 1.0 / t, 1.0 / p, 1.0 / n];
    if source {
        [
            vec![e * n * rt / v, e * v],
            energy
                .iter()
                .chain(&entropy)
                .chain(&intensive)
                .chain(&intensive)
                .map(|x| e * x)
                .collect(),
        ]
        .concat()
    } else {
        let chemical = [rt, 0.0, R, rt / p, rt / n];
        [
            vec![e * n * rt / v],
            volume
                .iter()
                .chain(&energy)
                .chain(&energy)
                .chain(&entropy)
                .chain(&energy)
                .chain(&energy)
                .chain(&chemical)
                .chain(&intensive)
                .chain(&intensive)
                .chain(&[n * R])
                .map(|x| e * x)
                .collect(),
            vec![0.0, e * n * R, 0.0, 2.0e-8],
        ]
        .concat()
    }
}

fn compare_vector(
    case: &str,
    label: &str,
    expected: &[f64],
    actual: &[f64],
    t: f64,
    source: bool,
    unqualified_source: bool,
) {
    let expected_len = if source { 22 } else { 51 };
    assert_eq!(
        expected.len(),
        expected_len,
        "{case} {label} expected length"
    );
    assert_eq!(actual.len(), expected_len, "{case} {label} actual length");
    assert!(
        expected.iter().chain(actual).all(|x| x.is_finite()),
        "{case} {label}: non-finite vector"
    );

    // Total pressure and volume are coordinates 0/1 in both protocols. The
    // source and physical oracle floors are scaled at the expected state.
    let n = case.split('|').nth(3).unwrap().parse::<f64>().unwrap();
    let p = expected[0];
    let v = expected[1];
    let f = floors(t, v, n, p, source);
    assert_eq!(f.len(), expected_len, "{case} {label} floor shape");

    let flag_indexes: Vec<usize> = if source {
        vec![3, 8, 13, 18]
    } else {
        vec![2, 7, 12, 17, 22, 27, 32, 37, 42, 47, 49]
    };
    for i in flag_indexes {
        if unqualified_source && source && [3, 13, 18].contains(&i) {
            for vector in [expected, actual] {
                assert!(
                    vector[i] == 0.0 || vector[i] == 1.0,
                    "{case} {label}: invalid unqualified flag {i}"
                );
                if vector[i] == 0.0 {
                    assert!(
                        vector[i + 1..i + 4].iter().all(|x| *x == 0.0),
                        "{case} {label}: None flag {i} must carry zero derivative coordinates"
                    );
                }
            }
            continue;
        }
        assert!(
            expected[i] == 0.0 || expected[i] == 1.0,
            "{case} {label}: invalid oracle flag {i}"
        );
        assert_eq!(actual[i], expected[i], "{case} {label}: optional flag {i}");
    }

    let ignored: Vec<usize> = if unqualified_source && source {
        (3..7).chain(13..17).chain(18..22).collect()
    } else {
        Vec::new()
    };
    for i in 0..expected_len {
        if ignored.contains(&i) {
            continue;
        }
        assert!(
            (actual[i] - expected[i]).abs()
                <= f[i] + RELATIVE * actual[i].abs().max(expected[i].abs()),
            "{case} {label} coordinate {i}: {} != {} (floor {})",
            actual[i],
            expected[i],
            f[i]
        );
    }
}

fn source_spinodal_unclosed(case: &str, original_density: &[f64]) -> bool {
    let fields: Vec<_> = case.split('|').collect();
    if fields[6] != "SOURCE"
        || original_density.len() != 4
        || !original_density.iter().all(|x| x.is_finite())
        || original_density[0] <= 0.0
    {
        return false;
    }
    let t: f64 = fields[1].parse().unwrap();
    let p: f64 = fields[7].parse().unwrap();
    original_density[3].abs() <= 64.0 * f64::EPSILON * R * t
        && (original_density[2] - p).abs() > 1e-3 + 1e-6 * original_density[2].abs().max(p.abs())
}

fn compare_density(case: &str, expected: &[f64], actual: &[f64]) {
    assert_eq!(expected.len(), 4, "{case} density oracle length");
    assert_eq!(actual.len(), 4, "{case} density native length");
    assert!(
        expected.iter().chain(actual).all(|x| x.is_finite()),
        "{case} density non-finite"
    );
    assert_eq!(actual[1], expected[1], "{case} density phase");
    for (i, floor) in [1e-12, 0.0, 1e-5, 1e-8].into_iter().enumerate() {
        let rel = if i == 1 { 0.0 } else { 2e-9 };
        assert!(
            (actual[i] - expected[i]).abs() <= floor + rel * actual[i].abs().max(expected[i].abs()),
            "{case} density field {i}: {} != {} (floor {floor}, relative {rel})",
            actual[i],
            expected[i]
        );
    }
}

fn main() {
    let mut rows = 0usize;
    let mut cohorts = [0usize; 4];
    let mut unqualified_source_rows = 0usize;
    for (line_no, line) in include_str!("fixtures/mbwr_total_tp.csv")
        .lines()
        .enumerate()
        .skip(1)
    {
        let columns: Vec<_> = line.split(';').collect();
        assert_eq!(columns.len(), 5, "fixture line {}", line_no + 1);
        let case = columns[0];
        let f: Vec<_> = case.split('|').collect();
        assert_eq!(f.len(), 10, "case {case}");
        let t: f64 = f[1].parse().unwrap();
        let v: f64 = f[2].parse().unwrap();
        let n: f64 = f[3].parse().unwrap();
        let p: f64 = f[7].parse().unwrap();
        let phase = match f[8] {
            "1" => Phase::Liquid,
            "2" => Phase::Vapor,
            _ => Phase::MinimumGibbs,
        };
        let mut model = if f[4] == "DEFAULT" {
            MbwrFluid::from_record(f[0]).unwrap()
        } else {
            let c = component(f[4]).unwrap();
            MbwrFluid::with_ideal_gas(
                Mbwr::from_record(f[0]).unwrap(),
                IdealGas::with_integration_constants(c, 0.0, 0.0).unwrap(),
                c.molar_mass * 1.0e-3,
            )
            .unwrap()
        };
        match f[5].parse::<usize>().unwrap() {
            1 | 3 => {
                model
                    .ideal_gas_mut()
                    .set_formation_enthalpy(-12345.678)
                    .unwrap();
                model
                    .ideal_gas_mut()
                    .set_standard_entropy(210.25, EntropyReference::Atmosphere)
                    .unwrap();
            }
            _ => {}
        }
        if f[5] == "2" || f[5] == "3" {
            model.ideal_gas_mut().set_heat_capacity(8, &[3.5]).unwrap();
        }

        let (baseline, physical, classical, density) = match f[6] {
            "SOURCE" => {
                cohorts[2] += 1;
                let a: Vec<f64> = f[9].split(',').map(|x| x.parse().unwrap()).collect();
                let state = model.baseline_caloric_tp_with_parameters(
                    t,
                    p,
                    n,
                    phase,
                    MbwrDensityParameters::new(a[0], a[1], a[2], a[3]).unwrap(),
                );
                match state {
                    Err(e) => (
                        Err(e.to_string()),
                        Err(e.to_string()),
                        Err(e.to_string()),
                        Err(e.to_string()),
                    ),
                    Ok(state) => {
                        let it = state.iteration;
                        let density = Ok(vec![
                            it.density,
                            if it.phase == Phase::Liquid { 1.0 } else { 2.0 },
                            it.pressure,
                            it.pressure_density_derivative,
                        ]);
                        let classical = state
                            .properties
                            .map(|q| protocol::source(&q))
                            .map_err(|e| e.to_string());
                        (Ok(Vec::new()), Ok(Vec::new()), classical, density)
                    }
                }
            }
            "BOUNDED" => {
                cohorts[3] += 1;
                let baseline = model
                    .baseline_caloric_tvp(t, v, n)
                    .map(|q| protocol::physical(&q))
                    .map_err(|e| e.to_string());
                let physical = model
                    .caloric_tp(t, p, n, phase, f[9].parse().unwrap())
                    .map(|q| protocol::physical(&q))
                    .map_err(|e| e.to_string());
                let classical = model
                    .baseline_caloric_tp_at_volume(t, p, v, n)
                    .map(|q| protocol::source(&q))
                    .map_err(|e| e.to_string());
                (baseline, physical, classical, Ok(Vec::new()))
            }
            "INVALID" => {
                cohorts[1] += 1;
                let b = model
                    .baseline_caloric_tvp(t, v, n)
                    .map(|q| protocol::physical(&q))
                    .map_err(|e| e.to_string());
                let ph = model
                    .caloric_tvp(t, v, n)
                    .map(|q| protocol::physical(&q))
                    .map_err(|e| e.to_string());
                let c = model
                    .baseline_caloric_tp_at_volume(t, p, v, n)
                    .map(|q| protocol::source(&q))
                    .map_err(|e| e.to_string());
                (b, ph, c, Ok(Vec::new()))
            }
            _ => {
                cohorts[0] += 1;
                let b = model
                    .baseline_caloric_tvp(t, v, n)
                    .map(|q| protocol::physical(&q))
                    .map_err(|e| e.to_string());
                let ph = model
                    .caloric_tvp(t, v, n)
                    .map(|q| protocol::physical(&q))
                    .map_err(|e| e.to_string());
                let c = model
                    .baseline_caloric_tp_at_volume(t, p, v, n)
                    .map(|q| protocol::source(&q))
                    .map_err(|e| e.to_string());
                (b, ph, c, Ok(Vec::new()))
            }
        };

        let original_density = parse_vector(columns[4]).unwrap_or_default();
        let unqualified_source = source_spinodal_unclosed(case, &original_density);
        if unqualified_source {
            unqualified_source_rows += 1;
        }
        let compare_cell =
            |cell: &str, label: &str, result: Result<Vec<f64>, String>, source: bool| {
                if cell.is_empty() {
                    return;
                }
                if cell.starts_with("ERROR") {
                    matches_error(cell, result, case, label);
                    return;
                }
                let expected = parse_vector(cell).unwrap();
                let actual =
                    result.unwrap_or_else(|e| panic!("{case} {label}: unexpected error: {e}"));
                compare_vector(
                    case,
                    label,
                    &expected,
                    &actual,
                    t,
                    source,
                    unqualified_source && label == "classical",
                );
            };
        compare_cell(columns[1], "baseline", baseline, false);
        compare_cell(columns[2], "physical", physical, false);
        compare_cell(columns[3], "classical", classical, true);
        if !columns[4].is_empty() {
            if columns[4].starts_with("ERROR") {
                matches_error(columns[4], density, case, "density");
            } else {
                let expected = parse_vector(columns[4]).unwrap();
                let actual =
                    density.unwrap_or_else(|e| panic!("{case} density: unexpected error: {e}"));
                compare_density(case, &expected, &actual);
            }
        }

        rows += 1;
    }
    assert_eq!(rows, 9752);
    assert_eq!(cohorts, [4959, 1035, 3434, 324]);
    assert_eq!(
        unqualified_source_rows, 93,
        "source-only conditioning domain derived from original density coordinates"
    );
    assert_eq!(
        rows - unqualified_source_rows,
        9659,
        "fully qualified requests exclude the source conditioning domain"
    );
}
