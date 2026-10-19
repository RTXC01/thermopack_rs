//! Coarse source grid comparisons run in the oracle's process-main context.
use thermopack_rs::{Mbwr, MbwrDensityParameters, mbwr_records};
fn main() {
    let models: Vec<_> = mbwr_records()
        .iter()
        .map(|r| Mbwr::from_record(r.source_name).unwrap())
        .collect();
    let model_for = |key: &str| -> &Mbwr {
        let record = key.split('|').nth(1).unwrap();
        &models[record
            .strip_prefix("mbwrcx")
            .unwrap()
            .parse::<usize>()
            .unwrap()
            - 1]
    };
    let mut statuses = [0; 3];
    for row in include_str!("fixtures/mbwr_initialization_status.csv")
        .lines()
        .skip(1)
    {
        let (case, status) = row.split_once(';').unwrap();
        let status: usize = status.parse().unwrap();
        let model = model_for(case);
        assert_eq!(model.source_initialization().is_ok(), status == 0, "{case}");
        assert_eq!(
            model.source_density_parameters().is_ok(),
            status == 0,
            "{case}"
        );
        statuses[status] += 1;
    }
    assert_eq!(statuses, [32, 16, 17]);
    let mut initialized = 0;
    for row in include_str!("fixtures/mbwr_initialization.csv")
        .lines()
        .skip(1)
    {
        let c: Vec<_> = row.split(';').collect();
        let model = model_for(c[0]);
        let s = model.source_initialization().unwrap();
        let p = s.parameters;
        let actual = [
            p.critical_temperature,
            p.critical_pressure,
            p.critical_density,
            s.compressibility,
            p.acentric_factor,
            s.srk_covolume,
            s.srk_alpha_coefficient,
            s.srk_attraction,
        ];
        let atol = [1e-8, 1e-3, 1e-6, 1e-10, 1e-14, 1e-14, 1e-14, 1e-12];
        for (i, a) in actual.into_iter().enumerate() {
            let b: f64 = c[i + 1].parse().unwrap();
            assert!(
                (a - b).abs() <= atol[i] + 2e-10 * a.abs().max(b.abs()),
                "{} field {i}: {a} != {b}",
                c[0]
            );
        }
        assert_eq!(
            s.from_critical_table,
            model.source_critical_point().is_some()
        );
        initialized += 1;
    }
    let mut grids = 0;
    for row in include_str!("fixtures/mbwr_critical_grid.csv")
        .lines()
        .skip(1)
    {
        let c: Vec<_> = row.split(';').collect();
        let f: Vec<_> = c[0].split('|').collect();
        let a: Vec<f64> = f[5]
            .split('@')
            .nth(1)
            .unwrap()
            .split(',')
            .map(|s| s.parse().unwrap())
            .collect();
        let result = model_for(c[0])
            .baseline_critical_grid(MbwrDensityParameters::new(a[0], a[1], a[2], 0.).unwrap())
            .unwrap();
        let actual = [result.temperature, result.pressure, 1. / result.volume];
        let atol = [1e-8, 1e-3, 1e-6];
        for (i, a) in actual.into_iter().enumerate() {
            let b: f64 = c[i + 1].parse().unwrap();
            assert!(
                (a - b).abs() <= atol[i] + 2e-10 * a.abs().max(b.abs()),
                "{} field {i}: {a} != {b}",
                c[0]
            );
        }
        grids += 1;
    }
    assert_eq!((initialized, grids), (32, 27));
    println!(
        "{initialized} source initializations, {grids} explicit critical grids, and {statuses:?} status classes verified"
    );
}
