//! Original safe-saturation results and independently certified MBWR states.
//! The fixture's `native_observed` column is classification metadata only; it
//! is never used as a numeric expected value.
use thermopack_rs::{Mbwr, MbwrDensityParameters, MbwrKind, R};

const FIXTURE: &str = include_str!("fixtures/mbwr_coexistence.csv");

fn close(a: f64, b: f64, floor: f64, relative: f64, label: &str) {
    assert!(
        (a - b).abs() <= floor + relative * a.abs().max(b.abs()),
        "{label}: {a} != {b} (floor={floor}, relative={relative})"
    );
}

fn read_values(text: &str, row: &str, column: &str) -> Result<Vec<f64>, String> {
    text.split(',')
        .map(|part| {
            part.parse::<f64>()
                .map_err(|e| format!("{row} {column} value {part:?}: {e}"))
        })
        .collect()
}

fn saturation(model: &Mbwr, request: &[&str]) -> thermopack_rs::Result<thermopack_rs::Saturation> {
    let known: f64 = request[2].parse().expect("finite fixture input");
    let max_density: f64 = request[3].parse().expect("finite fixture max density");
    let lower: f64 = request[4].parse().expect("finite fixture lower bound");
    let upper: f64 = request[5].parse().expect("finite fixture upper bound");
    match request[1] {
        "P" => model.saturation_pressure(known, max_density),
        "T" => model.saturation_temperature(known, lower, upper, max_density),
        mode => panic!("unknown saturation request mode {mode:?}"),
    }
}

fn source_pressure_roundoff(model: &Mbwr, t: f64, density: f64) -> f64 {
    // Critical metadata is not used to form the original MBWR pressure
    // polynomial. Valid positive values let us inspect its source-order terms
    // even for coefficient fits without source critical initialization data.
    let dummy = MbwrDensityParameters::new(1., 1., 1., 0.).unwrap();
    let prepared = model.prepare_baseline_density(t, dummy).unwrap();
    let a = prepared.source_parameters();
    let r = density * 1e-3; // the original polynomial uses mol/L
    let (polynomial_count, exponential_count) = match model.kind() {
        MbwrKind::Bender19 => (6, 2),
        MbwrKind::Mbwr32 => (9, 6),
    };
    let mut magnitude = 0.;
    for (n, coefficient) in a.iter().enumerate().take(polynomial_count + 1).skip(1) {
        magnitude += (coefficient * r.powi(n as i32)).abs();
    }
    let gamma = a[0];
    let exponential = (-gamma * r * r).exp();
    for m in 1..=exponential_count {
        let coefficient = a[polynomial_count + m];
        magnitude += (coefficient * r.powi((2 * m + 1) as i32) * exponential).abs();
    }
    64. * f64::EPSILON * magnitude
}

fn reduced_chemical_potential(model: &Mbwr, t: f64, density: f64) -> f64 {
    let alpha = model.residual_density(t, density).unwrap();
    density.ln() + alpha.value + density * alpha.drho
}

fn check_certified_properties(
    model: &Mbwr,
    sat: &thermopack_rs::Saturation,
    max_density: f64,
    label: &str,
) {
    let (t, p, liquid_volume, vapor_volume, dpdt) = (
        sat.temperature,
        sat.pressure,
        sat.liquid_volume,
        sat.vapor_volume,
        sat.dpdt,
    );
    assert!(
        p > 0. && t > 0. && liquid_volume > 0. && vapor_volume > liquid_volume,
        "{label}"
    );
    let rho_l = liquid_volume.recip();
    let rho_v = vapor_volume.recip();
    assert!(
        rho_l < max_density && rho_v < max_density,
        "{label}: root exceeds explicit density bound"
    );

    let pressure_l = model.pressure_density(t, rho_l).unwrap();
    let pressure_v = model.pressure_density(t, rho_v).unwrap();
    let pressure_tol_l = 2e-7 * p + source_pressure_roundoff(model, t, rho_l);
    let pressure_tol_v = 2e-7 * p + source_pressure_roundoff(model, t, rho_v);
    assert!(
        (pressure_l.value - p).abs() <= pressure_tol_l,
        "{label}: liquid pressure {} != {p}; tol {pressure_tol_l}",
        pressure_l.value
    );
    assert!(
        (pressure_v.value - p).abs() <= pressure_tol_v,
        "{label}: vapor pressure {} != {p}; tol {pressure_tol_v}",
        pressure_v.value
    );
    assert!(
        pressure_l.drho > 0. && pressure_v.drho > 0.,
        "{label}: unstable coexistence endpoint"
    );

    let mu_l = reduced_chemical_potential(model, t, rho_l);
    let mu_v = reduced_chemical_potential(model, t, rho_v);
    let mu_tol =
        2e-8 * (rho_l / rho_v).ln() + 128. * f64::EPSILON * mu_l.abs().max(mu_v.abs()).max(1.);
    assert!(
        (mu_l - mu_v).abs() <= mu_tol,
        "{label}: Δμ={} exceeds {mu_tol}",
        mu_l - mu_v
    );

    // Residual enthalpy is evaluated independently at each density. Ideal-gas
    // enthalpy cancels between the two phases at common temperature.
    let alpha_l = model.residual_density(t, rho_l).unwrap();
    let alpha_v = model.residual_density(t, rho_v).unwrap();
    let delta_volume = vapor_volume - liquid_volume;
    let latent_enthalpy = -R * t * t * (alpha_v.dt - alpha_l.dt) + p * delta_volume;
    let clapeyron = latent_enthalpy / (t * delta_volume);
    close(dpdt, clapeyron, 1e-7, 2e-8, &format!("{label}: Clapeyron"));

    let equilibrium = model.prepare_equilibrium(t, max_density).unwrap();
    let stable_roots = equilibrium.densities(p).unwrap();
    assert!(
        stable_roots
            .iter()
            .any(|rho| (rho - rho_l).abs() <= 2e-6 * rho_l),
        "{label}: liquid root missing from stable roots {stable_roots:?}"
    );
    assert!(
        stable_roots
            .iter()
            .any(|rho| (rho - rho_v).abs() <= 2e-6 * rho_v),
        "{label}: vapor root missing from stable roots {stable_roots:?}"
    );
    for rho in stable_roots {
        let mu = reduced_chemical_potential(model, t, rho);
        assert!(
            mu >= mu_l - mu_tol,
            "{label}: another stable root has lower G (ρ={rho}, μ={mu} < {mu_l})"
        );
    }
}

fn main() {
    let mut rows = 0;
    let mut original_gold = 0;
    let mut certified = 0;
    let mut unqualified_errors = 0;
    let mut contexts = [0usize; 3];
    for line in FIXTURE.lines().skip(1) {
        let columns: Vec<_> = line.split(';').collect();
        assert_eq!(columns.len(), 5, "fixture row: {line}");
        let case = columns[0];
        let context = case.split(':').next().unwrap();
        contexts[match context {
            "physical" => 0,
            "public" => 1,
            "internal" => 2,
            _ => panic!("{case}"),
        }] += 1;
        let request: Vec<_> = columns[1].split('|').collect();
        assert_eq!(request.len(), 8, "{case}");
        let model = Mbwr::from_record(request[0]).unwrap_or_else(|e| panic!("{case}: {e}"));
        let solved = saturation(&model, &request);
        let observed_error = columns[3].starts_with("ERROR:");
        if observed_error {
            unqualified_errors += 1;
        } else {
            certified += 1;
            assert!(
                solved.is_ok(),
                "{case}: independently certified state did not solve: {:?}",
                solved.err()
            );
        }

        match solved {
            Err(_) if observed_error => {}
            Err(e) => panic!("{case}: unexpected native error for certified state: {e}"),
            Ok(sat) => {
                let max_density: f64 = request[3].parse().unwrap();
                check_certified_properties(&model, &sat, max_density, case);
                if !columns[2].is_empty() {
                    original_gold += 1;
                    let expected = read_values(columns[2], case, "original solver").unwrap();
                    assert_eq!(expected.len(), 5, "{case}");
                    let actual = [
                        sat.temperature,
                        sat.pressure,
                        sat.liquid_volume,
                        sat.vapor_volume,
                        sat.dpdt,
                    ];
                    let floors = [1e-7, 1e-4, 1e-11, 1e-10, 1e-4];
                    for i in 0..5 {
                        close(
                            actual[i],
                            expected[i],
                            floors[i],
                            2e-6,
                            &format!("{case} coordinate {i}"),
                        );
                    }
                }
            }
        }
        rows += 1;
    }
    assert_eq!(rows, 1008);
    assert_eq!(original_gold, 404);
    assert_eq!(certified, 920);
    assert_eq!(unqualified_errors, 88);
    assert_eq!(contexts, [390, 276, 342]);
    println!(
        "MBWR coexistence: {rows} requests, {original_gold} original solver matches, {certified} certified returns, {unqualified_errors} unqualified error observations"
    );
}
