// Line-oriented reference protocol: KIND|COMP,COMP|T|P|x,x|L or V
use std::io::{self, BufRead};
use thermopack_rs::{Cubic, CubicKind, Phase};

fn translated_model(
    kind: CubicKind,
    components: &str,
    config: &[&str],
) -> thermopack_rs::Result<Cubic> {
    let mut model = if config.contains(&"qhelper") {
        Cubic::quantum(&components.split(',').collect::<Vec<_>>())?
    } else {
        Cubic::with_options(
            kind,
            &components.split(',').collect::<Vec<_>>(),
            thermopack_rs::CubicOptions {
                alpha: Some(thermopack_rs::AlphaCorrelation::from_name(config[1])?),
                parameter_reference: config[2],
                volume_translation: config[3] == "1",
            },
        )?
    };
    for operation in &config[4..] {
        let data: Vec<_> = operation.split(':').collect();
        match data[0] {
            "init" => model.init_peneloux_volume_translation(data[1])?,
            "disable" => model.disable_volume_translation(),
            "bounds" => {} // Applied to native caloric options below.
            "qhelper" => {}
            "rr" => {}      // RR iteration inputs, handled after model construction.
            "quality" => {} // Explicit pure two-phase acoustic state.
            "roots" => {}   // Optional fluid acoustic phase-root overrides.
            "amount" => {}  // Total moles in an extensive TV acoustic input.
            "guess" => {}   // UV/SV temperature and pressure guesses.
            "kij" => model.set_kij(
                data[1].parse().unwrap(),
                data[2].parse().unwrap(),
                data[3].parse().unwrap(),
            )?,
            "lij" => model.set_lij(
                data[1].parse().unwrap(),
                data[2].parse().unwrap(),
                data[3].parse().unwrap(),
            )?,
            "beta" => {
                let params = data[3]
                    .split(',')
                    .map(str::parse)
                    .collect::<Result<Vec<f64>, _>>()
                    .map_err(|_| thermopack_rs::Error::InvalidInput("beta parameters"))?;
                model.set_beta(
                    data[1].parse().unwrap(),
                    thermopack_rs::Beta::new(
                        thermopack_rs::BetaCorrelation::from_name(data[2])?,
                        params.try_into().map_err(|_| {
                            thermopack_rs::Error::InvalidInput("two beta parameters")
                        })?,
                    )?,
                )?;
            }
            "set" => {
                let coefficients = data[3]
                    .split(',')
                    .map(str::parse)
                    .collect::<Result<Vec<f64>, _>>()
                    .map_err(|_| thermopack_rs::Error::InvalidInput("translation coefficients"))?;
                model.set_volume_translation(
                    data[1]
                        .parse()
                        .map_err(|_| thermopack_rs::Error::InvalidInput("component index"))?,
                    thermopack_rs::VolumeTranslation::new(
                        thermopack_rs::TranslationKind::from_code(data[2].parse().unwrap())?,
                        coefficients.try_into().map_err(|_| {
                            thermopack_rs::Error::InvalidInput("six coefficients required")
                        })?,
                    )?,
                )?;
            }
            _ => return Err(thermopack_rs::Error::InvalidInput("translation operation")),
        }
    }
    Ok(model)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for line in io::stdin().lock().lines() {
        let line = line?;
        let fields: Vec<_> = line.split('|').collect();
        if fields.len() != 6 {
            return Err("expected six fields".into());
        }
        if fields[0] == "PETS" {
            let parameters: Vec<f64> = fields[1]
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?;
            let model = thermopack_rs::Pets::with_parameters(parameters[0], parameters[1])?;
            let result = (|| -> thermopack_rs::Result<Vec<f64>> {
                if fields[5] == "PETS_POTENTIAL" {
                    let r: Vec<f64> = fields[2].split(',').map(|v| v.parse().unwrap()).collect();
                    let mut values = vec![0.; r.len()];
                    let mut forces = vec![0.; r.len()];
                    model.potential_into(&r, &mut values, Some(&mut forces))?;
                    values.extend(forces);
                    return Ok(values);
                }
                let t = fields[2].parse().unwrap();
                let v = fields[3].parse().unwrap();
                let n = fields[4].parse().unwrap();
                let prepared = model.prepare(t)?;
                let diameter = prepared.diameter();
                let a = prepared.contributions(n / v)?;
                let f = prepared.residual_helmholtz(v, n)?;
                let c = model.lj_truncation_correction(t, v, n)?;
                let p = prepared.caloric_tv(v, n)?;
                let mut values = vec![diameter.value, diameter.dt, diameter.dtt];
                for a in [a.total, a.hard_sphere, a.dispersion] {
                    values.extend([a.value, a.dt, a.drho, a.dtt, a.dtrho, a.drhorho]);
                }
                for f in [f, c] {
                    values.extend([
                        f.value, f.dt, f.dv, f.dn, f.dtt, f.dtv, f.dtn, f.dvv, f.dvn, f.dnn,
                    ]);
                }
                values.extend([
                    p.pressure.value,
                    p.pressure.dt,
                    p.pressure.dv,
                    p.pressure.dn,
                    p.internal_energy,
                    p.entropy,
                    p.helmholtz_energy,
                    p.enthalpy,
                    p.heat_capacity_v,
                    p.chemical_potential,
                ]);
                let (b, c) = model.baseline_virial_coefficients(t)?;
                values.extend([
                    prepared.baseline_speed_of_sound_tv(v, n)?,
                    b,
                    c,
                    p.ln_fugacity_coefficient.unwrap_or(0.),
                    p.heat_capacity_p.unwrap_or(f64::NAN),
                ]);
                Ok(values)
            })();
            match result {
                Ok(values) => println!(
                    "{}",
                    values
                        .iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if fields[5] == "SYMMETRIC" {
            let mut a: Vec<f64> = fields[1]
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?;
            let mut b: Vec<f64> = fields[2]
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?;
            let mut pivots = vec![0; b.len()];
            let mut work = vec![0.; b.len()];
            match thermopack_rs::numerics::solve_symmetric_in_place(
                &mut a,
                &mut b,
                &mut pivots,
                &mut work,
            ) {
                Ok(success) => {
                    print!("{}", usize::from(success));
                    for value in pivots {
                        print!(",{value}");
                    }
                    for value in a.iter().chain(&b) {
                        print!(",{value:.17e}");
                    }
                    println!();
                }
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if fields[5] == "RR" {
            let z: Vec<f64> = fields[1]
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?;
            let k: Vec<f64> = fields[2]
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?;
            let mut state = thermopack_rs::numerics::RachfordRiceState::new(z.len());
            state.vapor_fraction = fields[3].parse()?;
            match thermopack_rs::numerics::rachford_rice_extended(
                &z,
                &k,
                &mut state,
                fields[4] == "1",
            ) {
                Ok(false) => println!("0"),
                Ok(true) => {
                    print!(
                        "1,{:.17e},{:.17e}",
                        state.vapor_fraction, state.liquid_fraction
                    );
                    for value in state.liquid.iter().chain(&state.vapor) {
                        print!(",{value:.17e}");
                    }
                    println!();
                }
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if fields[5] == "IDEAL" || fields[5] == "IDEAL_CONFIG" {
            let t: f64 = fields[2].parse()?;
            let p: f64 = fields[3].parse()?;
            let mut gas = thermopack_rs::IdealGas::new(thermopack_rs::component(fields[1])?)?;
            if fields[5] == "IDEAL_CONFIG" {
                for action in fields[4].split('/') {
                    let tokens: Vec<_> = action.split(':').collect();
                    match tokens[0] {
                        "cp" => gas.set_heat_capacity(
                            tokens[1].parse()?,
                            &tokens[2]
                                .split(',')
                                .filter(|v| !v.is_empty())
                                .map(str::parse)
                                .collect::<Result<Vec<_>, _>>()?,
                        )?,
                        "h" => gas.set_formation_enthalpy(tokens[1].parse()?)?,
                        "s" => gas.set_standard_entropy(
                            tokens[1].parse()?,
                            if tokens[2] == "101325.0" {
                                thermopack_rs::EntropyReference::Atmosphere
                            } else {
                                thermopack_rs::EntropyReference::Bar
                            },
                        )?,
                        _ => return Err("unknown ideal action".into()),
                    }
                }
            }
            match gas.properties(t, p) {
                Ok(properties) => {
                    print!(
                        "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                        properties.enthalpy,
                        properties.entropy,
                        properties.cp,
                        properties.cp / t,
                        -thermopack_rs::R / p
                    );
                    if fields[5] == "IDEAL_CONFIG" {
                        let (corr, params) = gas.heat_capacity_correlation();
                        print!(
                            ",{:.17e},{:.17e},{corr}",
                            gas.formation_enthalpy(),
                            gas.standard_entropy()
                        );
                        for value in params {
                            print!(",{value:.17e}");
                        }
                    }
                    println!();
                }
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        let kind = match fields[0] {
            "PR" => CubicKind::PengRobinson,
            "SRK" => CubicKind::SoaveRedlichKwong,
            "RK" => CubicKind::RedlichKwong,
            "VDW" => CubicKind::VanDerWaals,
            _ => return Err("unknown model".into()),
        };
        if fields[5] == "ALPHA" {
            let result = (|| -> Result<Vec<f64>, thermopack_rs::Error> {
                let selection: Vec<_> = fields[4].split('@').collect();
                let mut model = Cubic::with_alpha(
                    kind,
                    &fields[1].split(',').collect::<Vec<_>>(),
                    thermopack_rs::AlphaCorrelation::from_name(fields[3])?,
                    selection[0],
                )?;
                for setter in &selection[1..] {
                    let data: Vec<_> = setter.split(':').collect();
                    let parameters = data[2]
                        .split(',')
                        .map(str::parse)
                        .collect::<Result<Vec<f64>, _>>()
                        .map_err(|_| {
                            thermopack_rs::Error::InvalidInput("invalid alpha parameter")
                        })?;
                    model.set_alpha(
                        data[0].parse().map_err(|_| {
                            thermopack_rs::Error::InvalidInput("invalid component index")
                        })?,
                        thermopack_rs::AlphaCorrelation::from_name(data[1])?,
                        &parameters,
                    )?;
                }
                let mut values = Vec::new();
                for (i, component) in model.components().iter().enumerate() {
                    let alpha = model.alpha(i)?;
                    let value = alpha
                        .evaluate(fields[2].parse().unwrap(), component.critical_temperature)?;
                    values.extend([alpha.correlation as i32 as f64]);
                    values.extend(alpha.parameters);
                    values.extend([value.value, value.dt, value.dtt]);
                }
                Ok(values)
            })();
            match result {
                Ok(values) => println!(
                    "{}",
                    values
                        .iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        let property_config: Vec<_> = fields[5].split('@').collect();
        let alpha_properties = property_config[0].starts_with("A_");
        let translated = property_config[0].starts_with("V_");
        let mut model = if translated {
            match translated_model(kind, fields[1], &property_config) {
                Ok(model) => model,
                Err(error) => {
                    println!("ERROR: {error}");
                    continue;
                }
            }
        } else if alpha_properties {
            Cubic::with_alpha(
                kind,
                &fields[1].split(',').collect::<Vec<_>>(),
                thermopack_rs::AlphaCorrelation::from_name(property_config[1])?,
                property_config[2],
            )?
        } else {
            Cubic::new(kind, &fields[1].split(',').collect::<Vec<_>>())?
        };
        if alpha_properties {
            for setter in &property_config[3..] {
                let data: Vec<_> = setter.split(':').collect();
                let parameters = data[2]
                    .split(',')
                    .map(str::parse)
                    .collect::<Result<Vec<f64>, _>>()?;
                model.set_alpha(
                    data[0].parse()?,
                    thermopack_rs::AlphaCorrelation::from_name(data[1])?,
                    &parameters,
                )?;
            }
        }
        let t: f64 = fields[2].parse()?;
        if property_config[0] == "V_BETA" {
            let mut values = Vec::new();
            for (i, component) in model.components().iter().enumerate() {
                let beta = model.beta(i)?;
                let value = beta.evaluate(t, component.critical_temperature)?;
                values.push(beta.correlation as i32 as f64);
                values.extend(
                    if beta.correlation == thermopack_rs::BetaCorrelation::Classic {
                        [0.; 2]
                    } else {
                        beta.parameters
                    },
                );
                values.extend([
                    value.value,
                    value.dt,
                    value.dtt,
                    component.critical_temperature,
                    component.critical_pressure,
                    component.molar_mass,
                ]);
            }
            println!(
                "{}",
                values
                    .iter()
                    .map(|v| format!("{v:.17e}"))
                    .collect::<Vec<_>>()
                    .join(",")
            );
            continue;
        }
        if property_config[0] == "V_CI" {
            let mut values = Vec::new();
            for i in 0..model.components().len() {
                let ci = model.volume_translation(i)?;
                let c = ci.evaluate(t)?;
                values.extend([
                    if model.volume_translation_enabled() {
                        1.
                    } else {
                        0.
                    },
                    ci.kind as u8 as f64,
                ]);
                values.extend(ci.coefficients);
                values.extend([
                    c.value,
                    c.dt,
                    c.dtt,
                    if ci.kind == thermopack_rs::TranslationKind::Constant {
                        0.
                    } else {
                        1.
                    },
                ]);
            }
            println!(
                "{}",
                values
                    .iter()
                    .map(|v| format!("{v:.17e}"))
                    .collect::<Vec<_>>()
                    .join(",")
            );
            continue;
        }
        if fields[3] == "PSEUDO" {
            let x: Vec<f64> = fields[4]
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?;
            let critical = model.pseudo_critical(&x)?;
            println!(
                "{:.17e},{:.17e},{:.17e},{:.17e}",
                critical.temperature, critical.pressure, critical.volume, critical.compressibility
            );
            continue;
        }
        if property_config[0] == "V_SAT" {
            match model.saturation_pressure(t) {
                Ok(sat) => println!(
                    "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                    sat.temperature, sat.pressure, sat.liquid_volume, sat.vapor_volume, sat.dpdt
                ),
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if fields[3] == "SAT_T" {
            let sat = model.saturation_temperature(t)?;
            println!(
                "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                sat.temperature, sat.pressure, sat.liquid_volume, sat.vapor_volume, sat.dpdt
            );
            continue;
        }
        if fields[3] == "SAT" {
            let sat = model.saturation_pressure(t)?;
            println!(
                "{:.17e},{:.17e},{:.17e},{:.17e}",
                sat.pressure, sat.liquid_volume, sat.vapor_volume, sat.dpdt
            );
            continue;
        }
        let p: f64 = fields[3].parse()?;
        let x: Vec<f64> = fields[4]
            .split(',')
            .map(str::parse)
            .collect::<Result<_, _>>()?;
        if property_config[0] == "V_CALORIC_TV" {
            match model.prepare(t, &x).and_then(|state| state.caloric_tv(p)) {
                Ok(value) => println!(
                    "{:.17e},{:.17e},{:.17e}",
                    value.pressure, value.internal_energy, value.entropy
                ),
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if matches!(property_config[0], "V_UV" | "V_SV") {
            let mut options = model.volume_flash_options();
            if let Some(guess) = property_config.iter().find(|op| op.starts_with("guess:")) {
                let values: Vec<_> = guess.split(':').collect();
                options.temperature_guess = Some(values[1].parse()?);
                options.pressure_guess = Some(values[2].parse()?);
            }
            let result = if property_config[0] == "V_UV" {
                model.uv_flash_with_options(&x, t, p, options)
            } else {
                model.sv_flash_with_options(&x, t, p, options)
            };
            match result {
                Ok(flash) => {
                    let two = flash.phase == thermopack_rs::FlashPhase::TwoPhase;
                    let mut values = vec![
                        flash.temperature,
                        flash.pressure,
                        f64::from(two),
                        if two { flash.vapor_fraction } else { 0. },
                        if two { flash.liquid_fraction } else { 1. },
                    ];
                    values.extend(if two { flash.liquid } else { x.clone() });
                    values.extend(if two { flash.vapor } else { vec![0.; x.len()] });
                    println!(
                        "{}",
                        values
                            .iter()
                            .map(|v| format!("{v:.17e}"))
                            .collect::<Vec<_>>()
                            .join(",")
                    );
                }
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if property_config[0] == "V_SOUND_TV" {
            let amount: f64 = property_config
                .iter()
                .find(|op| op.starts_with("amount:"))
                .map(|op| op.split(':').nth(1).unwrap().parse().unwrap())
                .unwrap_or(1.);
            let moles: Vec<_> = x.iter().map(|z| z * amount).collect();
            match model.baseline_speed_of_sound_tv(t, p * amount, &moles) {
                Ok(value) => println!("{value:.17e}"),
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if property_config[0] == "V_SOUND" {
            let flash = if let Some(quality) =
                property_config.iter().find(|op| op.starts_with("quality:"))
            {
                let bv: f64 = quality.split(':').nth(1).unwrap().parse()?;
                Ok(thermopack_rs::FlashResult {
                    temperature: t,
                    pressure: p,
                    feed: x.clone(),
                    phase: thermopack_rs::FlashPhase::TwoPhase,
                    liquid: x.clone(),
                    vapor: x.clone(),
                    vapor_fraction: bv,
                    liquid_fraction: 1. - bv,
                    liquid_volume: None,
                    vapor_volume: None,
                    iterations: 0,
                    fugacity_residual: 0.,
                })
            } else {
                model.baseline_tp_flash(t, p, &x)
            };
            let mut roots = [Phase::Liquid, Phase::Vapor];
            if let Some(data) = property_config.iter().find(|op| op.starts_with("roots:")) {
                for (root, label) in roots.iter_mut().zip(data.split(':').skip(1)) {
                    *root = if label == "L" {
                        Phase::Liquid
                    } else {
                        Phase::Vapor
                    };
                }
            }
            match flash.and_then(|f| {
                Ok((
                    f.phase,
                    model.baseline_equilibrium_speed_of_sound_with_phases(&f, roots)?,
                ))
            }) {
                Ok((phase, value)) => println!("{},{value:.17e}", phase as i32),
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if property_config[0] == "V_RESP" {
            match model.baseline_tp_flash(t, p, &x).and_then(|f| {
                let r = model.baseline_equilibrium_response(&f)?.ok_or(
                    thermopack_rs::Error::NotConverged("undefined equilibrium TP response"),
                )?;
                Ok((f.phase, r))
            }) {
                Ok((phase, r)) => println!(
                    "{},{:.17e},{:.17e},{:.17e}",
                    phase as i32, r.heat_capacity, r.volume_dt, r.volume_dp
                ),
                Err(e) => println!("ERROR: {e}"),
            }
            continue;
        }
        if property_config[0] == "V_RR" {
            let inputs: Vec<_> = property_config
                .iter()
                .find(|v| v.starts_with("rr:"))
                .ok_or("missing RR inputs")?
                .split(':')
                .collect();
            let k = inputs[1]
                .split(',')
                .map(str::parse)
                .collect::<Result<Vec<f64>, _>>()?;
            let mut state = thermopack_rs::RachfordRiceIteration::new(x.len());
            state.split.vapor_fraction = inputs[2].parse()?;
            match model.rr_successive_substitution_iteration(
                t,
                p,
                &x,
                &k,
                &mut state,
                inputs[3] == "1",
            ) {
                Ok(()) => {
                    let mut values =
                        vec![f64::from(state.converged), f64::from(state.has_solution)];
                    if state.has_solution {
                        values.extend([
                            state.split.vapor_fraction,
                            state.split.liquid_fraction,
                            state.gibbs,
                        ]);
                        values.extend(
                            state
                                .split
                                .liquid
                                .iter()
                                .chain(&state.split.vapor)
                                .chain(&state.equilibrium_ratios)
                                .chain(&state.liquid_ln_fugacity)
                                .chain(&state.vapor_ln_fugacity),
                        );
                    }
                    println!(
                        "{}",
                        values
                            .iter()
                            .map(|v| format!("{v:.17e}"))
                            .collect::<Vec<_>>()
                            .join(",")
                    );
                }
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if matches!(property_config[0], "PH" | "PS" | "V_PH" | "V_PS") {
            let mut options = model.caloric_options();
            for operation in &property_config[1..] {
                if let Some(bounds) = operation.strip_prefix("bounds:") {
                    let data: Vec<_> = bounds.split(':').collect();
                    options.minimum_temperature = data[0].parse()?;
                    options.maximum_temperature = data[1].parse()?;
                }
            }
            let result = if property_config[0].ends_with("PH") {
                model.ph_flash_with_options(p, &x, t, options)
            } else {
                model.ps_flash_with_options(p, &x, t, options)
            };
            match result {
                Ok(flash) => {
                    print!(
                        "{:.17e},{},{:.17e},{:.17e}",
                        flash.temperature,
                        flash.phase as i32,
                        flash.vapor_fraction,
                        flash.liquid_fraction
                    );
                    for value in flash.liquid.iter().chain(&flash.vapor) {
                        print!(",{value:.17e}");
                    }
                    println!();
                }
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if matches!(
            property_config[0],
            "TP" | "V_TP" | "TP_BASELINE" | "V_TP_BASELINE"
        ) {
            let result = if property_config[0].ends_with("BASELINE") {
                model.baseline_tp_flash(t, p, &x)
            } else {
                model.tp_flash(t, p, &x)
            };
            match result {
                Ok(flash) => {
                    print!(
                        "{},{:.17e},{:.17e}",
                        flash.phase as i32, flash.vapor_fraction, flash.liquid_fraction
                    );
                    for value in flash.liquid.iter().chain(&flash.vapor) {
                        print!(",{value:.17e}");
                    }
                    println!();
                }
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if matches!(property_config[0], "EQ" | "V_EQ") {
            match model
                .tp_flash(t, p, &x)
                .and_then(|f| Ok((f.phase, model.equilibrium_properties(&f)?)))
            {
                Ok((phase, prop)) => println!(
                    "{},{:.17e},{:.17e},{:.17e},{:.17e}",
                    phase as i32,
                    prop.enthalpy,
                    prop.entropy,
                    prop.heat_capacity.unwrap_or(f64::NAN),
                    prop.volume
                ),
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        let state = model.prepare(t, &x)?;
        let v = state.volume(
            p,
            if property_config[0].ends_with('L') {
                Phase::Liquid
            } else {
                Phase::Vapor
            },
        )?;
        if fields[5].starts_with('D') {
            let (mut dt, mut dp, mut dn) = (
                vec![0.0; x.len()],
                vec![0.0; x.len()],
                vec![0.0; x.len() * x.len()],
            );
            state.fugacity_derivatives_into(v, &mut dt, &mut dp, &mut dn)?;
            let values: Vec<_> = dt.into_iter().chain(dp).chain(dn).collect();
            println!(
                "{}",
                values
                    .iter()
                    .map(|v| format!("{v:.17e}"))
                    .collect::<Vec<_>>()
                    .join(",")
            );
            continue;
        }
        let pressure = state.pressure(v)?;
        let residual = state.residual(v)?;
        let mut phi = vec![0.0; x.len()];
        state.ln_fugacity_into(v, &mut phi)?;
        print!(
            "{v:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            pressure.value,
            pressure.dt,
            pressure.dv,
            residual.enthalpy,
            residual.entropy,
            state.enthalpy(v)?,
            state.entropy(v)?,
            state.baseline_enthalpy_derivative_tp(v)?
        );
        for value in phi {
            print!(",{value:.17e}");
        }
        if alpha_properties || translated {
            let (mut dt, mut dp, mut dn) = (
                vec![0.; x.len()],
                vec![0.; x.len()],
                vec![0.; x.len() * x.len()],
            );
            state.fugacity_derivatives_into(v, &mut dt, &mut dp, &mut dn)?;
            for value in dt.iter().chain(&dp).chain(&dn) {
                print!(",{value:.17e}");
            }
        }
        println!();
    }
    Ok(())
}
