// Line-oriented reference protocol: KIND|COMP,COMP|T|P|x,x|L or V
use std::io::{self, BufRead};
use thermopack_rs::{Cubic, CubicKind, Phase};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for line in io::stdin().lock().lines() {
        let line = line?;
        let fields: Vec<_> = line.split('|').collect();
        if fields.len() != 6 {
            return Err("expected six fields".into());
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
        let model = Cubic::new(kind, &fields[1].split(',').collect::<Vec<_>>())?;
        let t: f64 = fields[2].parse()?;
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
        if fields[5] == "PH" || fields[5] == "PS" {
            let result = if fields[5] == "PH" {
                model.ph_flash(p, &x, t)
            } else {
                model.ps_flash(p, &x, t)
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
        if fields[5] == "TP" {
            match model.tp_flash(t, p, &x) {
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
        if fields[5] == "EQ" {
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
            if fields[5].ends_with('L') {
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
        println!();
    }
    Ok(())
}
