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
        let kind = match fields[0] {
            "PR" => CubicKind::PengRobinson,
            "SRK" => CubicKind::SoaveRedlichKwong,
            "RK" => CubicKind::RedlichKwong,
            "VDW" => CubicKind::VanDerWaals,
            _ => return Err("unknown model".into()),
        };
        let model = Cubic::new(kind, &fields[1].split(',').collect::<Vec<_>>())?;
        let t: f64 = fields[2].parse()?;
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
        let state = model.prepare(t, &x)?;
        let v = state.volume(
            p,
            if fields[5] == "L" {
                Phase::Liquid
            } else {
                Phase::Vapor
            },
        )?;
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
