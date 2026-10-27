use std::io::{self, BufRead};
use thermopack_rs::{
    EntropyReference, IdealGas, Mbwr, MbwrDensityParameters, MbwrFluid, Phase, component,
};
#[path = "../support/mbwr_total_tp_coordinates.rs"]
mod protocol;
fn emit(label: &str, result: thermopack_rs::Result<Vec<f64>>) {
    match result {
        Ok(v) => println!(
            "{label} {}",
            v.iter()
                .map(|v| format!("{v:.17e}"))
                .collect::<Vec<_>>()
                .join(",")
        ),
        Err(e) => println!("{label} ERROR: {e}"),
    }
}
fn main() {
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        let f: Vec<_> = line.split('|').collect();
        assert_eq!(f.len(), 10);
        let t = f[1].parse().unwrap();
        let v = f[2].parse().unwrap();
        let n = f[3].parse().unwrap();
        let p = f[7].parse().unwrap();
        let phase = if f[8] == "1" {
            Phase::Liquid
        } else if f[8] == "2" {
            Phase::Vapor
        } else {
            Phase::MinimumGibbs
        };
        let mut model = if f[4] == "DEFAULT" {
            MbwrFluid::from_record(f[0]).unwrap()
        } else {
            let c = component(f[4]).unwrap();
            MbwrFluid::with_ideal_gas(
                Mbwr::from_record(f[0]).unwrap(),
                IdealGas::with_integration_constants(c, 0., 0.).unwrap(),
                c.molar_mass * 1e-3,
            )
            .unwrap()
        };
        let profile: usize = f[5].parse().unwrap();
        if profile == 1 || profile == 3 {
            model
                .ideal_gas_mut()
                .set_formation_enthalpy(-12345.678)
                .unwrap();
            model
                .ideal_gas_mut()
                .set_standard_entropy(210.25, EntropyReference::Atmosphere)
                .unwrap();
        }
        if profile == 2 || profile == 3 {
            model.ideal_gas_mut().set_heat_capacity(8, &[3.5]).unwrap();
        }
        if f[6] == "SOURCE" {
            let a: Vec<f64> = f[9].split(',').map(|s| s.parse().unwrap()).collect();
            match model.baseline_caloric_tp_with_parameters(
                t,
                p,
                n,
                phase,
                MbwrDensityParameters::new(a[0], a[1], a[2], a[3]).unwrap(),
            ) {
                Ok(q) => {
                    let r = q.iteration;
                    emit(
                        "DENSITY",
                        Ok(vec![
                            r.density,
                            if r.phase == Phase::Liquid { 1. } else { 2. },
                            r.pressure,
                            r.pressure_density_derivative,
                        ]),
                    );
                    emit("CLASSICAL", q.properties.map(|p| protocol::source(&p)));
                }
                Err(e) => {
                    println!("DENSITY ERROR: {e}");
                    println!("CLASSICAL ERROR: {e}");
                }
            }
        } else if f[6] == "BOUNDED" {
            emit(
                "PHYSICAL",
                model
                    .caloric_tp(t, p, n, phase, f[9].parse().unwrap())
                    .map(|q| protocol::physical(&q)),
            );
        } else {
            emit(
                "BASELINE",
                model
                    .baseline_caloric_tvp(t, v, n)
                    .map(|q| protocol::physical(&q)),
            );
            emit(
                "PHYSICAL",
                model.caloric_tvp(t, v, n).map(|q| protocol::physical(&q)),
            );
            emit(
                "CLASSICAL",
                model
                    .baseline_caloric_tp_at_volume(t, p, v, n)
                    .map(|q| protocol::source(&q)),
            );
        }
    }
}
