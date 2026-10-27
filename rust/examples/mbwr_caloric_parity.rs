use std::io::{self, BufRead};
use thermopack_rs::{EntropyReference, IdealGas, Mbwr, MbwrFluid, component};
#[path = "../support/mbwr_caloric_coordinates.rs"]
mod protocol;

fn main() {
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        let f: Vec<_> = line.split('|').collect();
        assert_eq!(f.len(), 6);
        let t = f[1].parse().unwrap();
        let v = f[2].parse().unwrap();
        let n = f[3].parse().unwrap();
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
        for result in [
            model.baseline_caloric_tv(t, v, n),
            model.caloric_tv(t, v, n),
        ] {
            match result {
                Ok(q) => println!(
                    "{}",
                    protocol::coordinates(&q)
                        .iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
                Err(e) => println!("ERROR: {e}"),
            }
        }
    }
}
