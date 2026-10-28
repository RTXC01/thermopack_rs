use std::io::{self, BufRead};
use thermopack_rs::{Mbwr, Saturation};

fn main() {
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        let q: Vec<_> = line.split('|').collect();
        assert_eq!(q.len(), 8);
        let model = Mbwr::from_record(q[0]).unwrap();
        let spec: f64 = q[2].parse().unwrap();
        let maximum: f64 = q[3].parse().unwrap();
        if q[1] == "R" {
            match model
                .prepare_equilibrium(spec, maximum)
                .and_then(|e| e.densities(q[4].parse().unwrap()))
            {
                Ok(roots) => println!(
                    "ROOTS;{}",
                    roots
                        .iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
                Err(e) => println!("ERROR: {e}"),
            }
            continue;
        }
        let result = if q[1] == "T" {
            model.saturation_temperature(
                spec,
                q[4].parse().unwrap(),
                q[5].parse().unwrap(),
                maximum,
            )
        } else if q[6] != "0" {
            model.saturation_pressure_bounded(
                spec,
                q[6].parse().unwrap(),
                q[7].parse().unwrap(),
                maximum,
            )
        } else {
            model.saturation_pressure(spec, maximum)
        };
        match result {
            Ok(Saturation {
                temperature,
                pressure,
                liquid_volume,
                vapor_volume,
                dpdt,
            }) => {
                let roots = model
                    .prepare_equilibrium(temperature, maximum)
                    .unwrap()
                    .densities(pressure)
                    .unwrap();
                println!(
                    "{temperature:.17e},{pressure:.17e},{liquid_volume:.17e},{vapor_volume:.17e},{dpdt:.17e};{}",
                    roots
                        .iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                );
            }
            Err(e) => println!("ERROR: {e}"),
        }
    }
}
