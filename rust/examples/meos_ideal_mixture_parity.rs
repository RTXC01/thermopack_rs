//! Line-oriented residual mixture qualification protocol (not a user interface).
use std::io::{self, BufRead};
use thermopack_rs::{MultiparameterIdealMixture, Phase};
#[path = "../support/meos_ideal_mixture_coordinates.rs"]
mod protocol;

fn main() {
    for line in io::stdin().lock().lines() {
        let line = line.unwrap();
        let fields: Vec<_> = line.split('|').collect();
        assert_eq!(fields.len(), 5);
        let t = fields[0].parse().unwrap();
        let p = fields[1].parse().unwrap();
        let phase = match fields[2] {
            "1" => Phase::Liquid,
            "2" => Phase::Vapor,
            _ => panic!("phase must be 1 or 2"),
        };
        let components: Vec<_> = fields[3].split(',').collect();
        let x: Vec<_> = fields[4].split(',').map(|s| s.parse().unwrap()).collect();
        let model = MultiparameterIdealMixture::new(&components).unwrap();
        match model.residual_tp(t, p, &x, phase) {
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
