//! Per-case diagnostics for the shared MEOS flash benchmark manifest.
use std::{env, fs, hint::black_box, time::Instant};
use thermopack_rs::{CaloricFlashOptions, LennardJonesKind, Multiparameter};

fn main() {
    let args: Vec<_> = env::args().skip(1).collect();
    if args.is_empty() {
        eprintln!("usage: meos_caloric_profile MANIFEST [PROFILE=0] [REPEATS=50]");
        return;
    }
    let profile: usize = args.get(1).map_or(0, |s| s.parse().unwrap());
    let repeats: usize = args.get(2).map_or(50, |s| s.parse().unwrap());
    assert!(repeats > 0);
    let model = match profile {
        0 => Multiparameter::new("C3").unwrap(),
        1 => Multiparameter::new("R134A").unwrap(),
        2 => Multiparameter::lennard_jones("KR", LennardJonesKind::Full).unwrap(),
        3 => Multiparameter::lennard_jones("AR", LennardJonesKind::Full).unwrap(),
        4 => Multiparameter::lennard_jones("KR", LennardJonesKind::TruncatedShifted).unwrap(),
        5 => Multiparameter::reduced_lennard_jones_truncated_shifted(),
        _ => panic!("invalid profile"),
    };
    let input = fs::read_to_string(&args[0]).unwrap();
    let mut lines = input.lines();
    let declared: usize = lines.next().unwrap().parse().unwrap();
    let rows: Vec<Vec<f64>> = lines
        .map(|l| l.split_whitespace().map(|v| v.parse().unwrap()).collect())
        .collect();
    assert_eq!(rows.len(), declared);
    let mut counters = [0; 6];
    println!("case,ns,iterations,phase,T,P,target,guess,Tmin,Tmax,expected_T");
    for r in rows {
        assert_eq!(r.len(), 18);
        if r[0] as usize != profile {
            continue;
        }
        let lane = r[1] as usize;
        let index = counters[lane];
        counters[lane] += 1;
        if !(1..=4).contains(&lane) {
            continue;
        }
        let options = CaloricFlashOptions {
            temperature_guess: (r[7] != 0.).then_some(r[6]),
            minimum_temperature: r[8],
            maximum_temperature: r[9],
            ..Default::default()
        };
        let solve = || {
            if r[2] == 1. {
                model.ph_flash_with_options(r[4], r[5], options)
            } else {
                model.ps_flash_with_options(r[4], r[5], options)
            }
            .unwrap()
        };
        let result = solve();
        // Include the same returned-property evaluation as the shared benchmark.
        let start = Instant::now();
        for _ in 0..repeats {
            let q = black_box(solve());
            black_box(model.equilibrium_properties(&q).unwrap());
        }
        let elapsed = start.elapsed().as_nanos() as f64 / repeats as f64;
        println!(
            "check_{profile}_{lane}_{index},{elapsed:.3},{},{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            result.iterations,
            result.phase as i32,
            result.temperature,
            r[4],
            r[5],
            r[6],
            r[8],
            r[9],
            r[10],
        );
    }
}
