use std::{env, fs, hint::black_box, time::Instant};
use thermopack_rs::{FlashPhase, LennardJonesKind, Multiparameter};

// Manifest rows: profile mode kind T P specification guess guess_flag Tmin Tmax, then source outputs.
fn model(profile: usize) -> Multiparameter {
    match profile {
        0 => Multiparameter::new("C3").unwrap(),
        1 => Multiparameter::new("R134A").unwrap(),
        2 => Multiparameter::lennard_jones("KR", LennardJonesKind::Full).unwrap(),
        3 => Multiparameter::lennard_jones("AR", LennardJonesKind::Full).unwrap(),
        4 => Multiparameter::lennard_jones("KR", LennardJonesKind::TruncatedShifted).unwrap(),
        5 => Multiparameter::reduced_lennard_jones_truncated_shifted(),
        _ => panic!("invalid profile"),
    }
}
fn main() {
    let mut args = env::args().skip(1);
    let Some(path) = args.next() else {
        // `cargo test --all-targets` executes harness=false benches without args.
        // Keep that smoke invocation harmless; the Python driver supplies a manifest.
        eprintln!("meos_flash benchmark smoke invocation; pass a manifest to run workloads");
        return;
    };
    let rest: Vec<_> = args.collect();
    assert!(
        rest.is_empty() || rest == ["--check"],
        "unknown benchmark arguments"
    );
    let check = !rest.is_empty();
    let input = fs::read_to_string(path).unwrap();
    let mut lines = input.lines();
    let declared: usize = lines.next().expect("manifest row count").parse().unwrap();
    let rows: Vec<Vec<f64>> = lines
        .map(|l| l.split_whitespace().map(|x| x.parse().unwrap()).collect())
        .collect();
    assert_eq!(rows.len(), declared, "manifest row count mismatch");
    assert!(
        rows.iter()
            .all(|r| r.len() == 18 && r.iter().all(|x| x.is_finite())),
        "invalid manifest rows"
    );
    let mut begin = 0;
    while begin < rows.len() {
        let profile = rows[begin][0] as usize;
        let end = rows[begin..]
            .iter()
            .position(|r| r[0] as usize != profile)
            .map_or(rows.len(), |n| begin + n);
        let m = model(profile);
        let mut prepared = Vec::with_capacity(end - begin);
        for row in &rows[begin..end] {
            prepared.push(if row[2] as usize == 5 {
                Some(m.prepare(row[3]).unwrap())
            } else {
                None
            });
        }
        for mode in 0..6 {
            let selected: Vec<_> = (begin..end)
                .filter(|&i| rows[i][1] as usize == mode)
                .collect();
            if selected.is_empty() {
                continue;
            }
            let solve = |idx: usize| -> [f64; 8] {
                let r = &rows[idx];
                let (kind, t, p, spec, guess, has_guess, tmin, tmax) = (
                    r[2] as usize,
                    r[3],
                    r[4],
                    r[5],
                    r[6],
                    r[7] != 0.,
                    r[8],
                    r[9],
                );
                if kind == 5 {
                    let pr = prepared[idx - begin].as_ref().unwrap();
                    let vl = r[14];
                    let vv = r[15];
                    return [
                        t,
                        p,
                        FlashPhase::Single as i32 as f64,
                        -1.,
                        pr.enthalpy(vl).unwrap(),
                        pr.entropy(vl).unwrap(),
                        pr.enthalpy(vv).unwrap(),
                        pr.entropy(vv).unwrap(),
                    ];
                }
                let q = match kind {
                    0 => m.tp_flash(t, p).unwrap(),
                    1 => m
                        .ph_flash_with_options(
                            p,
                            spec,
                            thermopack_rs::CaloricFlashOptions {
                                temperature_guess: has_guess.then_some(guess),
                                minimum_temperature: tmin,
                                maximum_temperature: tmax,
                                ..Default::default()
                            },
                        )
                        .unwrap(),
                    2 => m
                        .ps_flash_with_options(
                            p,
                            spec,
                            thermopack_rs::CaloricFlashOptions {
                                temperature_guess: has_guess.then_some(guess),
                                minimum_temperature: tmin,
                                maximum_temperature: tmax,
                                ..Default::default()
                            },
                        )
                        .unwrap(),
                    _ => unreachable!(),
                };
                let props = m.equilibrium_properties(&q).unwrap();
                [
                    q.temperature,
                    q.pressure,
                    q.phase as i32 as f64,
                    q.vapor_fraction,
                    q.liquid_volume.unwrap_or(0.),
                    q.vapor_volume.unwrap_or(0.),
                    props.enthalpy,
                    props.entropy,
                ]
            };
            if check {
                for (n, &i) in selected.iter().enumerate() {
                    let q = solve(i);
                    println!(
                        "check_{profile}_{mode}_{n},{}",
                        q.map(|x| format!("{x:.17e}")).join(",")
                    );
                }
            } else {
                let iterations = 1000;
                let start = Instant::now();
                let mut sum = [0.; 8];
                for n in 0..iterations {
                    for (s, x) in sum
                        .iter_mut()
                        .zip(black_box(solve(selected[n % selected.len()])))
                    {
                        *s += x;
                    }
                }
                println!(
                    "native_{profile}_{mode},{:.6},{}",
                    start.elapsed().as_nanos() as f64 / iterations as f64,
                    sum.map(|x| format!("{x:.17e}")).join(",")
                );
            }
        }
        begin = end;
    }
}
