use std::hint::black_box;
use std::time::Instant;
use thermopack_rs::{Cubic, CubicKind, Phase};

fn main() {
    let model = Cubic::new(CubicKind::PengRobinson, &["CO2", "C1", "N2"]).unwrap();
    let state = model.prepare(300.0, &[0.7, 0.2, 0.1]).unwrap();
    let iterations = 1_000_000;
    for prepare in [false, true] {
        let mut checksum = [0.0; 3];
        let start = Instant::now();
        for i in 0..iterations {
            let prepared;
            let current = if prepare {
                prepared = model
                    .prepare(black_box(300.0), black_box(&[0.7, 0.2, 0.1]))
                    .unwrap();
                &prepared
            } else {
                &state
            };
            let p = current
                .pressure(black_box(0.0001 + (i % 1000) as f64 * 1e-6))
                .unwrap();
            checksum[0] += p.value;
            checksum[1] += p.dt;
            checksum[2] += p.dv;
        }
        println!(
            "pressure_{},{:.6},{:.17e},{:.17e},{:.17e}",
            if prepare { "with_mixing" } else { "prepared" },
            start.elapsed().as_nanos() as f64 / iterations as f64,
            checksum[0],
            checksum[1],
            checksum[2]
        );
    }
    for prepare in [false, true] {
        let mut checksum = 0.0;
        let start = Instant::now();
        for i in 0..iterations {
            let prepared;
            let current = if prepare {
                prepared = model
                    .prepare(black_box(300.0), black_box(&[0.7, 0.2, 0.1]))
                    .unwrap();
                &prepared
            } else {
                &state
            };
            checksum += current
                .volume(black_box(1e6 + (i % 1000) as f64 * 1000.0), Phase::Vapor)
                .unwrap();
        }
        println!(
            "volume_{},{:.6},{checksum:.17e}",
            if prepare { "with_mixing" } else { "prepared" },
            start.elapsed().as_nanos() as f64 / iterations as f64
        );
    }
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C2", "C3"]).unwrap();
    let flash_iterations = 10000;
    let mut checksum = [0.0; 3];
    let start = Instant::now();
    for i in 0..flash_iterations {
        let t = 180.0 + (i % 100) as f64 * 0.8;
        let result = model
            .tp_flash(black_box(t), black_box(1e6), black_box(&[0.5, 0.3, 0.2]))
            .unwrap();
        checksum[0] += result.vapor_fraction;
        checksum[1] += result.liquid_fraction;
        checksum[2] += result.phase as i32 as f64;
    }
    println!(
        "tp_flash,{:.6},{:.17e},{:.17e},{:.17e}",
        start.elapsed().as_nanos() as f64 / flash_iterations as f64,
        checksum[0],
        checksum[1],
        checksum[2]
    );
    let mut targets = Vec::new();
    for i in 0..100 {
        let t = 180.0 + i as f64 * 0.8;
        let properties = model
            .equilibrium_properties(&model.tp_flash(t, 1e6, &[0.5, 0.3, 0.2]).unwrap())
            .unwrap();
        let ph = model
            .ph_flash(1e6, &[0.5, 0.3, 0.2], properties.enthalpy)
            .unwrap();
        let ps = model
            .ps_flash(1e6, &[0.5, 0.3, 0.2], properties.entropy)
            .unwrap();
        assert!((ph.temperature - t).abs() < 1e-5 && (ps.temperature - t).abs() < 1e-5);
        targets.push((properties.enthalpy, properties.entropy));
    }
    let caloric_iterations = 1000;
    for enthalpy in [true, false] {
        let start = Instant::now();
        let mut checksum = [0.0; 3];
        for i in 0..caloric_iterations {
            let target = targets[i % 100];
            let flash = if enthalpy {
                model.ph_flash(
                    black_box(1e6),
                    black_box(&[0.5, 0.3, 0.2]),
                    black_box(target.0),
                )
            } else {
                model.ps_flash(
                    black_box(1e6),
                    black_box(&[0.5, 0.3, 0.2]),
                    black_box(target.1),
                )
            }
            .unwrap();
            checksum[0] += flash.temperature;
            checksum[1] += flash.vapor_fraction;
            checksum[2] += flash.phase as i32 as f64;
        }
        println!(
            "{}_flash,{:.6},{:.17e},{:.17e},{:.17e}",
            if enthalpy { "ph" } else { "ps" },
            start.elapsed().as_nanos() as f64 / caloric_iterations as f64,
            checksum[0],
            checksum[1],
            checksum[2]
        );
    }
}
