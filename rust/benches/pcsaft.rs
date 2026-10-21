use std::hint::black_box;
use std::time::Instant;

/// Measures the state of the PC-SAFT model, which the caloric properties
/// need. The rows report nanoseconds per evaluation and a checksum of the
/// resulting pressure.
fn main() {
    let iterations = 100_000;
    for (label, components, composition) in [
        ("pure", vec!["CO2"], vec![1.0]),
        ("mixture", vec!["CO2", "C1"], vec![0.5, 0.5]),
    ] {
        let model = thermopack_rs::PcSaft::new(&components).unwrap();
        let prepared = model.prepare(300.0, &composition).unwrap();
        let rho = 100.0;
        for _ in 0..1000 {
            let _ = prepared.state(black_box(rho)).unwrap();
        }
        let mut checksum = 0.0;
        let start = Instant::now();
        for _ in 0..iterations {
            let state = prepared.state(black_box(rho)).unwrap();
            checksum += state.pressure;
        }
        println!(
            "state_{label},{:.6},{checksum:.17e}",
            start.elapsed().as_nanos() as f64 / iterations as f64
        );

        for _ in 0..1000 {
            let _ = prepared.pressure(black_box(rho));
        }
        let mut checksum = 0.0;
        let start = Instant::now();
        for _ in 0..iterations {
            checksum += prepared.pressure(black_box(rho));
        }
        println!(
            "pressure_{label},{:.6},{checksum:.17e}",
            start.elapsed().as_nanos() as f64 / iterations as f64
        );
    }

    // The chain rule of the hard-sphere model of Boublik, which the state of
    // a mixture needs.
    let model = thermopack_rs::PcSaft::new(&["CO2", "C1"]).unwrap();
    let moles = [0.05, 0.05];
    for _ in 0..1000 {
        let _ = thermopack_rs::pc_saft::pcsaft_chain_energy(&model, 300.0, 1.0e-3, &moles).unwrap();
    }
    let mut checksum = 0.0;
    let start = Instant::now();
    for _ in 0..iterations {
        let energy =
            thermopack_rs::pc_saft::pcsaft_chain_energy(&model, 300.0, 1.0e-3, &moles).unwrap();
        checksum += energy.value;
    }
    println!(
        "chain_mixture,{:.6},{checksum:.17e}",
        start.elapsed().as_nanos() as f64 / iterations as f64
    );
}
