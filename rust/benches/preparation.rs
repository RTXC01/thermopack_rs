use std::hint::black_box;
use std::time::Instant;
use thermopack_rs::{Cubic, CubicKind, CubicOptions, MixingRule};

/// Measures the mixture preparation, which dominates the cost of a single
/// cubic property evaluation. The rows report nanoseconds per preparation and
/// a checksum of the resulting pressure.
fn main() {
    let cases = [
        (
            "pr_vdw",
            CubicKind::PengRobinson,
            MixingRule::VdW,
            "DEFAULT",
        ),
        (
            "srk_vdw",
            CubicKind::SoaveRedlichKwong,
            MixingRule::VdW,
            "DEFAULT",
        ),
        ("pr_tcpr", CubicKind::PengRobinson, MixingRule::VdW, "tcPR"),
        (
            "srk_hv1",
            CubicKind::SoaveRedlichKwong,
            MixingRule::HuronVidal1,
            "DEFAULT",
        ),
        (
            "srk_unifac",
            CubicKind::SoaveRedlichKwong,
            MixingRule::Unifac,
            "DEFAULT",
        ),
    ];
    let iterations = 200_000;
    for (label, kind, rule, reference) in cases {
        let model = Cubic::with_options(
            kind,
            &["CO2", "C1", "N2"],
            CubicOptions {
                mixing_rule: rule,
                parameter_reference: reference,
                ..Default::default()
            },
        )
        .unwrap();
        let composition = [0.7, 0.2, 0.1];
        for _ in 0..1000 {
            let state = model.prepare(300.0, &composition).unwrap();
            let _ = state.pressure(1.0e-3).unwrap();
        }
        let mut checksum = 0.0;
        let start = Instant::now();
        for _ in 0..iterations {
            let state = model
                .prepare(black_box(300.0), black_box(&composition))
                .unwrap();
            checksum += state.pressure(black_box(1.0e-3)).unwrap().value;
        }
        println!(
            "{label},{:.6},{checksum:.17e}",
            start.elapsed().as_nanos() as f64 / iterations as f64
        );
    }

    // The Wong-Sandler rule needs the binary parameters of the pairs.
    let mut data = thermopack_rs::wong_sandler::WongSandlerData::new(3);
    for (first, second) in [(0usize, 1usize), (0, 2), (1, 2)] {
        let parameter = thermopack_rs::wong_sandler::Fraction::constant(0.05);
        data.set_kij(first, second, parameter).unwrap();
        data.set_kij(second, first, parameter).unwrap();
    }
    let model = Cubic::with_options(
        CubicKind::PengRobinson,
        &["CO2", "C1", "N2"],
        CubicOptions {
            mixing_rule: MixingRule::WongSandler,
            wong_sandler: Some(&data),
            ..Default::default()
        },
    )
    .unwrap();
    let composition = [0.7, 0.2, 0.1];
    for _ in 0..1000 {
        let state = model.prepare(300.0, &composition).unwrap();
        let _ = state.pressure(1.0e-3).unwrap();
    }
    let mut checksum = 0.0;
    let start = Instant::now();
    for _ in 0..iterations {
        let state = model
            .prepare(black_box(300.0), black_box(&composition))
            .unwrap();
        checksum += state.pressure(black_box(1.0e-3)).unwrap().value;
    }
    println!(
        "pr_ws,{:.6},{checksum:.17e}",
        start.elapsed().as_nanos() as f64 / iterations as f64
    );
}
