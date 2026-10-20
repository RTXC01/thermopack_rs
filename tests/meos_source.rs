//! Source-rounding snapshot runs in the oracle process-main context.
//! Windows CRT exp may differ by an ulp on new threads; critical Cp amplifies it.
use thermopack_rs::helmholtz::PureHelmholtz;
use thermopack_rs::{HelmholtzJet, LennardJonesKind, Multiparameter, Phase};

fn model(config: &str) -> Multiparameter {
    let c: Vec<_> = config.split(':').collect();
    match c[0] {
        "LJ_MEOS" => Multiparameter::lennard_jones(c[1], LennardJonesKind::Full).unwrap(),
        "LJTS_MEOS" => {
            Multiparameter::lennard_jones(c[1], LennardJonesKind::TruncatedShifted).unwrap()
        }
        "LJTS_REDUCED" => Multiparameter::reduced_lennard_jones_truncated_shifted(),
        _ => Multiparameter::new(c[1]).unwrap(),
    }
}
fn jet(a: HelmholtzJet) -> [f64; 6] {
    [
        a.alpha,
        a.d_tau,
        a.d_tautau,
        a.d_delta,
        a.d_deltadelta,
        a.d_deltatau,
    ]
}
fn extensive(a: PureHelmholtz) -> [f64; 10] {
    [
        a.value, a.dt, a.dv, a.dn, a.dtt, a.dtv, a.dtn, a.dvv, a.dvn, a.dnn,
    ]
}
fn main() {
    let atol = [
        1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-11, 1e-5,
        1e-7, 1e-2, 1e-6, 1e-9, 1e-6, 1e-6, 1e-6, 1e-9, 1e-8, 1e-10, 1e-8, 1e-10, 1e-10, 1e-6,
        1e-10, 1e-10, 1e-6, 1e-10, 1e-1, 1e-6, 1e-10, 1e-10, 1e-10, 1e-6, 1e-10, 1e-10, 1e-6,
        1e-10, 1e-1, 1e-6, 1e-10, 1e-10, 1e-10,
    ];
    let mut cache = None;
    let mut previous = "";
    let mut count = 0;
    let mut absent = [0; 4];
    for row in include_str!("fixtures/meos.csv").lines().skip(1) {
        let c: Vec<_> = row.split(';').collect();
        let f: Vec<_> = c[0].split('|').collect();
        if f[1] != previous {
            cache = Some(model(f[1]));
            previous = f[1];
        }
        let eos = cache.as_ref().unwrap();
        let t = f[2].parse().unwrap();
        let v: f64 = f[3].parse().unwrap();
        let n: f64 = f[4].parse().unwrap();
        let prep = eos.prepare(t).unwrap();
        let s = prep.state(v).unwrap();
        let mut actual = jet(s.ideal).to_vec();
        actual.extend(jet(s.residual));
        actual.extend([
            s.pressure.value,
            s.pressure.dt,
            s.pressure.dv,
            s.internal_energy,
            s.entropy,
            s.enthalpy,
            s.helmholtz_energy,
            s.gibbs_energy,
            s.heat_capacity_v,
            s.heat_capacity_p,
            s.ln_fugacity.unwrap_or(f64::NAN),
            s.baseline_speed_of_sound.unwrap_or(f64::NAN),
        ]);
        actual.extend(extensive(prep.residual_helmholtz(v * n, n).unwrap()));
        actual.extend(extensive(prep.ideal_helmholtz(v * n, n).unwrap()));
        actual.extend([
            prep.sat_delta_estimate(Phase::Liquid).unwrap_or(f64::NAN),
            prep.sat_delta_estimate(Phase::Vapor).unwrap_or(f64::NAN),
        ]);
        for (i, a) in actual.into_iter().enumerate() {
            if c[i + 1].is_empty() {
                assert!(!a.is_finite(), "{} field {i}", c[0]);
                absent[match i {
                    22 => 0,
                    23 => 1,
                    44 => 2,
                    45 => 3,
                    _ => panic!("unexpected non-finite source field"),
                }] += 1;
            } else {
                let b: f64 = c[i + 1].parse().unwrap();
                assert!(
                    (a - b).abs() <= atol[i] + 2e-10 * a.abs().max(b.abs()),
                    "{} field {i}: {a} != {b}",
                    c[0]
                );
            }
        }
        count += 1;
    }
    assert_eq!(count, 5082);
    assert_eq!(absent, [465, 231, 330, 330]);
    println!(
        "5082 pure MEOS states / 46 fields verified in process-main context; non-finite source fields retained"
    );
    if cfg!(windows) {
        let cp = || {
            let m = Multiparameter::new("C3").unwrap();
            m.prepare(m.tc())
                .unwrap()
                .state(1. / m.rc())
                .unwrap()
                .heat_capacity_p
        };
        println!(
            "Windows critical Cp rounding audit: main={}, worker={}",
            cp(),
            std::thread::spawn(cp).join().unwrap()
        );
    }
}
