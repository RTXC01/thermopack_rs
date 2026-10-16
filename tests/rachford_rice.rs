use thermopack_rs::numerics::{RachfordRiceState, rachford_rice, rachford_rice_extended};

#[test]
fn negative_flash_and_independent_liquid_fraction() {
    for (z, expected) in [([0.2, 0.8], -0.4), ([0.8, 0.2], 1.4)] {
        let mut state = RachfordRiceState::new(2);
        assert!(rachford_rice_extended(&z, &[2., 0.5], &mut state, false).unwrap());
        assert!((state.vapor_fraction - expected).abs() < 1e-12);
        assert!((state.liquid.iter().sum::<f64>() - 1.).abs() < 1e-12);
        assert!((state.vapor.iter().sum::<f64>() - 1.).abs() < 1e-12);
        for (i, &zi) in z.iter().enumerate() {
            assert!(
                (state.liquid_fraction * state.liquid[i] + state.vapor_fraction * state.vapor[i]
                    - zi)
                    .abs()
                    < 1e-12
            );
        }
        assert_eq!(
            rachford_rice(&z, &[2., 0.5]).unwrap(),
            if expected < 0. { 0. } else { 1. }
        );
    }
    let mut state = RachfordRiceState::new(2);
    assert!(rachford_rice_extended(&[1e-20, 1.], &[1e-20, 2.], &mut state, false).unwrap());
    assert_eq!(state.vapor_fraction, 1.);
    assert!(state.liquid_fraction > 0. && state.liquid_fraction < 1e-19);
    assert!((state.liquid[0] - 0.5).abs() < 1e-12);
}

#[test]
fn buffers_reused_and_no_root_is_not_a_solution() {
    let mut state = RachfordRiceState::new(2);
    state.liquid = [0.3, 0.7].to_vec();
    state.vapor = [0.8, 0.2].to_vec();
    assert!(!rachford_rice_extended(&[0.5, 0.5], &[1., 1.], &mut state, false).unwrap());
    assert_eq!(state.liquid, [0.3, 0.7]);
    assert_eq!(state.vapor, [0.8, 0.2]);
    assert!(rachford_rice_extended(&[0.5, 0.5], &[f64::NAN, 1.], &mut state, false).is_err());
    assert!(rachford_rice_extended(&[0.5, 0.5], &[1., 1.], &mut state, true).is_err());
}

#[test]
fn extended_and_one_step_rr_match_exact_source() {
    let mut count = 0;
    for row in include_str!("fixtures/rachford_rice.csv").lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let z: Vec<f64> = fields[1].split(',').map(|v| v.parse().unwrap()).collect();
        let k: Vec<f64> = fields[2].split(',').map(|v| v.parse().unwrap()).collect();
        let mut state = RachfordRiceState::new(z.len());
        state.vapor_fraction = fields[3].parse().unwrap();
        let success = rachford_rice_extended(&z, &k, &mut state, fields[4] == "1").unwrap();
        let mut actual = vec![if success { 1. } else { 0. }];
        if success {
            actual.extend([state.vapor_fraction, state.liquid_fraction]);
            actual.extend(state.liquid.iter().chain(&state.vapor));
        }
        assert_eq!(actual.len(), columns.len() - 1);
        for (i, (a, b)) in actual.iter().zip(&columns[1..]).enumerate() {
            let b: f64 = b.parse().unwrap();
            let tolerance = if i == 0 {
                0.
            } else {
                1e-12 + 2e-10 * a.abs().max(b.abs())
            };
            assert!(
                (a - b).abs() <= tolerance,
                "{} field {i}: {a}, {b}",
                columns[0]
            );
        }
        count += 1;
    }
    assert_eq!(count, 698);
}
