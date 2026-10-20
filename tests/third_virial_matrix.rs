//! Third virial coefficient matrix comparisons.
//! Copyright thermotools and Contributors. Apache-2.0.

#![allow(clippy::excessive_precision)]

use thermopack_rs::Cubic;
use thermopack_rs::CubicKind;

#[test]
fn the_pure_entries_match_the_third_virial_coefficient() {
    // The diagonal of the matrix holds the third virial coefficient of the
    // pure fluids.
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C3"]).unwrap();
    let matrix = model.binary_third_virial_coefficient_matrix(300.0).unwrap();
    for index in 0..2 {
        let mut composition = vec![0.0; 2];
        composition[index] = 1.0;
        let state = model.prepare(300.0, &composition).unwrap();
        let (_, third) = state.virial_coefficients().unwrap();
        assert!(
            (matrix[index * 2 + index] - third).abs() < 1.0e-18,
            "{index}: {} {third}",
            matrix[index * 2 + index]
        );
    }
}

#[test]
fn the_cross_entries_differ() {
    // The linear system of the original is not symmetric, so the two cross
    // coefficients of a pair differ.
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C3"]).unwrap();
    let matrix = model.binary_third_virial_coefficient_matrix(300.0).unwrap();
    assert!((matrix[1] - matrix[2]).abs() > 0.0);
}

#[test]
fn the_matrix_needs_a_binary_model() {
    let model = Cubic::new(CubicKind::PengRobinson, &["C1", "C3", "N2"]).unwrap();
    let error = model
        .binary_third_virial_coefficient_matrix(300.0)
        .unwrap_err();
    assert!(matches!(error, thermopack_rs::Error::InvalidInput(_)));
}
