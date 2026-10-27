use thermopack_rs::{MultiparameterIdealMixture, Phase};
#[path = "../rust/support/meos_ideal_mixture_coordinates.rs"]
mod protocol;

#[test]
fn original_internal_mixture_kernels_replay_without_fortran() {
    let mut cases = 0;
    for row in include_str!("fixtures/meos_ideal_mixture.csv")
        .lines()
        .skip(1)
    {
        let row: Vec<_> = row.split(';').collect();
        let case: Vec<_> = row[0].split('|').collect();
        let components: Vec<_> = case[3].split(',').collect();
        let x: Vec<_> = case[4].split(',').map(|s| s.parse().unwrap()).collect();
        let phase = if case[2] == "1" {
            Phase::Liquid
        } else {
            Phase::Vapor
        };
        let q = MultiparameterIdealMixture::new(&components)
            .unwrap()
            .residual_tp(
                case[0].parse().unwrap(),
                case[1].parse().unwrap(),
                &x,
                phase,
            )
            .unwrap_or_else(|e| panic!("{}: {e}", row[0]));
        let actual = protocol::coordinates(&q);
        let expected: Vec<f64> = row[1].split(',').map(|s| s.parse().unwrap()).collect();
        let floors: Vec<f64> = row[2].split(',').map(|s| s.parse().unwrap()).collect();
        assert_eq!(actual.len(), expected.len());
        assert_eq!(actual.len(), floors.len());
        assert_eq!(actual.len(), 13 + 8 * x.len() + x.len() * x.len());
        for (i, ((a, b), floor)) in actual.into_iter().zip(expected).zip(floors).enumerate() {
            assert!(
                a.is_finite() && (a - b).abs() <= floor + 2e-8 * a.abs().max(b.abs()),
                "{} coordinate {i}: {a} != {b}",
                row[0]
            );
        }
        cases += 1;
    }
    assert_eq!(cases, 734);
}
