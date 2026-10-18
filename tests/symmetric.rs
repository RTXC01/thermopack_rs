use thermopack_rs::numerics::solve_symmetric_in_place;

#[test]
fn upper_factorization_pivots_and_solutions_match_pinned_lapack() {
    let mut counts = [0; 2];
    for row in include_str!("fixtures/symmetric.csv").lines().skip(1) {
        let columns: Vec<_> = row.split(';').collect();
        let fields: Vec<_> = columns[0].split('|').collect();
        let mut matrix: Vec<f64> = fields[1].split(',').map(|s| s.parse().unwrap()).collect();
        let mut rhs: Vec<f64> = fields[2].split(',').map(|s| s.parse().unwrap()).collect();
        let original = matrix.clone();
        let original_rhs = rhs.clone();
        let n = rhs.len();
        let mut pivots = vec![0; n];
        let mut work = vec![0.; n];
        let success =
            solve_symmetric_in_place(&mut matrix, &mut rhs, &mut pivots, &mut work).unwrap();
        counts[usize::from(success)] += 1;
        let mut actual = vec![f64::from(success)];
        actual.extend(pivots.iter().map(|p| *p as f64));
        actual.extend(matrix.iter().chain(&rhs));
        assert_eq!(actual.len(), columns.len() - 1);
        for (i, (a, b)) in actual.iter().zip(&columns[1..]).enumerate() {
            let b: f64 = b.parse().unwrap();
            let tolerance = if i <= n {
                0.
            } else {
                1e-12 + 2e-12 * a.abs().max(b.abs())
            };
            assert!(
                (a - b).abs() <= tolerance,
                "{} field {i}: {a}, {b}",
                columns[0]
            );
        }
        if success {
            // A small componentwise backward error is meaningful even for the
            // scaled matrices whose solution has large forward sensitivity.
            for i in 0..n {
                let mut value = 0.;
                let mut scale = original_rhs[i].abs();
                for j in 0..n {
                    let coefficient = original[i.min(j) * n + i.max(j)];
                    value += coefficient * rhs[j];
                    scale += (coefficient * rhs[j]).abs();
                }
                assert!(
                    (value - original_rhs[i]).abs() <= 1e-12 + 5e-12 * scale,
                    "{} row {i}",
                    columns[0]
                );
            }
        } else {
            assert_eq!(rhs, original_rhs);
        }
        for i in 0..n {
            for j in 0..i {
                assert_eq!(matrix[i * n + j], original[i * n + j]);
            }
        }
    }
    assert_eq!(counts, [4, 222]);
}

#[test]
fn unused_triangle_is_ignored_and_invalid_workspace_is_rejected() {
    let mut a = [0., 1., f64::NAN, 0.];
    let mut b = [2., 3.];
    let mut pivots = [0; 2];
    let mut work = [0.; 2];
    assert!(solve_symmetric_in_place(&mut a, &mut b, &mut pivots, &mut work).unwrap());
    assert_eq!(b, [3., 2.]);
    assert_eq!(pivots, [-1, -1]);
    assert!(a[2].is_nan());
    assert!(solve_symmetric_in_place(&mut a, &mut b, &mut pivots, &mut work[..1]).is_err());
    a[0] = f64::NAN;
    assert!(solve_symmetric_in_place(&mut a, &mut b, &mut pivots, &mut work).is_err());
    assert!(solve_symmetric_in_place(&mut [], &mut [], &mut [], &mut []).is_err());
}
