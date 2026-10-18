// SPDX-License-Identifier: BSD-3-Clause
// Derived from LAPACK DSYTF2, DSYTRS2, DSYCONV and reference BLAS DTRSM/DSYR.
// Copyright and redistribution terms: rust/LAPACK-LICENSE.

/// Solve a symmetric system using upper-triangle Bunch–Kaufman factorization.
///
/// The matrix is row-major; only its upper triangle is read or overwritten.
/// `pivots` must have exactly n entries and `work` at least n. Both are reusable
/// caller-owned buffers, so the solve allocates no memory. Pivots retain the
/// source's signed one-based convention, including paired negative 2x2 blocks.
///
/// True means the solution overwrote `rhs`. False means a zero pivot was found:
/// `rhs` is unchanged, while the factorization and pivots are available. Invalid
/// inputs and non-finite computed results return errors. The factorization is
/// unblocked at every size; LAPACK's workspace-dependent blocked path is not
/// reproduced. This routine does not estimate the condition number.
pub fn solve_symmetric_in_place(
    a: &mut [f64],
    rhs: &mut [f64],
    pivots: &mut [isize],
    work: &mut [f64],
) -> crate::Result<bool> {
    let n = rhs.len();
    if n == 0
        || n.checked_mul(n) != Some(a.len())
        || pivots.len() != n
        || work.len() < n
        || rhs.iter().any(|x| !x.is_finite())
        || (0..n).any(|i| (i..n).any(|j| !a[i * n + j].is_finite()))
    {
        return Err(crate::Error::InvalidInput(
            "finite upper symmetric system and matching pivot/work buffers required",
        ));
    }
    let alpha = (1. + 17.0_f64.sqrt()) / 8.;
    let mut remaining = n;
    let mut singular = false;
    while remaining > 0 {
        let k = remaining - 1;
        let mut step = 1;
        let mut imax = 0;
        let mut colmax = 0.;
        for i in 0..k {
            if a[i * n + k].abs() > colmax {
                colmax = a[i * n + k].abs();
                imax = i;
            }
        }
        let diagonal = a[k * n + k].abs();
        let mut pivot = k;
        if diagonal.max(colmax) == 0. {
            singular = true;
        } else {
            if diagonal < alpha * colmax {
                let mut rowmax = 0.0_f64;
                for j in imax + 1..=k {
                    rowmax = rowmax.max(a[imax * n + j].abs());
                }
                for i in 0..imax {
                    rowmax = rowmax.max(a[i * n + imax].abs());
                }
                if diagonal < alpha * colmax * (colmax / rowmax) {
                    pivot = imax;
                    if a[imax * n + imax].abs() < alpha * rowmax {
                        step = 2;
                    }
                }
            }
            let kk = k + 1 - step;
            if pivot != kk {
                for i in 0..pivot {
                    a.swap(i * n + kk, i * n + pivot);
                }
                for i in pivot + 1..kk {
                    a.swap(i * n + kk, pivot * n + i);
                }
                a.swap(kk * n + kk, pivot * n + pivot);
                if step == 2 {
                    a.swap((k - 1) * n + k, pivot * n + k);
                }
            }
            if step == 1 {
                let inverse = 1. / a[k * n + k];
                for j in 0..k {
                    if a[j * n + k] != 0. {
                        let temp = -inverse * a[j * n + k];
                        for i in 0..=j {
                            a[i * n + j] += a[i * n + k] * temp;
                        }
                    }
                }
                for i in 0..k {
                    a[i * n + k] *= inverse;
                }
            } else if k > 1 {
                let d12 = a[(k - 1) * n + k];
                let d22 = a[(k - 1) * n + k - 1] / d12;
                let d11 = a[k * n + k] / d12;
                let d12 = (1. / (d11 * d22 - 1.)) / d12;
                for j in (0..k - 1).rev() {
                    let previous = d12 * (d11 * a[j * n + k - 1] - a[j * n + k]);
                    let current = d12 * (d22 * a[j * n + k] - a[j * n + k - 1]);
                    for i in (0..=j).rev() {
                        a[i * n + j] =
                            a[i * n + j] - a[i * n + k] * current - a[i * n + k - 1] * previous;
                    }
                    a[j * n + k] = current;
                    a[j * n + k - 1] = previous;
                }
            }
        }
        let signed = (pivot + 1) as isize;
        pivots[k] = if step == 1 { signed } else { -signed };
        if step == 2 {
            pivots[k - 1] = -signed;
        }
        remaining -= step;
    }
    if (0..n).any(|i| (i..n).any(|j| !a[i * n + j].is_finite())) {
        return Err(crate::Error::NotConverged(
            "non-finite symmetric factorization",
        ));
    }
    if singular {
        return Ok(false);
    }
    // DSYCONV: extract off-diagonal blocks and apply trailing permutations.
    work[..n].fill(0.);
    remaining = n;
    while remaining > 1 {
        let k = remaining - 1;
        if pivots[k] < 0 {
            work[k] = a[(k - 1) * n + k];
            a[(k - 1) * n + k] = 0.;
            remaining -= 2;
        } else {
            remaining -= 1;
        }
    }
    remaining = n;
    while remaining > 0 {
        let k = remaining - 1;
        let pivot = pivots[k].unsigned_abs() - 1;
        let row = if pivots[k] > 0 { k } else { k - 1 };
        for j in k + 1..n {
            a.swap(pivot * n + j, row * n + j);
        }
        remaining -= if pivots[k] > 0 { 1 } else { 2 };
    }
    // P^T B, followed by the reference upper-unit triangular solve.
    remaining = n;
    while remaining > 0 {
        let k = remaining - 1;
        let step = if pivots[k] > 0 { 1 } else { 2 };
        rhs.swap(k + 1 - step, pivots[k].unsigned_abs() - 1);
        remaining -= step;
    }
    for k in (0..n).rev() {
        if rhs[k] != 0. {
            for i in 0..k {
                rhs[i] -= rhs[k] * a[i * n + k];
            }
        }
    }
    remaining = n;
    while remaining > 0 {
        let k = remaining - 1;
        if pivots[k] > 0 {
            rhs[k] *= 1. / a[k * n + k];
            remaining -= 1;
        } else {
            let off = work[k];
            let previous = a[(k - 1) * n + k - 1] / off;
            let current = a[k * n + k] / off;
            let denominator = previous * current - 1.;
            let b_previous = rhs[k - 1] / off;
            let b_current = rhs[k] / off;
            rhs[k - 1] = (current * b_previous - b_current) / denominator;
            rhs[k] = (previous * b_current - b_previous) / denominator;
            remaining -= 2;
        }
    }
    for i in 0..n {
        let mut temp = rhs[i];
        for k in 0..i {
            temp -= a[k * n + i] * rhs[k];
        }
        rhs[i] = temp;
    }
    let mut k = 0;
    while k < n {
        rhs.swap(k, pivots[k].unsigned_abs() - 1);
        k += if pivots[k] > 0 { 1 } else { 2 };
    }
    // Restore the factorization returned by DSYTF2, as DSYTRS2 does.
    k = 0;
    while k < n {
        let pivot = pivots[k].unsigned_abs() - 1;
        let last = if pivots[k] > 0 { k } else { k + 1 };
        for j in last + 1..n {
            a.swap(pivot * n + j, k * n + j);
        }
        k = last + 1;
    }
    remaining = n;
    while remaining > 1 {
        k = remaining - 1;
        if pivots[k] < 0 {
            a[(k - 1) * n + k] = work[k];
            remaining -= 2;
        } else {
            remaining -= 1;
        }
    }
    if rhs.iter().any(|x| !x.is_finite()) {
        return Err(crate::Error::NotConverged("non-finite symmetric solution"));
    }
    Ok(true)
}
