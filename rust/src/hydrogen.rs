//! Leachman/Valenti ideal-H2 correlations, ported from idealh2.f90.
//! The legacy H/S primitives have a factor 1000 absent from Cp. Keep that
//! behavior explicit; ideal::IdealProperties also reports the actual dH/dT.
use crate::{Error, Result};

const U: [[f64; 7]; 4] = [
    [
        4.30256, 13.0289, -47.7365, 50.0013, -18.6261, 0.993973, 0.53678,
    ],
    [1.616, -0.4117, -0.792, 0.758, 1.217, 0.0, 0.0],
    [2.54151, -2.3661, 1.00365, 1.22447, 0.0, 0.0, 0.0],
    [
        29.5750, -83.1270, 54.4050, 1.84050, -1.80600, 2.69390, -2.09290,
    ],
];
const V: [[f64; 7]; 4] = [
    [499.0, 826.5, 970.8, 1166.2, 1341.4, 5395.0, 10185.0],
    [531.0, 751.0, 1989.0, 2484.0, 6859.0, 0.0, 0.0],
    [856.0, 1444.0, 2194.0, 6968.0, 0.0, 0.0, 0.0],
    [201.9, 263.5, 294.7, 1046.3, 1241.5, 5394.8, 10185.0],
];
const COUNT: [usize; 4] = [7, 5, 4, 7];
const H_REF: [f64; 4] = [-1413.6, -354.4, 0.0, -1413.6];
const S_REF: [f64; 4] = [0.0, 13.7009, 18.2678, 0.0];
// Preserve the fitting gas constant, which differs from the library R.
const FIT_R: f64 = 8.3144772;

pub(crate) fn primitives(ident: &str, t: f64) -> Result<(f64, f64, f64)> {
    let i = match ident {
        "P-H2" => 0,
        "N-H2" => 1,
        "O-H2" => 2,
        "E-H2" => 3,
        _ => {
            return Err(Error::InvalidInput(
                "hydrogen correlation requires an H2 isomer",
            ));
        }
    };
    let count = COUNT[i] - if t <= 30.0 { 2 } else { 0 };
    let (mut cp, mut h, mut s) = (2.5, 2.5 * t, 2.5 * t.ln());
    for k in 0..count {
        let x = V[i][k] / t;
        let q = (-x).exp();
        let d = -(-x).exp_m1();
        cp += U[i][k] * x * x * q / (d * d);
        h += U[i][k] * V[i][k] * q / d;
        s += U[i][k] * (x * q / d - d.ln());
    }
    Ok((
        cp * FIT_R,
        H_REF[i] + 1000.0 * FIT_R * h,
        S_REF[i] + 1000.0 * FIT_R * s,
    ))
}
