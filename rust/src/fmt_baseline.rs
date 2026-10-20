//! Source arithmetic trace for the FMT compatibility API.
//! Retains hyperdual operand order, integer powers and inverse derivative floor.
use super::{FmtDensity, FmtModel};
use crate::{Error, Result};
use std::f64::consts::PI;

fn reciprocal(hi: f64, lo: f64) -> f64 {
    let r = 1. / hi;
    r + ((-r).mul_add(hi, 1.) - r * lo) / hi
}

#[derive(Clone, Copy, Default)]
struct D {
    v: f64,
    x: f64,
    y: f64,
    xy: f64,
}
impl D {
    fn c(v: f64) -> Self {
        Self {
            v,
            ..Self::default()
        }
    }
    fn scale(self, k: f64) -> Self {
        Self {
            v: self.v * k,
            x: self.x * k,
            y: self.y * k,
            xy: self.xy * k,
        }
    }
    fn pow(self, n: usize) -> Self {
        let mut result = Self::c(1.);
        for _ in 0..n {
            result = result * self;
        }
        result
    }
    fn inverse(self) -> Self {
        let z = if self.v.abs() < 1e-15 {
            if self.v >= 0. { 1e-15_f64 } else { -1e-15_f64 }
        } else {
            self.v
        };
        // Correct the rounded square/cube before reciprocal division. The
        // Windows Rust math runtime's pow(x,-2/-3) loses a bit here relative
        // to the source runtime; White Bear cancellation amplifies it.
        let square = z * z;
        let square_lo = z.mul_add(z, -square);
        let cube = square * z;
        let cube_lo = square.mul_add(z, -cube) + square_lo * z;
        let first = -reciprocal(square, square_lo);
        let second = 2. * reciprocal(cube, cube_lo);
        Self {
            v: 1. / self.v,
            x: self.x * first,
            y: self.y * first,
            xy: self.xy * first + self.x * self.y * second,
        }
    }
    fn log(self) -> Self {
        let inv = 1. / self.v;
        let x = self.x * inv;
        let y = self.y * inv;
        Self {
            v: self.v.ln(),
            x,
            y,
            xy: self.xy * inv - x * y,
        }
    }
}
impl std::ops::Add for D {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            v: self.v + b.v,
            x: self.x + b.x,
            y: self.y + b.y,
            xy: self.xy + b.xy,
        }
    }
}
impl std::ops::Sub for D {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        Self {
            v: self.v - b.v,
            x: self.x - b.x,
            y: self.y - b.y,
            xy: self.xy - b.xy,
        }
    }
}
impl std::ops::Mul for D {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self {
            v: self.v * b.v,
            x: self.v * b.x + self.x * b.v,
            y: self.v * b.y + self.y * b.v,
            xy: self.v * b.xy + self.x * b.y + self.y * b.x + self.xy * b.v,
        }
    }
}

fn function(model: FmtModel, n: [D; 6]) -> D {
    let raw = D::c(1.) - n[3];
    let u = if raw.v > 1e-12 { raw } else { D::c(1e-12) };
    let log = u.log();
    let first = n[0].scale(-1.) * log;
    let s = n[1] * n[2] - n[4] * n[5];
    let q = n[2].pow(3) - n[2].scale(3.) * n[5].pow(2);
    match model {
        FmtModel::Rosenfeld => first + s * u.inverse() + q * u.pow(2).scale(24. * PI).inverse(),
        FmtModel::WhiteBear => {
            first
                + s * u.inverse()
                + (q * (n[3] + u.pow(2) * log))
                    * ((n[3].pow(2).scale(36. * PI)) * u.pow(2)).inverse()
        }
        FmtModel::WhiteBearII => {
            let eta2 = n[3].pow(2);
            let phi2 = (D::c(2.) - n[3] + ((u.scale(2.) * log) * n[3].inverse())).scale(1. / 3.);
            let phi3 = (n[3].scale(2.) - eta2.scale(3.)
                + n[3].scale(2.) * eta2
                + u.pow(2).scale(2.) * log)
                .scale(1. / 3.)
                * eta2.inverse();
            first
                + (s * (D::c(1.) + phi2)) * u.inverse()
                + (q * (D::c(1.) - phi3)) * u.pow(2).scale(24. * PI).inverse()
        }
    }
}

pub(super) fn evaluate(model: FmtModel, n: [f64; 6]) -> Result<FmtDensity> {
    if n.iter().any(|v| !v.is_finite()) {
        return Err(Error::InvalidInput(
            "finite FMT weighted densities required",
        ));
    }
    let mut variables = n.map(D::c);
    let mut result = FmtDensity {
        value: 0.,
        gradient: [0.; 6],
        hessian: [[0.; 6]; 6],
    };
    for i in 0..6 {
        variables[i].x = 1.;
        for j in i..6 {
            variables[j].y = 1.;
            let q = function(model, variables);
            result.hessian[i][j] = q.xy;
            result.hessian[j][i] = q.xy;
            result.gradient[i] = q.x;
            result.value = q.v;
            variables[j].y = 0.;
        }
        variables[i].x = 0.;
    }
    if !result.value.is_finite()
        || result
            .gradient
            .iter()
            .chain(result.hessian.iter().flatten())
            .any(|v| !v.is_finite())
    {
        return Err(Error::InvalidInput("non-finite FMT energy or derivatives"));
    }
    Ok(result)
}
