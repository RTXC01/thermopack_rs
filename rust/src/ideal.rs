//! Analytic ideal-gas caloric correlations from ideal.f90.
use crate::error::positive;
use crate::{Component, Error, R, Result};

#[derive(Debug, Clone, Copy)]
pub struct IdealProperties {
    pub cp: f64,
    pub enthalpy: f64,
    pub entropy: f64,
}

impl Component {
    /// Pure ideal gas at T [K], P [Pa], with the fork's reference state.
    pub fn ideal(&self, temperature: f64, pressure: f64) -> Result<IdealProperties> {
        positive(temperature, "temperature")?;
        positive(pressure, "pressure")?;
        let (cp, h, s) = self.ideal_primitives(temperature)?;
        let valid_reference = (self.reference_enthalpy != 0.0 || self.reference_entropy != 0.0)
            && self.reference_enthalpy.abs() < 1e15
            && self.reference_entropy.abs() < 1e15;
        let (enthalpy, entropy) = if valid_reference {
            let (_, h0, s0) = self.ideal_primitives(298.15)?;
            (
                self.reference_enthalpy + h - h0,
                self.reference_entropy + s - s0 - R * (pressure / self.reference_pressure).ln(),
            )
        } else {
            (h, s - R * pressure.ln())
        };
        Ok(IdealProperties {
            cp,
            enthalpy,
            entropy,
        })
    }

    fn ideal_primitives(&self, t: f64) -> Result<(f64, f64, f64)> {
        let c = &self.cp;
        match self.cp_type {
            1 | 2 | 3 | 4 | 5 | 6 | 8 | 11 => {
                let mut p = [0.0; 5];
                let mut inv = 0.0;
                let (scale, ts) = match self.cp_type {
                    1 => {
                        p[..4].copy_from_slice(&c[..4]);
                        (4.1868, t)
                    }
                    2 => {
                        for i in 0..5 {
                            p[i] = (i + 1) as f64 * c[i + 1];
                        }
                        (self.molar_mass, t)
                    }
                    3 => {
                        p[..3].copy_from_slice(&c[..3]);
                        (4.1868 * self.molar_mass, t * 1.8)
                    }
                    4 => {
                        p[..4].copy_from_slice(&c[..4]);
                        (1.0, t)
                    }
                    5 => {
                        p[..4].copy_from_slice(&c[..4]);
                        inv = c[4];
                        (self.molar_mass, t)
                    }
                    6 => {
                        p.copy_from_slice(&c[..5]);
                        (self.molar_mass, t)
                    }
                    8 => {
                        p.copy_from_slice(&c[..5]);
                        (R, t)
                    }
                    11 => {
                        p[..4].copy_from_slice(&c[..4]);
                        inv = c[4];
                        (1.0, t * 1e-3)
                    }
                    _ => unreachable!(),
                };
                let cp = p.iter().rev().fold(0.0, |acc, ci| acc * ts + ci) + inv / (ts * ts);
                let mut h = -inv / ts;
                let mut s = p[0] * ts.ln() - inv / (2.0 * ts * ts);
                for (i, pi) in p.iter().enumerate() {
                    h += pi * ts.powi((i + 1) as i32) / (i + 1) as f64;
                    if i > 0 {
                        s += pi * ts.powi(i as i32) / i as f64;
                    }
                }
                let hscale = match self.cp_type {
                    3 => 1.0 / 1.8,
                    11 => 1e3,
                    _ => 1.0,
                };
                Ok((scale * cp, scale * h * hscale, scale * s))
            }
            7 => {
                let x = c[2] / t;
                let y = c[4] / t;
                let sinh_ratio = if x.abs() < 1e-8 { 1.0 } else { x / x.sinh() };
                let cosh_ratio = y / y.cosh();
                let cp = c[0] + c[1] * sinh_ratio.powi(2) + c[3] * cosh_ratio.powi(2);
                let h = c[0] * t
                    + if x == 0.0 {
                        c[1] * t
                    } else {
                        c[1] * c[2] / x.tanh()
                    }
                    - c[3] * c[4] * y.tanh();
                let s = c[0] * t.ln()
                    + c[1] * (x / x.tanh() - log_sinh(x))
                    + c[3] * (log_cosh(y) - y * y.tanh());
                Ok((cp * 1e-3, h * 1e-3, s * 1e-3))
            }
            9 => {
                positive(t + c[3], "Mogensen logarithm domain")?;
                let cp = c[0] + c[1] * t + c[2] / (t + c[3]);
                let h = c[0] * t + c[1] * t * t / 2.0 + c[2] * (t + c[3]).ln();
                let s = c[0] * t.ln()
                    + c[1] * t
                    + if c[3] == 0.0 {
                        -c[2] / t
                    } else {
                        c[2] / c[3] * (t.ln() - (t + c[3]).ln())
                    };
                Ok((cp, h, s))
            }
            12 => {
                let (mut cp, mut h, mut s) = (c[0], c[0] * t, c[0] * t.ln());
                for pair in c[1..].as_chunks::<2>().0 {
                    let (amplitude, theta) = (pair[0], pair[1]);
                    if amplitude.abs() < 1e-12 {
                        break;
                    }
                    let x = theta / t;
                    let e = (-x).exp();
                    let denom = -(-x).exp_m1();
                    cp += amplitude * x * x * e / (denom * denom);
                    h += amplitude * theta * e / denom;
                    s += amplitude * (x * e / denom - denom.ln());
                }
                Ok((cp * R, h * R, s * R))
            }
            10 => Err(Error::Unsupported("hydrogen ideal heat capacity")),
            _ => Err(Error::Unsupported(
                "component has no default ideal heat-capacity correlation",
            )),
        }
    }
}

fn log_sinh(x: f64) -> f64 {
    x.abs() + (-(-2.0 * x.abs()).exp()).ln_1p() - 2.0_f64.ln()
}
fn log_cosh(x: f64) -> f64 {
    x.abs() + (-2.0 * x.abs()).exp().ln_1p() - 2.0_f64.ln()
}
