//! Pure-argon PeTS Helmholtz model from pets.f90.
//! Temperature preparation shares diameter/packing terms across density sweeps.
//!
//! ```
//! use thermopack_rs::Pets;
//! # fn main() -> Result<(), thermopack_rs::Error> {
//! let model = Pets::new()?;
//! let prepared = model.prepare(300.0)?;
//! let state = prepared.caloric_tv(1e-3, 1.0)?; // total m3 and mol
//! assert!(state.pressure.value > 0.0);
//! let derivatives = prepared.residual_helmholtz(1e-3, 1.0)?;
//! assert!(derivatives.value.is_finite());
//! # Ok(())
//! # }
//! ```
use crate::{
    DensityHelmholtz, Error, IdealGas, PureHelmholtz, PurePressure, R, Result, component,
    error::positive,
};
use std::f64::consts::PI;
const NA: f64 = 6.02214076e23;
const A: [f64; 7] = [
    0.690603404,
    1.189317012,
    1.265604153,
    -24.34554201,
    93.67300357,
    -157.8773415,
    96.93736697,
];
const B: [f64; 7] = [
    0.664852128,
    2.10733079,
    -9.597951213,
    -17.37871193,
    30.17506222,
    209.3942909,
    -353.2743581,
];

#[path = "pets_equilibrium.rs"]
mod equilibrium;
pub use equilibrium::{PetsEquilibrium, PetsVolumeIteration, PureCriticalPoint};

#[path = "pets_flash.rs"]
mod flash;
pub use flash::PetsEquilibriumTv;

/// Pure-argon PeTS parameters and owned ideal-gas reference state.
#[derive(Debug, Clone)]
pub struct Pets {
    sigma: f64,
    epsilon_k: f64,
    ideal: IdealGas,
    molar_mass: f64,
}

/// Effective hard-sphere diameter (m), differentiated in temperature (K).
#[derive(Debug, Clone, Copy)]
pub struct PetsDiameter {
    pub value: f64,
    pub dt: f64,
    pub dtt: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct PetsContributions {
    pub total: DensityHelmholtz,
    pub hard_sphere: DensityHelmholtz,
    pub dispersion: DensityHelmholtz,
}

/// Extensive caloric quantities use total mole number and total volume.
/// Energies are J, entropy/Cv/Cp J/K, chemical potential J/mol.
#[derive(Debug, Clone, Copy)]
pub struct PetsTvProperties {
    pub pressure: PurePressure,
    pub helmholtz_energy: f64,
    pub internal_energy: f64,
    pub enthalpy: f64,
    pub entropy: f64,
    pub heat_capacity_v: f64,
    /// None for a singular/non-finite isobaric response. Phase stability is
    /// not tested; finite negative heat capacities remain observable.
    pub heat_capacity_p: Option<f64>,
    pub chemical_potential: f64,
    /// None when pressure is nonpositive.
    pub ln_fugacity_coefficient: Option<f64>,
}

/// Temperature-dependent terms borrowed from an immutable model.
#[derive(Debug, Clone, Copy)]
pub struct PreparedPets<'a> {
    model: &'a Pets,
    t: f64,
    diameter: PetsDiameter,
    packing: f64,
    packing_t: f64,
    packing_tt: f64,
    q: f64,
    dispersion_factor: f64,
}

impl Pets {
    /// Source defaults: sigma=3.42e-10 m, epsilon/k=136 K, argon ideal gas.
    pub fn new() -> Result<Self> {
        Self::with_parameters(3.42e-10, 136.)
    }
    pub fn with_parameters(sigma: f64, epsilon_k: f64) -> Result<Self> {
        positive(sigma, "PeTS segment diameter")?;
        positive(epsilon_k, "PeTS epsilon/k")?;
        let argon = component("AR")?;
        Ok(Self {
            sigma,
            epsilon_k,
            ideal: IdealGas::new(argon)?,
            molar_mass: argon.molar_mass * 1e-3,
        })
    }
    pub fn parameters(&self) -> (f64, f64) {
        (self.sigma, self.epsilon_k)
    }
    /// Matches the actual Fortran setter. The fork's Python setter is miswired.
    pub fn set_parameters(&mut self, sigma: f64, epsilon_k: f64) -> Result<()> {
        positive(sigma, "PeTS segment diameter")?;
        positive(epsilon_k, "PeTS epsilon/k")?;
        self.sigma = sigma;
        self.epsilon_k = epsilon_k;
        Ok(())
    }
    pub fn ideal_gas(&self) -> &IdealGas {
        &self.ideal
    }
    pub fn ideal_gas_mut(&mut self) -> &mut IdealGas {
        &mut self.ideal
    }
    pub fn prepare(&self, t: f64) -> Result<PreparedPets<'_>> {
        positive(t, "temperature")?;
        let b = 3.052785558 * self.epsilon_k;
        let inv = 1. / t;
        let inv2 = inv * inv;
        let b_inv2 = b * inv2;
        let expo = -0.127112544 * (-b * inv).exp();
        let diameter = PetsDiameter {
            value: self.sigma * (1. + expo),
            dt: self.sigma * expo * b_inv2,
            dtt: self.sigma * expo * (b_inv2 * b_inv2 - 2. * b * (inv2 * inv)),
        };
        let d2 = diameter.value * diameter.value;
        let packing = PI / 6. * (d2 * diameter.value) * NA;
        let packing_t = packing * 3. * diameter.dt / diameter.value;
        let packing_tt =
            packing / d2 * (6. * (diameter.dt * diameter.dt) + 3. * diameter.value * diameter.dtt);
        let q = self.epsilon_k / t;
        let dispersion_factor = -PI * NA * (self.sigma * self.sigma * self.sigma);
        if [
            diameter.value,
            diameter.dt,
            diameter.dtt,
            packing,
            packing_t,
            packing_tt,
            q,
            dispersion_factor,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || packing <= 0.
        {
            return Err(Error::InvalidInput("non-finite PeTS temperature terms"));
        }
        Ok(PreparedPets {
            model: self,
            t,
            diameter,
            packing,
            packing_t,
            packing_tt,
            q,
            dispersion_factor,
        })
    }
    pub fn diameter(&self, t: f64) -> Result<PetsDiameter> {
        Ok(self.prepare(t)?.diameter())
    }
    pub fn residual_density(&self, t: f64, rho: f64) -> Result<DensityHelmholtz> {
        self.prepare(t)?.residual_density(rho)
    }
    pub fn residual_helmholtz(&self, t: f64, volume: f64, moles: f64) -> Result<PureHelmholtz> {
        self.prepare(t)?.residual_helmholtz(volume, moles)
    }
    pub fn pressure_tv(&self, t: f64, volume: f64, moles: f64) -> Result<PurePressure> {
        self.prepare(t)?.pressure_tv(volume, moles)
    }
    pub fn caloric_tv(&self, t: f64, volume: f64, moles: f64) -> Result<PetsTvProperties> {
        self.prepare(t)?.caloric_tv(volume, moles)
    }
    /// B (m3/mol), C (m6/mol2) from the analytic zero-density limit.
    pub fn virial_coefficients(&self, t: f64) -> Result<(f64, f64)> {
        let f = self.residual_density(t, 0.)?;
        Ok((f.drho, f.drhorho))
    }
    /// This fork estimates virials at n=1 mol, V=1e7 m3. Its C expression
    /// subtracts nearly equal terms, so it has appreciable rounding noise.
    pub fn baseline_virial_coefficients(&self, t: f64) -> Result<(f64, f64)> {
        let v = 1e7_f64;
        let (ar, arr) = self.prepare(t)?.source_virial_density_derivatives(1. / v);
        let inv_v = 1. / v;
        let v2 = v * v;
        let v4 = v2 * v2;
        let fv = -(inv_v * inv_v) * ar;
        let fvv = 2. / (v2 * v) * ar + 1. / v4 * arr;
        let result = (-v2 * fv, v4 * (fvv + 2. * fv / v));
        if !result.0.is_finite() || !result.1.is_finite() {
            return Err(Error::InvalidInput("non-finite baseline virial estimate"));
        }
        Ok(result)
    }
    pub fn speed_of_sound_tv(&self, t: f64, v: f64, n: f64) -> Result<f64> {
        self.prepare(t)?.speed_of_sound_tv(v, n)
    }
    pub fn baseline_speed_of_sound_tv(&self, t: f64, v: f64, n: f64) -> Result<f64> {
        self.prepare(t)?.baseline_speed_of_sound_tv(v, n)
    }

    /// Pair potential divided by k (K), and radial force divided by k (K/m).
    /// At/above 2.5 sigma both are zero.
    pub fn potential(&self, distance: f64) -> Result<(f64, f64)> {
        positive(distance, "pair separation")?;
        if distance >= 2.5 * self.sigma {
            return Ok((0., 0.));
        }
        let r = self.sigma / distance;
        let r6 = crate::numerics::ipow(r, 6);
        let r12 = r6 * r6;
        let cut_inv = 1_f64 / 2.5;
        let cut = 4.
            * self.epsilon_k
            * (crate::numerics::ipow(cut_inv, 12) - crate::numerics::ipow(cut_inv, 6));
        let potential = 4. * self.epsilon_k * (r12 - r6) - cut;
        let force = 24. * self.epsilon_k / distance * (2. * r12 - r6);
        if !potential.is_finite() || !force.is_finite() {
            return Err(Error::InvalidInput("non-finite pair potential"));
        }
        Ok((potential, force))
    }
    /// Source vector convention: evaluation stops at the first cutoff entry.
    /// Later entries remain zero even if the distances are unsorted.
    pub fn potential_into(
        &self,
        distances: &[f64],
        values: &mut [f64],
        mut forces: Option<&mut [f64]>,
    ) -> Result<()> {
        if values.len() != distances.len()
            || forces.as_ref().is_some_and(|v| v.len() != distances.len())
            || distances.iter().any(|r| !r.is_finite() || *r <= 0.)
        {
            return Err(Error::InvalidInput(
                "positive separations and matching potential buffers required",
            ));
        }
        values.fill(0.);
        if let Some(f) = forces.as_deref_mut() {
            f.fill(0.);
        }
        for (i, &r) in distances.iter().enumerate() {
            if r >= 2.5 * self.sigma {
                break;
            }
            let (value, force) = self.potential(r)?;
            values[i] = value;
            if let Some(f) = forces.as_deref_mut() {
                f[i] = force;
            }
        }
        Ok(())
    }
    /// Reduced extensive A(full LJ)-A(truncated and shifted), including the
    /// two correction contributions and all first/second T,V,n derivatives.
    pub fn lj_truncation_correction(&self, t: f64, v: f64, n: f64) -> Result<PureHelmholtz> {
        positive(t, "temperature")?;
        positive(v, "total volume")?;
        positive(n, "mole number")?;
        let inv = 1_f64 / 2.5;
        let inv3 = inv * inv * inv;
        let inv9 = crate::numerics::ipow(inv, 9);
        let kc = 8. * PI * NA * (self.sigma * self.sigma * self.sigma) * self.epsilon_k;
        let a = kc * (n * n / v) / t * (inv9 / 9. - inv3 / 3.)
            + kc / 3. * (n * n / v) / t * (inv9 - inv3);
        let result = PureHelmholtz {
            value: a,
            dt: -a / t,
            dv: -a / v,
            dn: 2. * a / n,
            dtt: 2. * a / (t * t),
            dtv: a / (t * v),
            dtn: -2. * a / (t * n),
            dvv: 2. * a / (v * v),
            dvn: -2. * a / (v * n),
            dnn: 2. * a / (n * n),
        };
        if result.values().iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput("non-finite LJ truncation correction"));
        }
        Ok(result)
    }
}

impl PreparedPets<'_> {
    // Pressure-root work needs only density derivatives. Avoid computing and
    // validating the unused temperature/mixed derivatives at every grid node.
    // The public full Helmholtz kernel retains its complete derivative checks.
    fn pressure_density_derivatives(&self, rho: f64) -> (f64, f64, f64) {
        let e = rho * self.packing;
        let free = 1. - e;
        let free2 = free * free;
        let free4 = free2 * free2;
        let em1 = e - 1.;
        let em1_2 = em1 * em1;
        let hs_e = 2. * (e - 2.) / (em1_2 * em1);
        let hs_ee = 2. * (5. - 2. * e) / (em1_2 * em1_2);
        let (i1, i1e, i1ee) = polynomial(e, &A);
        let (i2, i2e, i2ee) = polynomial(e, &B);
        let o = (8. * e - 2. * e * e) / free4;
        let c = 1. / (1. + o);
        let oe = (-4. * e * e + 20. * e + 8.) / crate::numerics::ipow(free, 5);
        let oee = ((-8. * e + 20.) * free + 5. * (-4. * e * e + 20. * e + 8.))
            / crate::numerics::ipow(free, 6);
        let ce = -oe * c * c;
        let cee = 2. * ce * ce / c - oee * c * c;
        let d = c * i2;
        let de = ce * i2 + c * i2e;
        let dee = cee * i2 + 2. * ce * i2e + c * i2ee;
        let q = self.q;
        let w = 2. * i1 * q + d * q * q;
        let we = 2. * i1e * q + de * q * q;
        let wee = 2. * i1ee * q + dee * q * q;
        let wr = we * self.packing;
        let wrr = wee * self.packing * self.packing;
        let hs_r = hs_e * self.packing;
        let dispersion_r = self.dispersion_factor * (w + rho * wr);
        (
            hs_r + dispersion_r,
            hs_ee * self.packing * self.packing + self.dispersion_factor * (2. * wr + rho * wrr),
            hs_r.abs() + dispersion_r.abs(),
        )
    }
    // The source finite-volume third virial estimate loses many digits to
    // cancellation. Preserve its operation order in this compatibility path.
    fn source_virial_density_derivatives(&self, rho: f64) -> (f64, f64) {
        let e = rho * self.packing;
        let er = self.packing;
        let mut i1 = A[0];
        let mut i2 = B[0];
        let (mut i1e, mut i1ee, mut i2e, mut i2ee) = (0., 0., 0., 0.);
        for i in 1..=6 {
            let k = i as i32;
            let x = i as f64;
            i1 += A[i] * crate::numerics::ipow(e, k);
            i2 += B[i] * crate::numerics::ipow(e, k);
            i1e += A[i] * x * crate::numerics::ipow(e, k - 1);
            i2e += B[i] * x * crate::numerics::ipow(e, k - 1);
            i1ee += x * (x - 1.) * A[i] * crate::numerics::ipow(e, k - 2);
            i2ee += x * (x - 1.) * B[i] * crate::numerics::ipow(e, k - 2);
        }
        let er2 = er * er;
        let (i1r, i1rr, i2r, i2rr) = (i1e * er, i1ee * er2, i2e * er, i2ee * er2);
        let om = 1. - e;
        let om2 = om * om;
        let o = (8. * e - 2. * e * e) / (om2 * om2);
        let c = 1. / (1. + o);
        let oe = (-4. * e * e + 20. * e + 8.) / crate::numerics::ipow(om, 5);
        let ce = -oe * (c * c);
        let oee = ((-8. * e + 20.) * om + (-4. * e * e + 20. * e + 8.) * 5.)
            / crate::numerics::ipow(om, 6);
        let cee = 2. * (ce * ce) / c - oee * (c * c);
        let cr = ce * er;
        let crr = cee * er2;
        let m1 = self.model.epsilon_k / self.t
            * (self.model.sigma * self.model.sigma * self.model.sigma);
        let m2 = m1 * self.model.epsilon_k / self.t;
        let ar = (-2. * PI * (i1 + rho * i1r) * m1
            - PI * (c * i2 + rho * cr * i2 + rho * c * i2r) * m2)
            * NA;
        let arr = (-2. * PI * (2. * i1r + rho * i1rr) * m1
            - PI * (2. * cr * i2
                + 2. * c * i2r
                + 2. * rho * cr * i2r
                + rho * crr * i2
                + rho * c * i2rr)
                * m2)
            * NA;
        let em1 = e - 1.;
        let em1_2 = em1 * em1;
        let hs_r = 2. * (e - 2.) / (em1_2 * em1) * er;
        let hs_rr = 2. * (5. - 2. * e) / (em1_2 * em1_2) * er2;
        (hs_r + ar, hs_rr + arr)
    }
    pub fn diameter(&self) -> PetsDiameter {
        self.diameter
    }
    pub fn residual_density(&self, rho: f64) -> Result<DensityHelmholtz> {
        Ok(self.contributions(rho)?.total)
    }
    pub fn contributions(&self, rho: f64) -> Result<PetsContributions> {
        if !rho.is_finite() || rho < 0. {
            return Err(Error::InvalidInput("nonnegative molar density required"));
        }
        let e = rho * self.packing;
        if !(0. ..1.).contains(&e) {
            return Err(Error::InvalidInput(
                "PeTS packing fraction must be below one",
            ));
        }
        let et = rho * self.packing_t;
        let ett = rho * self.packing_tt;
        let free = 1. - e;
        let free2 = free * free;
        let free4 = free2 * free2;
        let em1 = e - 1.;
        let em1_2 = em1 * em1;
        let hs = (4. * e - 3. * e * e) / free2;
        let hs_e = 2. * (e - 2.) / (em1_2 * em1);
        let hs_ee = 2. * (5. - 2. * e) / (em1_2 * em1_2);
        let hard_sphere = DensityHelmholtz {
            value: hs,
            dt: hs_e * et,
            drho: hs_e * self.packing,
            dtt: hs_e * ett + hs_ee * et * et,
            dtrho: hs_e * self.packing_t + hs_ee * et * self.packing,
            drhorho: hs_ee * self.packing * self.packing,
        };
        let (i1, i1e, i1ee) = polynomial(e, &A);
        let (i2, i2e, i2ee) = polynomial(e, &B);
        let o = (8. * e - 2. * e * e) / free4;
        let c = 1. / (1. + o);
        let oe = (-4. * e * e + 20. * e + 8.) / crate::numerics::ipow(free, 5);
        let oee = ((-8. * e + 20.) * free + 5. * (-4. * e * e + 20. * e + 8.))
            / crate::numerics::ipow(free, 6);
        let ce = -oe * c * c;
        let cee = 2. * ce * ce / c - oee * c * c;
        let d = c * i2;
        let de = ce * i2 + c * i2e;
        let dee = cee * i2 + 2. * ce * i2e + c * i2ee;
        let q = self.q;
        let qt = -q / self.t;
        let qtt = 2. * q / (self.t * self.t);
        let w = 2. * i1 * q + d * q * q;
        let we = 2. * i1e * q + de * q * q;
        let wee = 2. * i1ee * q + dee * q * q;
        let wq = 2. * i1 + 2. * d * q;
        let weq = 2. * i1e + 2. * de * q;
        let wr = we * self.packing;
        let wt = we * et + wq * qt;
        let wrr = wee * self.packing * self.packing;
        let wrt = wee * et * self.packing + we * self.packing_t + weq * self.packing * qt;
        let wtt = wee * et * et + we * ett + 2. * weq * et * qt + 2. * d * qt * qt + wq * qtt;
        let factor = self.dispersion_factor;
        let dispersion = DensityHelmholtz {
            value: factor * rho * w,
            dt: factor * rho * wt,
            drho: factor * (w + rho * wr),
            dtt: factor * rho * wtt,
            dtrho: factor * (wt + rho * wrt),
            drhorho: factor * (2. * wr + rho * wrr),
        };
        let total = hard_sphere.add(dispersion);
        if total
            .values()
            .iter()
            .chain(hard_sphere.values().iter())
            .chain(dispersion.values().iter())
            .any(|v| !v.is_finite())
        {
            return Err(Error::InvalidInput("non-finite PeTS Helmholtz derivatives"));
        }
        Ok(PetsContributions {
            total,
            hard_sphere,
            dispersion,
        })
    }
    pub fn residual_helmholtz(&self, volume: f64, moles: f64) -> Result<PureHelmholtz> {
        positive(volume, "total volume")?;
        positive(moles, "mole number")?;
        self.residual_density(moles / volume)?
            .extensive(volume, moles)
    }
    pub fn pressure_tv(&self, v: f64, n: f64) -> Result<PurePressure> {
        let f = self.residual_helmholtz(v, n)?;
        self.pressure_from(f, v, n)
    }
    fn pressure_from(&self, f: PureHelmholtz, v: f64, n: f64) -> Result<PurePressure> {
        let p = PurePressure {
            value: R * self.t * (n / v - f.dv),
            dt: R * (n / v - f.dv) - R * self.t * f.dtv,
            dv: -R * self.t * (n / (v * v) + f.dvv),
            dn: R * self.t * (1. / v - f.dvn),
        };
        if [p.value, p.dt, p.dv, p.dn].iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput("non-finite PeTS pressure"));
        }
        Ok(p)
    }
    pub fn caloric_tv(&self, v: f64, n: f64) -> Result<PetsTvProperties> {
        let f = self.residual_helmholtz(v, n)?;
        let p = self.pressure_from(f, v, n)?;
        let t = self.t;
        let ideal = self.model.ideal.properties(t, n * R * t / v)?;
        let u = n * (ideal.enthalpy - R * t) - R * t * t * f.dt;
        let s = n * ideal.entropy - R * (f.value + t * f.dt);
        let cv = n * (ideal.enthalpy_dt - R) - R * (2. * t * f.dt + t * t * f.dtt);
        let a = u - t * s;
        let h = u + p.value * v;
        let mu = ideal.enthalpy - t * ideal.entropy + R * t * f.dn;
        let cp = cv - t * p.dt * p.dt / p.dv;
        if [u, s, cv, a, h, mu].iter().any(|v| !v.is_finite()) {
            return Err(Error::InvalidInput("non-finite PeTS caloric state"));
        }
        let phi = if p.value > 0. {
            Some(f.dn - (p.value * v / (n * R * t)).ln())
        } else {
            None
        };
        Ok(PetsTvProperties {
            pressure: p,
            helmholtz_energy: a,
            internal_energy: u,
            enthalpy: h,
            entropy: s,
            heat_capacity_v: cv,
            heat_capacity_p: if cp.is_finite() { Some(cp) } else { None },
            chemical_potential: mu,
            ln_fugacity_coefficient: phi,
        })
    }
    fn isentropic_stiffness(&self, v: f64, n: f64) -> Result<f64> {
        let f = self.residual_helmholtz(v, n)?;
        let p = self.pressure_from(f, v, n)?;
        let cp = self.model.ideal.properties(self.t, 1e5)?.enthalpy_dt;
        let cv = n * (cp - R) - R * (2. * self.t * f.dt + self.t * self.t * f.dtt);
        Ok(p.dv - self.t * p.dt * p.dt / cv)
    }
    /// Physical single-phase sound speed (m/s). Invalid stiffness is an error.
    pub fn speed_of_sound_tv(&self, v: f64, n: f64) -> Result<f64> {
        let squared = -v * v / (n * self.model.molar_mass) * self.isentropic_stiffness(v, n)?;
        positive(squared, "sound speed squared")?;
        Ok(squared.sqrt())
    }
    /// Source convention: nonnegative/NaN isentropic stiffness returns zero.
    pub fn baseline_speed_of_sound_tv(&self, v: f64, n: f64) -> Result<f64> {
        let stiffness = self.isentropic_stiffness(v, n)?;
        if stiffness < 0. {
            let squared = -v * v / (n * self.model.molar_mass) * stiffness;
            positive(squared, "sound speed squared")?;
            Ok(squared.sqrt())
        } else {
            Ok(0.)
        }
    }
}

// Simultaneous Horner evaluation of a polynomial and its two derivatives.
fn polynomial(e: f64, c: &[f64; 7]) -> (f64, f64, f64) {
    let (mut p, mut dp, mut ddp) = (c[6], 0., 0.);
    for &coefficient in c[..6].iter().rev() {
        ddp = ddp * e + 2. * dp;
        dp = dp * e + p;
        p = p * e + coefficient;
    }
    (p, dp, ddp)
}
