//! Bounded pure MBWR pressure branches and analytic critical points.
//!
//! ```
//! use thermopack_rs::{Mbwr, MbwrKind, Phase};
//! # fn main() -> Result<(), thermopack_rs::Error> {
//! let model = Mbwr::new(MbwrKind::Mbwr32, "C1")?;
//! let critical = model.critical_point()?;
//! assert!((critical.temperature - 190.56).abs() < 0.1);
//! let branches = model.prepare_equilibrium(150., 30000.)?;
//! let vapor_volume = branches.volume(1e6, Phase::Vapor)?;
//! assert!((model.pressure_tv(150., vapor_volume, 1.)?.value - 1e6).abs() < 0.01);
//! # Ok(())
//! # }
//! ```
use super::*;
use crate::{Phase, PureCriticalPoint};

/// Prepared stable pressure branches in an explicitly bounded density domain.
/// The upper limit is a search boundary, not a fitted EOS validity claim.
#[derive(Debug, Clone)]
pub struct MbwrEquilibrium<'a> {
    prepared: PreparedMbwr<'a>,
    maximum_density: f64,
    spinodals: Vec<f64>,
    intervals: Vec<Interval>,
}
#[derive(Debug, Clone, Copy)]
struct Interval {
    low: f64,
    high: f64,
    low_pressure: f64,
    high_pressure: f64,
}

impl Mbwr {
    /// Original tabulated approximate critical properties. Alternate fits with
    /// the same component identifier inherit the same source table entry.
    /// None means no explicit entry, rather than a fabricated critical point.
    pub fn source_critical_point(&self) -> Option<PureCriticalPoint> {
        let record = self.record?;
        let family = self.kind.coefficient_count() as u8;
        data::STORED_PROPERTIES.iter().find_map(|&(kind, id, a)| {
            (kind == family && id == record.identifier).then(|| PureCriticalPoint {
                temperature: a[0],
                pressure: a[1],
                volume: 1e-3 / a[2],
                compressibility: a[1] / (1000. * a[2] * R * a[0]),
            })
        })
    }
    /// Original triple-point (K, Pa) metadata, where explicitly stored.
    pub fn source_triple_point(&self) -> Option<(f64, f64)> {
        let record = self.record?;
        data::STORED_PROPERTIES.iter().find_map(|&(kind, id, a)| {
            (kind == 0 && id == record.identifier).then_some((a[0], a[1]))
        })
    }
    /// Solve P_rho=P_rhorho=0 from source-table or component initial values.
    /// These local equations do not select among multiple critical points of
    /// an extrapolated or custom EOS. Use explicit guesses for another basin.
    pub fn critical_point(&self) -> Result<PureCriticalPoint> {
        let (t, rho) = if let Some(c) = self.source_critical_point() {
            (c.temperature, 1. / c.volume)
        } else {
            let record = self.record.ok_or(Error::InvalidInput(
                "custom MBWR critical point requires explicit initial T and density",
            ))?;
            let seed = data::COMPONENT_CRITICAL_SEEDS
                .iter()
                .find(|(identifier, _)| *identifier == record.identifier)
                .ok_or_else(|| Error::UnknownComponent(record.identifier.to_owned()))?
                .1;
            (seed[0], seed[1])
        };
        self.critical_point_from(t, rho)
    }
    /// Damped Newton in log(T), log(rho), with analytic pressure derivatives
    /// through third density order and second density/first temperature order.
    /// Initial density is SI mol/m3. Positive P and ordinary P_rhorhorho>0
    /// are required at the returned local vapor/liquid critical point.
    pub fn critical_point_from(&self, mut t: f64, mut rho: f64) -> Result<PureCriticalPoint> {
        positive(t, "critical initial temperature")?;
        positive(rho, "critical initial molar density")?;
        for _ in 0..80 {
            let prepared = self.prepare(t)?;
            let q = prepared.critical_jet(rho)?;
            let norm = (q[1] / (R * t)).hypot(rho * q[2] / (R * t));
            if norm < 2e-11 && q[0] > 0. && q[3] > 0. {
                let result = PureCriticalPoint {
                    temperature: t,
                    pressure: q[0],
                    volume: 1. / rho,
                    compressibility: q[0] / (rho * R * t),
                };
                if result.volume.is_finite() && result.compressibility.is_finite() {
                    return Ok(result);
                }
                break;
            }
            let (a, b, c, d) = (t * q[4], rho * q[2], t * q[5], rho * q[3]);
            let determinant = a * d - b * c;
            if !determinant.is_finite() || determinant == 0. {
                break;
            }
            let st = (-q[1] * d + b * q[2]) / determinant;
            let sr = (-a * q[2] + q[1] * c) / determinant;
            let mut length = 1_f64.min(0.4 / st.abs().max(sr.abs()));
            let mut accepted = false;
            for _ in 0..24 {
                let nt = t * (length * st).exp();
                let nr = rho * (length * sr).exp();
                if let Ok(nq) = self.prepare(nt).and_then(|p| p.critical_jet(nr)) {
                    let nn = (nq[1] / (R * nt)).hypot(nr * nq[2] / (R * nt));
                    if nn < norm {
                        t = nt;
                        rho = nr;
                        accepted = true;
                        break;
                    }
                }
                length *= 0.5;
            }
            if !accepted {
                break;
            }
        }
        Err(Error::NotConverged("MBWR local critical point"))
    }
    pub fn prepare_equilibrium(&self, t: f64, maximum_density: f64) -> Result<MbwrEquilibrium<'_>> {
        self.prepare(t)?.equilibrium(maximum_density)
    }
}

// Simultaneous value and derivatives through third order; ascending powers.
fn jet(coefficients: &[f64], x: f64) -> [f64; 4] {
    let mut q = [0.; 4];
    for &c in coefficients.iter().rev() {
        q[3] = q[3] * x + 3. * q[2];
        q[2] = q[2] * x + 2. * q[1];
        q[1] = q[1] * x + q[0];
        q[0] = q[0] * x + c;
    }
    q
}
fn derivative(a: &[f64]) -> Vec<f64> {
    a.iter()
        .enumerate()
        .skip(1)
        .map(|(i, &v)| i as f64 * v)
        .collect()
}
fn trim(a: &mut Vec<f64>) {
    while a.last() == Some(&0.) {
        a.pop();
    }
}
fn bisect(f: &impl Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> f64 {
    let negative = f(lo) < 0.;
    for _ in 0..180 {
        let mid = lo + (hi - lo) * 0.5;
        if mid == lo || mid == hi {
            break;
        }
        let value = f(mid);
        if value == 0. {
            return mid;
        }
        if (value < 0.) == negative {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo + (hi - lo) * 0.5
}
// If zeros of f' are known, Rolle's theorem partitions f into monotone
// intervals. For g=exp(x*x)*A+B we evaluate exp(-x*x)*g to avoid overflow;
// this positive scaling preserves signs, zeros and the interval ordering.
fn between(f: impl Fn(f64) -> f64, stationary: &[f64], maximum: f64) -> Vec<f64> {
    let mut roots = Vec::new();
    let mut left = 0.;
    let mut a = f(left);
    for right in stationary.iter().copied().chain(std::iter::once(maximum)) {
        let b = f(right);
        if a == 0. {
            roots.push(left);
        }
        if a != 0. && b != 0. && (a < 0.) != (b < 0.) {
            roots.push(bisect(&f, left, right));
        }
        left = right;
        a = b;
    }
    if a == 0. {
        roots.push(maximum);
    }
    roots.dedup_by(|a, b| (*a - *b).abs() <= 32. * f64::EPSILON * a.abs().max(b.abs()));
    roots
}
fn polynomial_roots(a: &[f64], maximum: f64) -> Vec<f64> {
    if a.len() <= 1 {
        return Vec::new();
    }
    let stationary = polynomial_roots(&derivative(a), maximum);
    between(|x| jet(a, x)[0], &stationary, maximum)
}
fn gaussian_roots(mut a: Vec<f64>, mut b: Vec<f64>, maximum: f64) -> Result<Vec<f64>> {
    trim(&mut a);
    trim(&mut b);
    let mut levels = Vec::new();
    // Differentiate exp(x*x)*A+B until B vanishes. The remaining roots are
    // polynomial roots, then recover earlier derivative zeros in reverse.
    while !b.is_empty() {
        levels.push((a.clone(), b.clone()));
        let mut next = derivative(&a);
        next.resize(a.len() + 1, 0.);
        for (i, &c) in a.iter().enumerate() {
            next[i + 1] += 2. * c;
        }
        a = next;
        b = derivative(&b);
        trim(&mut a);
        trim(&mut b);
        let scale = a.iter().chain(&b).map(|v| v.abs()).fold(0., f64::max);
        if !scale.is_finite() {
            return Err(Error::InvalidInput("MBWR branch polynomial overflow"));
        }
        if scale > 0. {
            for c in a.iter_mut().chain(&mut b) {
                *c /= scale;
            }
        }
    }
    // All subsequent bisections stay inside this interval. Check absolute
    // Horner bounds at its upper end before relying on polynomial signs.
    // Finite EOS input alone does not guarantee finite hierarchy arithmetic.
    for (aa, bb) in levels
        .iter()
        .chain(std::iter::once(&(a.clone(), Vec::new())))
    {
        for coefficients in [aa, bb] {
            let bound = coefficients
                .iter()
                .rev()
                .fold(0_f64, |s, c| s * maximum + c.abs());
            if !bound.is_finite() {
                return Err(Error::InvalidInput(
                    "MBWR branch search range overflows polynomial hierarchy",
                ));
            }
        }
    }
    let mut roots = polynomial_roots(&a, maximum);
    for (a, b) in levels.iter().rev() {
        roots = between(
            |x| jet(a, x)[0] + (-x * x).exp() * jet(b, x)[0],
            &roots,
            maximum,
        );
    }
    Ok(roots)
}

impl<'a> PreparedMbwr<'a> {
    // P, P_r, P_rr, P_rrr, P_Tr, P_Trr in SI units.
    fn critical_jet(&self, rho: f64) -> Result<[f64; 6]> {
        positive(rho, "critical molar density")?;
        let r = rho * 1e-3;
        let g = self.model.gamma;
        let exponential = (-g * r * r).exp();
        let evaluate = |k: usize| {
            let mut a = [0.; 10];
            let mut b = [0.; 14];
            a[1] = 1000. * R * if k == 0 { self.temperature } else { 1. };
            for (i, c) in self.pressure_polynomial.iter().enumerate() {
                a[i + 2] = c[k];
            }
            for (i, c) in self.pressure_exponential.iter().enumerate() {
                b[2 * i + 3] = c[k];
            }
            let p = jet(&a, r);
            let q = jet(&b, r);
            [
                p[0] + exponential * q[0],
                (p[1] + exponential * (q[1] - 2. * g * r * q[0])) * 1e-3,
                (p[2]
                    + exponential
                        * (q[2] - 4. * g * r * q[1] + (4. * g * g * r * r - 2. * g) * q[0]))
                    * 1e-6,
                (p[3]
                    + exponential
                        * (q[3] - 6. * g * r * q[2]
                            + (12. * g * g * r * r - 6. * g) * q[1]
                            + (12. * g * g * r - 8. * g.powi(3) * r.powi(3)) * q[0]))
                    * 1e-9,
            ]
        };
        let p = evaluate(0);
        let pt = evaluate(1);
        let q = [p[0], p[1], p[2], p[3], pt[1], pt[2]];
        if q.iter().any(|x| !x.is_finite()) {
            return Err(Error::InvalidInput("non-finite MBWR critical derivatives"));
        }
        Ok(q)
    }
    /// Isolate pressure branches over [0, maximum_density] in mol/m3.
    /// A Gaussian-polynomial derivative hierarchy replaces sampling, so narrow
    /// spinodal pairs are resolved without a fixed density grid.
    pub fn equilibrium(&self, maximum_density: f64) -> Result<MbwrEquilibrium<'a>> {
        positive(maximum_density, "maximum MBWR search density")?;
        let scale = 1000. / self.model.gamma.sqrt();
        let maximum = maximum_density / scale;
        positive(maximum, "scaled maximum MBWR density")?;
        let mut a = vec![0.; 9];
        let mut b = vec![0.; 15];
        a[0] = R * self.temperature * scale;
        let rscale = scale * 1e-3;
        for (i, c) in self.pressure_polynomial.iter().enumerate() {
            a[i + 1] = (i + 2) as f64 * c[0] * rscale.powi((i + 2) as i32);
        }
        for (i, c) in self.pressure_exponential.iter().enumerate() {
            let power = 2 * i + 3;
            let value = c[0] * rscale.powi(power as i32);
            b[power - 1] += power as f64 * value;
            b[power + 1] -= 2. * value;
        }
        if a.iter().chain(&b).any(|x| !x.is_finite()) {
            return Err(Error::InvalidInput("MBWR branch coefficient overflow"));
        }
        let mut spinodals = gaussian_roots(a, b, maximum)?;
        spinodals.retain(|r| *r > 0. && *r < maximum);
        for r in &mut spinodals {
            *r *= scale;
        }
        let mut intervals = Vec::new();
        let mut low = 0.;
        for high in spinodals
            .iter()
            .copied()
            .chain(std::iter::once(maximum_density))
        {
            if self.pressure_density(low + (high - low) * 0.5)?.drho >= 0. {
                intervals.push(Interval {
                    low,
                    high,
                    low_pressure: self.pressure_density(low)?.value,
                    high_pressure: self.pressure_density(high)?.value,
                });
            }
            low = high;
        }
        Ok(MbwrEquilibrium {
            prepared: self.clone(),
            maximum_density,
            spinodals,
            intervals,
        })
    }
}

impl MbwrEquilibrium<'_> {
    pub fn maximum_density(&self) -> f64 {
        self.maximum_density
    }
    /// Positive spinodal densities in ascending order, in mol/m3.
    pub fn spinodal_densities(&self) -> &[f64] {
        &self.spinodals
    }
    fn root(&self, p: f64, interval: &Interval) -> Result<f64> {
        let (mut lo, mut hi) = (interval.low, interval.high);
        let ideal = p / (R * self.prepared.temperature);
        let mut rho = if ideal > lo && ideal < hi {
            ideal
        } else {
            lo + (hi - lo) * 0.5
        };
        for _ in 0..160 {
            let q = self.prepared.pressure_density(rho)?;
            let error = q.value - p;
            let tolerance = 2e-10 * p + 64. * f64::EPSILON * R * self.prepared.temperature * rho;
            let step = error / q.drho;
            if error.abs() <= tolerance
                && (step.abs() <= 32. * f64::EPSILON * rho || hi - lo <= 64. * f64::EPSILON * rho)
            {
                return Ok(rho);
            }
            if error > 0. {
                hi = rho;
            } else {
                lo = rho;
            }
            let next = rho - step;
            let next = if q.drho > 0. && next > lo && next < hi {
                next
            } else {
                lo + (hi - lo) * 0.5
            };
            if next == rho || next == lo || next == hi {
                if error.abs() <= tolerance {
                    return Ok(rho);
                }
                break;
            }
            rho = next;
        }
        Err(Error::NotConverged("MBWR bounded pressure branch"))
    }
    /// Mechanically stable roots, ascending in density. An empty vector means
    /// no root in the requested domain, not a result outside that domain.
    pub fn densities(&self, p: f64) -> Result<Vec<f64>> {
        positive(p, "pressure")?;
        let mut roots = Vec::new();
        for interval in &self.intervals {
            if p >= interval.low_pressure && p <= interval.high_pressure {
                roots.push(self.root(p, interval)?);
            }
        }
        roots.dedup_by(|a, b| (*a - *b).abs() <= 64. * f64::EPSILON * a.abs().max(b.abs()));
        Ok(roots)
    }
    /// Lowest/highest density for vapor/liquid, or minimum residual Gibbs
    /// energy among the mechanically stable roots in the specified domain.
    pub fn volume(&self, p: f64, phase: Phase) -> Result<f64> {
        let roots = self.densities(p)?;
        let rho = match phase {
            Phase::Vapor => *roots.first().ok_or(Error::NoPhysicalRoot)?,
            Phase::Liquid => *roots.last().ok_or(Error::NoPhysicalRoot)?,
            Phase::MinimumGibbs => {
                let mut best = None;
                for rho in roots {
                    let g = self.prepared.residual_tp(1. / rho, 1.)?.gibbs_energy;
                    if best.is_none_or(|(_, bg)| g < bg) {
                        best = Some((rho, g));
                    }
                }
                best.ok_or(Error::NoPhysicalRoot)?.0
            }
        };
        Ok(1. / rho)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn derivative_hierarchy_resolves_narrow_pairs_and_tangencies() {
        let a = vec![1. - 1e-12, -2., 1.];
        for b in [vec![], vec![0., 0., 0.]] {
            let roots = gaussian_roots(a.clone(), b, 2.).unwrap();
            assert_eq!(roots.len(), 2);
            assert!((roots[0] - (1. - 1e-6)).abs() < 1e-9);
            assert!((roots[1] - (1. + 1e-6)).abs() < 1e-9);
        }
        let roots = gaussian_roots(vec![1., -2., 1.], vec![], 2.).unwrap();
        assert_eq!(roots, vec![1.]);
        let roots = gaussian_roots(vec![-0.5], vec![0., 0., 1.], 5.).unwrap();
        assert!(roots.is_empty()); // max x^2 exp(-x^2) = 1/e < 1/2
        let roots = gaussian_roots(vec![-0.3], vec![0., 0., 1.], 5.).unwrap();
        assert_eq!(roots.len(), 2);
        for x in roots {
            assert!((x * x * (-x * x).exp() - 0.3).abs() < 1e-14);
        }
    }
    #[test]
    fn analytic_critical_jet_matches_pressure_and_finite_differences() {
        for record in mbwr_records() {
            let model = Mbwr::from_record(record.source_name).unwrap();
            let prepared = model.prepare(800.).unwrap();
            let r = 5000.;
            let h = 0.1;
            let q = prepared.critical_jet(r).unwrap();
            let p = prepared.pressure_density(r).unwrap();
            for (a, b) in [
                (q[0], p.value),
                (q[1], p.drho),
                (q[2], p.drhorho),
                (q[4], p.dtrho),
            ] {
                assert!((a - b).abs() < 1e-9 * a.abs().max(b.abs()).max(1.));
            }
            let hi = prepared.pressure_density(r + h).unwrap();
            let lo = prepared.pressure_density(r - h).unwrap();
            for (a, b) in [
                (q[3], (hi.drhorho - lo.drhorho) / (2. * h)),
                (q[5], (hi.dtrho - lo.dtrho) / (2. * h)),
            ] {
                assert!(
                    (a - b).abs() < 1e-7 * a.abs().max(b.abs()).max(1e-8),
                    "{} {a} {b}",
                    record.source_name
                );
            }
        }
    }
}
