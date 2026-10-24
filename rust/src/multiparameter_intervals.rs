//! Conservatively bounded pressure-branch discovery, without a fixed density grid.
use super::*;
use std::cell::RefCell;

const MAX_TERMS: usize = 32;
const NATURAL_CACHE_SIZE: usize = 8;

#[derive(Debug, Clone, Copy)]
struct Interval {
    lo: f64,
    hi: f64,
}
impl Interval {
    const ZERO: Self = Self { lo: 0., hi: 0. };
    fn point(x: f64) -> Self {
        Self { lo: x, hi: x }
    }
    fn add(self, b: Self) -> Self {
        if self.lo == 0. && self.hi == 0. {
            return b;
        }
        if b.lo == 0. && b.hi == 0. {
            return self;
        }
        Self {
            lo: (self.lo + b.lo).next_down(),
            hi: (self.hi + b.hi).next_up(),
        }
    }
    fn negative(self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }
    fn multiply(self, b: Self) -> Self {
        if (self.lo == 0. && self.hi == 0.) || (b.lo == 0. && b.hi == 0.) {
            return Self::ZERO;
        }
        let q = [
            self.lo * b.lo,
            self.lo * b.hi,
            self.hi * b.lo,
            self.hi * b.hi,
        ];
        if q.iter().any(|v| v.is_nan()) {
            return Self {
                lo: f64::NEG_INFINITY,
                hi: f64::INFINITY,
            };
        }
        Self {
            lo: q.into_iter().fold(f64::INFINITY, f64::min).next_down(),
            hi: q.into_iter().fold(f64::NEG_INFINITY, f64::max).next_up(),
        }
    }
    fn divide(self, b: Self) -> Self {
        debug_assert!(b.lo > 0. || b.hi < 0.);
        self.multiply(Self {
            lo: (1. / b.hi).next_down(),
            hi: (1. / b.lo).next_up(),
        })
    }
    fn square(self) -> Self {
        let lo = if self.lo <= 0. && self.hi >= 0. {
            0.
        } else {
            let minimum = self.lo.abs().min(self.hi.abs());
            (minimum * minimum).next_down().max(0.)
        };
        Self {
            lo,
            hi: {
                let maximum = self.lo.abs().max(self.hi.abs());
                (maximum * maximum).next_up()
            },
        }
    }
    fn powi(self, n: i32) -> Self {
        let mut value = Self::point(1.);
        for _ in 0..n {
            value = value.multiply(self);
        }
        value
    }
    fn exp(self) -> Self {
        let (mut lo, mut hi) = (self.lo.exp(), self.hi.exp());
        // Scalar exp is platform supplied. Four outward ulps are an explicit
        // transcendental rounding guard, not a correctly-rounded-libm proof.
        for _ in 0..4 {
            lo = lo.next_down();
            hi = hi.next_up();
        }
        Self { lo: lo.max(0.), hi }
    }
    fn sign(self) -> i32 {
        if self.lo > 0. {
            1
        } else if self.hi < 0. {
            -1
        } else {
            0
        }
    }
}

#[derive(Debug, Clone)]
struct BoundTerm {
    polynomial: [[Interval; COEFFICIENTS]; 4],
    lengths: [usize; 4],
    exponent: DensityExponent,
    coefficient: f64,
}
impl BoundTerm {
    fn new(term: &DensityTerm) -> Self {
        let mut q = [Interval::ZERO; COEFFICIENTS];
        let degree = match term.exponent {
            DensityExponent::Plain => 0,
            DensityExponent::Power(l) => {
                q[l as usize - 1] = Interval::point(-f64::from(l));
                l as usize - 1
            }
            DensityExponent::Gaussian { eta, epsilon, .. } => {
                q[0] = Interval::point(2.)
                    .multiply(Interval::point(eta))
                    .multiply(Interval::point(epsilon));
                q[1] = Interval::point(-2.).multiply(Interval::point(eta));
                1
            }
        };
        let d = term.density_power;
        let mut polynomial = [[Interval::ZERO; COEFFICIENTS]; 4];
        polynomial[0][d + 1] = Interval::point(d as f64);
        for (i, &v) in q[..=degree].iter().enumerate() {
            polynomial[0][d + 2 + i] = polynomial[0][d + 2 + i].add(v);
        }
        for k in 1..4 {
            for i in 0..term.lengths[k - 1] {
                if i > 0 {
                    polynomial[k][i - 1] = polynomial[k][i - 1]
                        .add(Interval::point(i as f64).multiply(polynomial[k - 1][i]));
                }
                for (j, &v) in q[..=degree].iter().enumerate() {
                    polynomial[k][i + j] =
                        polynomial[k][i + j].add(v.multiply(polynomial[k - 1][i]));
                }
            }
        }
        Self {
            polynomial,
            lengths: term.lengths,
            exponent: term.exponent,
            coefficient: term.coefficient,
        }
    }
    fn exponential(&self, x: Interval) -> Interval {
        let exponent = match self.exponent {
            DensityExponent::Plain => Interval::ZERO,
            DensityExponent::Power(l) => x.powi(l).negative(),
            DensityExponent::Gaussian {
                eta,
                epsilon,
                temperature,
            } => x
                .add(Interval::point(-epsilon))
                .square()
                .multiply(Interval::point(-eta))
                .add(Interval::point(-temperature)),
        };
        exponent.exp()
    }
    fn bound_with_exponential(&self, x: Interval, k: usize, exponential: Interval) -> Interval {
        let poly = self.polynomial[k][..self.lengths[k]]
            .iter()
            .rev()
            .fold(Interval::ZERO, |v, &a| v.multiply(x).add(a));
        poly.multiply(exponential)
            .multiply(Interval::point(self.coefficient))
    }
    #[cfg(test)]
    fn bound(&self, x: Interval, k: usize) -> Interval {
        self.bound_with_exponential(x, k, self.exponential(x))
    }
}

/// A density box containing a spinodal or an unresolved nearly stationary
/// region. Pressures within its conservative pressure box are rejected by
/// complete-root queries, rather than silently dropping marginal roots.
#[derive(Debug, Clone, Copy)]
pub struct MultiparameterSpinodalBracket {
    pub minimum_density: f64,
    pub maximum_density: f64,
    pub minimum_pressure: f64,
    pub maximum_pressure: f64,
    /// The derivative bounds and endpoint signs isolate one simple spinodal.
    /// False denotes a region whose stationary-point count is unresolved.
    pub isolated: bool,
}
#[derive(Debug, Clone, Copy)]
struct StableInterval {
    low: f64,
    high: f64,
    low_pressure: f64,
    high_pressure: f64,
}

/// Reusable mechanically stable pressure branches in an explicit density
/// domain. Discovery allocates; scalar selection and caller-buffer queries
/// allocate no memory. Root completeness is conditional on the documented
/// outward arithmetic/exp guard; unresolved pressure windows remain errors.
#[derive(Debug, Clone)]
pub struct MultiparameterEquilibrium<'a> {
    prepared: PreparedMultiparameter<'a>,
    maximum_density: f64,
    brackets: Vec<MultiparameterSpinodalBracket>,
    intervals: Vec<StableInterval>,
}

#[derive(Clone, Copy)]
struct NaturalEntry {
    key: Option<[u64; 2]>,
    exponentials: [Interval; MAX_TERMS],
    values: [Option<Interval>; 4],
}
impl NaturalEntry {
    const EMPTY: Self = Self {
        key: None,
        exponentials: [Interval::ZERO; MAX_TERMS],
        values: [None; 4],
    };
}
struct NaturalCache {
    entries: [NaturalEntry; NATURAL_CACHE_SIZE],
    next: usize,
}
struct Bounds {
    terms: Vec<BoundTerm>,
    cache: RefCell<NaturalCache>,
}
impl Bounds {
    fn new(terms: Vec<BoundTerm>) -> Result<Self> {
        if terms.len() > MAX_TERMS {
            return Err(Error::InvalidInput(
                "multiparameter bound cache term capacity",
            ));
        }
        Ok(Self {
            terms,
            cache: RefCell::new(NaturalCache {
                entries: [NaturalEntry::EMPTY; NATURAL_CACHE_SIZE],
                next: 0,
            }),
        })
    }
    fn natural(&self, x: Interval, k: usize) -> Interval {
        let key = [x.lo.to_bits(), x.hi.to_bits()];
        let mut cache = self.cache.borrow_mut();
        let slot = if let Some(i) = cache
            .entries
            .iter()
            .position(|entry| entry.key == Some(key))
        {
            i
        } else {
            let i = cache.next;
            cache.next = (i + 1) % NATURAL_CACHE_SIZE;
            let entry = &mut cache.entries[i];
            entry.key = Some(key);
            entry.values = [None; 4];
            for (out, term) in entry.exponentials.iter_mut().zip(&self.terms) {
                *out = term.exponential(x);
            }
            i
        };
        let entry = &mut cache.entries[slot];
        if let Some(value) = entry.values[k] {
            return value;
        }
        let mut sum = match k {
            0 => x,
            1 => Interval::point(1.),
            _ => Interval::ZERO,
        };
        for (term, &exponential) in self.terms.iter().zip(&entry.exponentials) {
            sum = sum.add(term.bound_with_exponential(x, k, exponential));
        }
        entry.values[k] = Some(sum);
        sum
    }
    fn evaluate(&self, x: Interval, k: usize) -> Interval {
        let natural = self.natural(x, k);
        if x.lo == x.hi || k >= 3 {
            return natural;
        }
        let mid = x.lo + (x.hi - x.lo) * 0.5;
        let center = Interval::point(mid);
        let dx = x.add(center.negative());
        let slope = if k == 2 {
            self.natural(x, k + 1)
        } else {
            self.natural(center, k + 1)
        };
        let mut taylor = self.natural(center, k).add(slope.multiply(dx));
        if k < 2 {
            taylor = taylor.add(
                self.natural(x, k + 2)
                    .multiply(dx.square())
                    .multiply(Interval::point(0.5)),
            );
        }
        // The first-order E bound uses F over the whole interval. Pressure
        // and stiffness use second-order centered forms. Intersecting these
        // with Horner bounds removes density-dependency inflation near flat
        // critical regions without excluding any enclosed function value.
        let result = Interval {
            lo: natural.lo.max(taylor.lo),
            hi: natural.hi.min(taylor.hi),
        };
        debug_assert!(
            result.lo <= result.hi,
            "inconsistent density derivative bounds"
        );
        result
    }
    fn isolate(&self, mut lo: f64, mut hi: f64) -> [f64; 2] {
        for _ in 0..64 {
            let mid = lo + (hi - lo) * 0.5;
            let derivative = self.evaluate(Interval { lo, hi }, 2);
            let step = self.evaluate(Interval::point(mid), 1).divide(derivative);
            let next = Interval::point(mid).add(step.negative());
            let nl = lo.max(next.lo);
            let nh = hi.min(next.hi);
            if nl > nh || nh - nl >= 0.95 * (hi - lo) {
                break;
            }
            lo = nl;
            hi = nh;
            if hi - lo < 64. * f64::EPSILON * hi.abs().max(1.) {
                break;
            }
        }
        [lo, hi]
    }
    fn spinodals(&self, maximum: f64) -> Result<Vec<([f64; 2], bool)>> {
        let mut pending = vec![(0., maximum, 0)];
        let mut found = Vec::new();
        let mut visited = 0;
        while let Some((lo, hi, depth)) = pending.pop() {
            visited += 1;
            if visited > 100000 {
                return Err(Error::NotConverged(
                    "multiparameter interval discovery budget",
                ));
            }
            let x = Interval { lo, hi };
            let d = self.evaluate(x, 1);
            if d.sign() != 0 {
                continue;
            }
            let e = self.evaluate(x, 2);
            if e.sign() != 0 {
                let a = self.evaluate(Interval::point(lo), 1).sign();
                let b = self.evaluate(Interval::point(hi), 1).sign();
                if a != 0 && b != 0 {
                    if a != b {
                        found.push((self.isolate(lo, hi), true));
                    }
                    continue;
                }
            }
            let width = hi - lo;
            let mid = lo + width * 0.5;
            if depth >= 60
                || mid == lo
                || mid == hi
                || width <= 1e-12 * hi.abs().max(1.)
                || (width <= 1e-5 * hi.abs().max(1.)
                    && self.evaluate(Interval::point(mid), 1).sign() == 0)
            {
                found.push(([lo, hi], false));
                continue;
            }
            pending.push((mid, hi, depth + 1));
            pending.push((lo, mid, depth + 1));
        }
        found.sort_by(|a, b| a.0[0].total_cmp(&b.0[0]));
        let mut merged: Vec<([f64; 2], bool)> = Vec::new();
        for (x, isolated) in found {
            if let Some(previous) = merged.last_mut()
                && x[0] <= previous.0[1]
            {
                previous.0[1] = previous.0[1].max(x[1]);
                previous.1 = false;
            } else {
                merged.push((x, isolated));
            }
        }
        Ok(merged)
    }
}

impl<'a> PreparedMultiparameter<'a> {
    /// Discover stable branches over [0, maximum_density]. The boundary is a
    /// search limit, rather than an EOS validity or phase-stability guarantee.
    pub fn equilibrium(&self, maximum_density: f64) -> Result<MultiparameterEquilibrium<'a>> {
        positive(maximum_density, "maximum multiparameter search density")?;
        // Reject an undefined outer EOS endpoint before doing interval work.
        self.pressure_density(maximum_density)?;
        let maximum = maximum_density / self.eos.rc();
        positive(maximum, "maximum reduced density")?;
        let mut terms = Vec::new();
        self.visit_density_terms(|t| {
            terms.push(BoundTerm::new(&t));
            Ok(())
        })?;
        let bounds = Bounds::new(terms)?;
        let boxes = bounds.spinodals(maximum)?;
        let scale = self.eos.gas_constant() * self.temperature * self.eos.rc();
        let mut brackets = Vec::new();
        let mut intervals = Vec::new();
        let mut low = 0.;
        for (x, isolated) in boxes {
            if x[0] > low
                && self
                    .pressure_density((low + (x[0] - low) * 0.5) * self.eos.rc())?
                    .1
                    > 0.
            {
                intervals.push(StableInterval {
                    low: low * self.eos.rc(),
                    high: x[0] * self.eos.rc(),
                    low_pressure: self.pressure_density(low * self.eos.rc())?.0,
                    high_pressure: self.pressure_density(x[0] * self.eos.rc())?.0,
                });
            }
            let p = bounds
                .evaluate(Interval { lo: x[0], hi: x[1] }, 0)
                .multiply(Interval::point(scale));
            if !p.lo.is_finite() || !p.hi.is_finite() {
                return Err(Error::InvalidInput(
                    "nonfinite multiparameter spinodal pressure bounds",
                ));
            }
            brackets.push(MultiparameterSpinodalBracket {
                minimum_density: x[0] * self.eos.rc(),
                maximum_density: x[1] * self.eos.rc(),
                minimum_pressure: p.lo,
                maximum_pressure: p.hi,
                isolated,
            });
            low = x[1];
        }
        if low < maximum
            && self
                .pressure_density((low + (maximum - low) * 0.5) * self.eos.rc())?
                .1
                > 0.
        {
            intervals.push(StableInterval {
                low: low * self.eos.rc(),
                high: maximum_density,
                low_pressure: self.pressure_density(low * self.eos.rc())?.0,
                high_pressure: self.pressure_density(maximum_density)?.0,
            });
        }
        Ok(MultiparameterEquilibrium {
            prepared: self.clone(),
            maximum_density,
            brackets,
            intervals,
        })
    }
}
impl Multiparameter {
    /// Prepare reusable stable branches in an explicit molar-density domain.
    ///
    /// ```
    /// use thermopack_rs::{Multiparameter, Phase};
    /// let model = Multiparameter::new("C3")?;
    /// let branches = model.prepare_equilibrium(280.0, 8.0 * model.rc())?;
    /// let mut buffer = vec![0.0; branches.maximum_root_count()];
    /// let n = branches.densities_into(1e6, &mut buffer)?;
    /// assert!(n > 0);
    /// let volume = branches.volume(1e6, Phase::MinimumGibbs)?;
    /// assert!(volume > 0.0);
    /// # Ok::<(), thermopack_rs::Error>(())
    /// ```
    pub fn prepare_equilibrium(
        &self,
        temperature: f64,
        maximum_density: f64,
    ) -> Result<MultiparameterEquilibrium<'_>> {
        self.prepare(temperature)?.equilibrium(maximum_density)
    }
}
impl MultiparameterEquilibrium<'_> {
    pub fn maximum_density(&self) -> f64 {
        self.maximum_density
    }
    pub fn spinodal_brackets(&self) -> &[MultiparameterSpinodalBracket] {
        &self.brackets
    }
    pub fn maximum_root_count(&self) -> usize {
        self.intervals.len()
    }
    fn check_pressure(&self, p: f64) -> Result<()> {
        positive(p, "multiparameter pressure")?;
        if self
            .brackets
            .iter()
            .any(|b| p >= b.minimum_pressure && p <= b.maximum_pressure)
        {
            return Err(Error::NotConverged(
                "multiparameter pressure intersects unresolved spinodal box",
            ));
        }
        Ok(())
    }
    fn root(&self, p: f64, interval: StableInterval) -> Result<f64> {
        let (mut lo, mut hi) = (interval.low, interval.high);
        let rt = self.prepared.eos.gas_constant() * self.prepared.temperature;
        let ideal = p / rt;
        let mut rho = if ideal > lo && ideal < hi {
            ideal
        } else {
            lo + (hi - lo) * 0.5
        };
        for _ in 0..160 {
            let (value, derivative) = self.prepared.pressure_density(rho)?;
            let error = value - p;
            let tolerance = 4e-12 * p + 64. * f64::EPSILON * (rt * rho + p);
            let step = error / derivative;
            if derivative > 0.
                && error.abs() <= tolerance
                && (step.abs() <= 64. * f64::EPSILON * rho || hi - lo <= 128. * f64::EPSILON * rho)
            {
                return Ok(rho);
            }
            if error > 0. {
                hi = rho;
            } else {
                lo = rho;
            }
            let next = rho - step;
            let next = if derivative > 0. && next > lo && next < hi {
                next
            } else {
                lo + (hi - lo) * 0.5
            };
            if next == rho || next == lo || next == hi {
                let cancellation_tolerance = 4e-12 * p + self.prepared.pressure_roundoff(rho)?;
                let cap = 2e-8 * p + 256. * f64::EPSILON * rt * rho;
                let allowed = tolerance.max(cancellation_tolerance.min(cap));
                let mut best = (rho, error.abs());
                for trial in [lo, hi, rho.next_down(), rho.next_up()] {
                    if trial >= interval.low && trial <= interval.high && trial > 0. {
                        let (value, stiffness) = self.prepared.pressure_density(trial)?;
                        if stiffness > 0. && (value - p).abs() < best.1 {
                            best = (trial, (value - p).abs());
                        }
                    }
                }
                if derivative > 0. && best.1 <= allowed {
                    return Ok(best.0);
                }
                return Err(Error::NotConverged(
                    "multiparameter pressure unresolved at density precision",
                ));
            }
            rho = next;
        }
        Err(Error::NotConverged(
            "multiparameter bounded pressure branch",
        ))
    }
    /// Write all resolved stable densities in ascending order. The buffer must
    /// hold maximum_root_count() entries. Returns zero when no bounded root exists.
    pub fn densities_into(&self, p: f64, output: &mut [f64]) -> Result<usize> {
        self.check_pressure(p)?;
        if output.len() < self.maximum_root_count() {
            return Err(Error::InvalidInput("multiparameter root buffer too small"));
        }
        let mut count = 0;
        for &i in &self.intervals {
            if p >= i.low_pressure && p <= i.high_pressure {
                output[count] = self.root(p, i)?;
                count += 1;
            }
        }
        Ok(count)
    }
    pub fn densities(&self, p: f64) -> Result<Vec<f64>> {
        let mut roots = vec![0.; self.maximum_root_count()];
        let n = self.densities_into(p, &mut roots)?;
        roots.truncate(n);
        Ok(roots)
    }
    /// Lowest/highest density for vapor/liquid, or minimum Gibbs among all
    /// resolved stable roots within the domain. A scalar query allocates nothing.
    pub fn volume(&self, p: f64, phase: Phase) -> Result<f64> {
        self.check_pressure(p)?;
        if phase == Phase::Liquid {
            for &i in self.intervals.iter().rev() {
                if p >= i.low_pressure && p <= i.high_pressure {
                    return Ok(1. / self.root(p, i)?);
                }
            }
            return Err(Error::NoPhysicalRoot);
        }
        let mut best = None;
        for &i in &self.intervals {
            if p < i.low_pressure || p > i.high_pressure {
                continue;
            }
            if phase == Phase::Vapor {
                return Ok(1. / self.root(p, i)?);
            }
            let rho = self.root(p, i)?;
            let a = self
                .prepared
                .residual_at_delta(rho / self.prepared.eos.rc())?;
            let g = a.alpha
                + (rho / self.prepared.eos.rc()).ln()
                + p / (rho * self.prepared.eos.gas_constant() * self.prepared.temperature);
            if best.is_none_or(|(_, value)| g < value) {
                best = Some((rho, g));
            }
        }
        Ok(1. / best.ok_or(Error::NoPhysicalRoot)?.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn uncached_natural(bounds: &Bounds, x: Interval, k: usize) -> Interval {
        let mut value = match k {
            0 => x,
            1 => Interval::point(1.),
            _ => Interval::ZERO,
        };
        for term in &bounds.terms {
            value = value.add(term.bound(x, k));
        }
        value
    }
    fn uncached_evaluate(bounds: &Bounds, x: Interval, k: usize) -> Interval {
        let value = uncached_natural(bounds, x, k);
        if x.lo == x.hi || k >= 3 {
            return value;
        }
        let mid = x.lo + (x.hi - x.lo) * 0.5;
        let center = Interval::point(mid);
        let dx = x.add(center.negative());
        let mut taylor = uncached_natural(bounds, center, k)
            .add(uncached_natural(bounds, center, k + 1).multiply(dx));
        if k < 2 {
            taylor = taylor.add(
                uncached_natural(bounds, x, k + 2)
                    .multiply(dx.square())
                    .multiply(Interval::point(0.5)),
            );
        }
        if k == 2 {
            taylor = uncached_natural(bounds, center, k)
                .add(uncached_natural(bounds, x, k + 1).multiply(dx));
        }
        Interval {
            lo: value.lo.max(taylor.lo),
            hi: value.hi.min(taylor.hi),
        }
    }
    #[test]
    fn interval_cache_preserves_bounds_bitwise_across_derivatives_and_eviction() {
        for model in [
            Multiparameter::new("C3").unwrap(),
            Multiparameter::new("R134A").unwrap(),
            Multiparameter::new("O-H2").unwrap(),
            Multiparameter::lennard_jones("AR", LennardJonesKind::Full).unwrap(),
            Multiparameter::reduced_lennard_jones_truncated_shifted(),
        ] {
            for tr in [0.6, 0.999999, 1., 1.5] {
                let prepared = model.prepare(tr * model.tc()).unwrap();
                let mut terms = Vec::new();
                prepared
                    .visit_density_terms(|t| {
                        terms.push(BoundTerm::new(&t));
                        Ok(())
                    })
                    .unwrap();
                let bounds = Bounds::new(terms).unwrap();
                let mut intervals = vec![
                    Interval::point(0.),
                    Interval { lo: -0., hi: 0. },
                    Interval { lo: 0., hi: 8. },
                ];
                for i in 1..=20 {
                    let d = f64::from(i) / 5.;
                    intervals.push(Interval::point(d));
                    intervals.push(Interval {
                        lo: d - 1e-5,
                        hi: d + 1e-5,
                    });
                }
                for x in intervals
                    .iter()
                    .copied()
                    .chain(intervals.iter().rev().copied())
                {
                    for k in [2, 1, 3, 0, 1, 2] {
                        let actual = bounds.evaluate(x, k);
                        let expected = uncached_evaluate(&bounds, x, k);
                        assert_eq!(
                            [actual.lo.to_bits(), actual.hi.to_bits()],
                            [expected.lo.to_bits(), expected.hi.to_bits()],
                            "{} {tr} {x:?} derivative {k}",
                            model.record.ident
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn centered_bounds_resolve_narrow_pairs_and_report_unresolved_contact() {
        for epsilon in [1e-8, 1e-12, 0.] {
            let a = DensityTerm::new(1, DensityExponent::Plain, -1., 0., 0).unwrap();
            let c = (1. - epsilon) / 6.;
            let b = DensityTerm::new(2, DensityExponent::Plain, c, 0., 0).unwrap();
            let bounds = Bounds::new(vec![BoundTerm::new(&a), BoundTerm::new(&b)]).unwrap();
            let roots = bounds.spinodals(2.).unwrap();
            if epsilon > 0. {
                assert_eq!(roots.len(), 2);
                let s = (1. - 6. * c).sqrt();
                let expected = [1. / (1. + s), 1. / (1. - s)];
                for ((interval, isolated), root) in roots.iter().zip(expected) {
                    assert!(*isolated && interval[1] - interval[0] < 1e-7);
                    assert!(interval[0] - 1e-9 <= root && root <= interval[1] + 1e-9);
                }
            } else {
                assert!(
                    roots
                        .iter()
                        .any(|(x, isolated)| !*isolated && x[0] < 1. && x[1] > 1.)
                );
            }
        }
    }
    #[test]
    fn interval_bounds_enclose_scalar_derivatives_and_hidden_narrow_pairs() {
        for m in [
            Multiparameter::new("C3").unwrap(),
            Multiparameter::new("R134A").unwrap(),
            Multiparameter::new("O-H2").unwrap(),
            Multiparameter::lennard_jones("AR", LennardJonesKind::Full).unwrap(),
            Multiparameter::reduced_lennard_jones_truncated_shifted(),
        ] {
            for tr in [0.6, 0.95, 0.999999, 1., 1.5] {
                let p = m.prepare(tr * m.tc()).unwrap();
                let mut terms = Vec::new();
                p.visit_density_terms(|t| {
                    terms.push(BoundTerm::new(&t));
                    Ok(())
                })
                .unwrap();
                let b = Bounds::new(terms).unwrap();
                for (lo, hi) in [(0., 0.1), (0.9, 1.1), (1.999, 2.001), (3., 4.)] {
                    let intervals = [
                        b.evaluate(Interval { lo, hi }, 0),
                        b.evaluate(Interval { lo, hi }, 1),
                        b.evaluate(Interval { lo, hi }, 2),
                    ];
                    for j in 0..=10 {
                        let d = lo + (hi - lo) * f64::from(j) / 10.;
                        let q = p.pressure_density_jet(d * m.rc()).unwrap();
                        let values = [
                            q.value / (m.gas_constant() * p.temperature * m.rc()),
                            q.drho / (m.gas_constant() * p.temperature),
                            q.drhorho * m.rc() / (m.gas_constant() * p.temperature),
                        ];
                        for k in 0..3 {
                            assert!(
                                values[k] >= intervals[k].lo - 1e-12
                                    && values[k] <= intervals[k].hi + 1e-12
                            );
                        }
                    }
                }
                let roots = p.equilibrium(6. * m.rc()).unwrap_or_else(|e| {
                    panic!(
                        "{} {tr}: {e}; D(1)={:?}",
                        m.record.ident,
                        b.evaluate(Interval::point(1.), 1)
                    )
                });
                if tr < 0.99999 {
                    assert!(roots.brackets.len() >= 2, "{} {tr}", m.record.ident);
                }
            }
        }
    }
}
