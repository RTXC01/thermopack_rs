//! Spinodal grids and property intersections from `src/spinodal.f90`.
use crate::mixture_spinodal::{
    SpinodalCurve, SpinodalTraceEnd, spinodal_temperature_at_volume, spinodal_volume_at_temperature,
};
use crate::{Cubic, Error, Phase, Result, error::positive, spinodal::SpinodalPoint};

/// Scalar molar property specified on a local spinodal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpinodalProperty {
    Entropy,
    Volume,
    Enthalpy,
    Temperature,
    Pressure,
    InternalEnergy,
}

/// The two candidates selected and ordered by the fork's pure extrema mapper.
/// These are not guaranteed global extrema. In particular, its minimum
/// pressure is evaluated using the maximum candidate's temperature.
#[derive(Debug, Clone, Copy)]
pub struct BaselineSpinodalExtrema {
    pub minimum: SpinodalPoint,
    pub maximum: SpinodalPoint,
    /// Pressure candidates at the common initial pressure cannot be reliably
    /// ordered by their residuals. Coordinates are candidates in either order.
    pub ordering_ambiguous: bool,
}

impl SpinodalProperty {
    /// The source `spin_locate_from_*` identifier.
    pub fn code(self) -> usize {
        match self {
            Self::Entropy => 1,
            Self::Volume => 2,
            Self::Enthalpy => 3,
            Self::Temperature => 4,
            Self::Pressure => 5,
            Self::InternalEnergy => 6,
        }
    }
    pub fn evaluate(self, model: &Cubic, z: &[f64], point: SpinodalPoint) -> Result<f64> {
        Ok(match self {
            Self::Temperature => point.temperature,
            Self::Volume => point.volume,
            Self::Pressure => {
                model
                    .prepare(point.temperature, z)?
                    .pressure(point.volume)?
                    .value
            }
            p => {
                let c = model
                    .prepare(point.temperature, z)?
                    .caloric_tv(point.volume)?;
                match p {
                    Self::Entropy => c.entropy,
                    Self::Enthalpy => c.internal_energy + c.pressure * point.volume,
                    _ => c.internal_energy,
                }
            }
        })
    }
}

/// The fork's pure-fluid grid: 100 liquid points followed by 100 vapor points
/// in reverse order. The two critical entries retain the component's stored
/// T/P and phase volume-root convention rather than imposing a zero eigenvalue.
/// Optional start temperatures change the branch ends before grid construction.
pub fn baseline_pure_spinodal_curve(
    model: &Cubic,
    z: &[f64],
    initial_pressure: f64,
    liquid_temperature: Option<f64>,
    vapor_temperature: Option<f64>,
) -> Result<SpinodalCurve> {
    positive(initial_pressure, "spinodal initial pressure")?;
    model.prepare(300.0, z)?;
    if z.iter().filter(|&&x| x > 0.0).count() != 1 {
        return Err(Error::InvalidInput(
            "a pure spinodal composition is required",
        ));
    }
    // singleCompStabilityLimit unconditionally reads critical component 1.
    let tc = model.components()[0].critical_temperature;
    let pc = model.components()[0].critical_pressure;
    let critical = model.pseudo_critical(z)?;
    let mut branches = Vec::with_capacity(2);
    for (phase, temperature) in [
        (Phase::Liquid, liquid_temperature),
        (Phase::Vapor, vapor_temperature),
    ] {
        let mut start = crate::spinodal::spinodal_point(model, initial_pressure, phase, z)?;
        if let Some(t) = temperature {
            start = spinodal_volume_at_temperature(model, t, z, start.volume)?;
        }
        let mut points = vec![start];
        let dt = (tc - start.temperature) / 99.0;
        let mut t = start.temperature;
        for _ in 1..99 {
            t += dt;
            let v = crate::spinodal::limit_volume(
                model,
                t,
                critical.volume,
                phase == Phase::Liquid,
                z,
            )?;
            points.push(SpinodalPoint {
                temperature: t,
                volume: v,
                pressure: model.prepare(t, z)?.pressure(v)?.value,
            });
        }
        let volume = model.prepare(tc, z)?.baseline_volume(pc, phase)?;
        points.push(SpinodalPoint {
            temperature: tc,
            volume,
            pressure: pc,
        });
        branches.push(points);
    }
    let mut points = branches.remove(0);
    points.extend(branches.remove(0).into_iter().rev());
    Ok(SpinodalCurve {
        points,
        end: if vapor_temperature.is_some() {
            SpinodalTraceEnd::MinimumTemperature
        } else {
            SpinodalTraceEnd::InitialPressure
        },
    })
}

/// All bracketed intersections along an already traced spinodal. Reusing a
/// curve avoids remapping it for every property target. Each returned state
/// closes the property and minimum Helmholtz eigenvalue equations. A discrete
/// trace cannot establish the absence of unbracketed tangent crossings.
pub fn spinodal_intersections(
    model: &Cubic,
    z: &[f64],
    curve: &SpinodalCurve,
    property: SpinodalProperty,
    target: f64,
) -> Result<Vec<SpinodalPoint>> {
    if !target.is_finite() {
        return Err(Error::InvalidInput("finite spinodal target required"));
    }
    model.prepare(300.0, z)?;
    let mut roots: Vec<SpinodalPoint> = Vec::new();
    for pair in curve.points.windows(2) {
        let f0 = property.evaluate(model, z, pair[0])? - target;
        let f1 = property.evaluate(model, z, pair[1])? - target;
        if f0 * f1 > 0.0 {
            continue;
        }
        if (pair[0].volume - pair[1].volume).abs() <= 1e-14 * pair[0].volume.abs() {
            continue;
        }
        let point = refine_property(
            model,
            z,
            pair[0],
            pair[1],
            property,
            target,
            Refinement::Physical,
        )?;
        if roots
            .iter()
            .all(|p| (p.volume / point.volume - 1.0).abs() > 1e-7)
        {
            roots.push(point);
        }
    }
    Ok(roots)
}

fn property_gradient(
    model: &Cubic,
    z: &[f64],
    p: SpinodalPoint,
    property: SpinodalProperty,
) -> Result<f64> {
    if property == SpinodalProperty::Pressure {
        return Ok(1.0);
    }
    let h = 1e-5;
    let state = model.prepare(p.temperature, z)?;
    let pressure = state.pressure(p.volume)?;
    let pvt = (model
        .prepare(p.temperature * (1.0 + h), z)?
        .pressure(p.volume)?
        .dv
        - model
            .prepare(p.temperature * (1.0 - h), z)?
            .pressure(p.volume)?
            .dv)
        / (2.0 * h * p.temperature);
    let pvv = (state.pressure(p.volume * (1.0 + h))?.dv - state.pressure(p.volume * (1.0 - h))?.dv)
        / (2.0 * h * p.volume);
    let vt = -pvt / pvv;
    let pt = pressure.dt + pressure.dv * vt;
    let value = |t, v| {
        property.evaluate(
            model,
            z,
            SpinodalPoint {
                temperature: t,
                volume: v,
                pressure: 0.0,
            },
        )
    };
    let dt = (value(p.temperature * (1.0 + h), p.volume)?
        - value(p.temperature * (1.0 - h), p.volume)?)
        / (2.0 * h * p.temperature);
    let dv = (value(p.temperature, p.volume * (1.0 + h))?
        - value(p.temperature, p.volume * (1.0 - h))?)
        / (2.0 * h * p.volume);
    let scale = match property {
        SpinodalProperty::Entropy => 1e4,
        SpinodalProperty::Enthalpy | SpinodalProperty::InternalEnergy => {
            1e5 / model.components()[0].critical_temperature
        }
        _ => 1.0,
    };
    Ok((dt + dv * vt) / pt * scale)
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Refinement {
    Physical,
    SourceValue,
    SourceGradient,
}

fn refine_property(
    model: &Cubic,
    z: &[f64],
    left: SpinodalPoint,
    right: SpinodalPoint,
    property: SpinodalProperty,
    target: f64,
    mode: Refinement,
) -> Result<SpinodalPoint> {
    let baseline = mode != Refinement::Physical;
    let gradient = mode == Refinement::SourceGradient;
    let h = 1e-5;
    let pvt = (model
        .prepare(left.temperature * (1.0 + h), z)?
        .pressure(left.volume)?
        .dv
        - model
            .prepare(left.temperature * (1.0 - h), z)?
            .pressure(left.volume)?
            .dv)
        / (2.0 * h * left.temperature);
    let state = model.prepare(left.temperature, z)?;
    let pvv = (state.pressure(left.volume * (1.0 + h))?.dv
        - state.pressure(left.volume * (1.0 - h))?.dv)
        / (2.0 * h * left.volume);
    let vt = -pvt / pvv;
    let eval = |v: f64| -> Result<(f64, SpinodalPoint)> {
        let guess = if baseline {
            left.temperature + (v - left.volume) / vt
        } else {
            left.temperature
                + (right.temperature - left.temperature) * (v - left.volume)
                    / (right.volume - left.volume)
        };
        let actual = if baseline {
            let state = model.prepare(guess, z)?;
            if state.mixture_stability(v)?.eigenvalue.abs() < 2e-8 {
                SpinodalPoint {
                    temperature: guess,
                    volume: v,
                    pressure: state.pressure(v)?.value,
                }
            } else {
                spinodal_volume_at_temperature(model, guess, z, v)?
            }
        } else {
            spinodal_temperature_at_volume(model, v, z, guess)?
        };
        let f = if gradient {
            property_gradient(model, z, actual, property)?
        } else {
            (property.evaluate(model, z, actual)? - target)
                / if !baseline && property == SpinodalProperty::Volume {
                    target.abs().max(1e-12)
                } else {
                    target.abs().max(1.0)
                }
        };
        let reported = if baseline {
            SpinodalPoint {
                temperature: guess,
                volume: v,
                pressure: model.prepare(guess, z)?.pressure(v)?.value,
            }
        } else {
            actual
        };
        Ok((f, reported))
    };
    let (mut a, mut b) = (left.volume.min(right.volume), left.volume.max(right.volume));
    let (mut fa, pa) = eval(a)?;
    let (fb, pb) = eval(b)?;
    if !baseline {
        if fa.abs() < 1e-10 {
            return Ok(pa);
        }
        if fb.abs() < 1e-10 {
            return Ok(pb);
        }
    }
    if fa * fb > 0.0 {
        return Err(Error::NotConverged("spinodal property bracket"));
    }
    if baseline {
        // Source King/Pegasus iteration and absolute objective tolerance.
        // In particular the source volume objective is normalized by 1, not V.
        let (mut x0, mut x1, mut f0, mut f1) = (a, b, fa, fb);
        let (mut x2, mut f2) = (0.0, 0.0);
        let mut pass = false;
        for _ in 0..100 {
            if !pass {
                x2 = x1 - (x0 - x1) * f1 / (f0 - f1);
                let (f, p) = eval(x2)?;
                f2 = f;
                if f2.abs() < 1e-8 {
                    return Ok(p);
                }
                if f2 * f1 < 0.0 {
                    std::mem::swap(&mut x0, &mut x1);
                    std::mem::swap(&mut f0, &mut f1);
                }
            }
            pass = false;
            if f2 * f1 > 0.0 {
                f0 = f0 * f1 / (f1 + f2);
                x1 = x2;
                f1 = f2;
                x2 = x1 - (x0 - x1) * f1 / (f0 - f1);
                let (f, p) = eval(x2)?;
                f2 = f;
                if f2.abs() < 1e-8 {
                    return Ok(p);
                }
            }
            if f2 * f1 < 0.0 {
                x0 = x1;
                f0 = f1;
                x1 = x2;
                f1 = f2;
            } else {
                pass = true;
            }
        }
        return Err(Error::NotConverged("source Pegasus spinodal intersection"));
    }
    for _ in 0..90 {
        let v = 0.5 * (a + b);
        let (f, point) = eval(v)?;
        if f.abs() < if baseline { 1e-8 } else { 1e-10 } {
            return Ok(point);
        }
        if fa * f <= 0.0 {
            b = v;
        } else {
            a = v;
            fa = f;
        }
    }
    Err(Error::NotConverged("spinodal property intersection"))
}

/// Source pure branch selection, including vapor precedence when both branches
/// are requested. The source evaluates a corrected volume during refinement
/// but returns the uncorrected trial volume; this explicit API preserves it.
/// `None` selects both branches, `Some` selects liquid or vapor.
pub fn baseline_pure_spinodal_intersection(
    model: &Cubic,
    z: &[f64],
    property: SpinodalProperty,
    target: f64,
    initial_pressure: f64,
    phase: Option<Phase>,
    liquid_temperature: Option<f64>,
) -> Result<Option<SpinodalPoint>> {
    if !target.is_finite() || phase == Some(Phase::MinimumGibbs) {
        return Err(Error::InvalidInput("spinodal intersection specification"));
    }
    let curve = baseline_pure_spinodal_curve(model, z, initial_pressure, liquid_temperature, None)?;
    locate_on_pure_curve(model, z, &curve, property, target, phase, None)
}

/// Port `locate_spinodal_prop_min_max_pure_fluid`, including the last turning
/// segment, monotone-branch first-point fallback and mixed-coordinate P_min.
pub fn baseline_pure_spinodal_extrema(
    model: &Cubic,
    z: &[f64],
    property: SpinodalProperty,
    initial_pressure: f64,
) -> Result<BaselineSpinodalExtrema> {
    let curve = baseline_pure_spinodal_curve(model, z, initial_pressure, None, None)?;
    extrema_on_pure_curve(model, z, &curve, property, None)
}

/// Immutable reusable pure spinodal grid and all six property arrays. The
/// borrowed model cannot be modified while the search is alive. Preparation
/// is excluded from subsequent queries but remains explicit to the caller.
pub struct PreparedPureSpinodal<'a> {
    model: &'a Cubic,
    composition: Vec<f64>,
    curve: SpinodalCurve,
    properties: [Vec<f64>; 6],
}

impl Cubic {
    pub fn prepare_pure_spinodal_search(
        &self,
        z: &[f64],
        initial_pressure: f64,
        liquid_temperature: Option<f64>,
    ) -> Result<PreparedPureSpinodal<'_>> {
        let curve =
            baseline_pure_spinodal_curve(self, z, initial_pressure, liquid_temperature, None)?;
        let mut properties: [Vec<f64>; 6] = std::array::from_fn(|_| Vec::with_capacity(200));
        for p in &curve.points {
            let c = self.prepare(p.temperature, z)?.caloric_tv(p.volume)?;
            let values = [
                c.entropy,
                p.volume,
                c.internal_energy + c.pressure * p.volume,
                p.temperature,
                c.pressure,
                c.internal_energy,
            ];
            for (buffer, value) in properties.iter_mut().zip(values) {
                buffer.push(value);
            }
        }
        Ok(PreparedPureSpinodal {
            model: self,
            composition: z.to_vec(),
            curve,
            properties,
        })
    }
}

impl PreparedPureSpinodal<'_> {
    pub fn curve(&self) -> &SpinodalCurve {
        &self.curve
    }
    /// Reuse the stored grid/property arrays for source-selected extrema.
    ///
    /// ```
    /// use thermopack_rs::{Cubic, CubicKind};
    /// use thermopack_rs::spinodal_properties::SpinodalProperty;
    /// # fn main() -> Result<(), thermopack_rs::Error> {
    /// let model = Cubic::new(CubicKind::PengRobinson, &["C1"])?;
    /// let search = model.prepare_pure_spinodal_search(&[1.0], 1e6, None)?;
    /// let extrema = search.extrema(SpinodalProperty::Entropy)?;
    /// assert!(!extrema.ordering_ambiguous);
    /// # Ok(())
    /// # }
    /// ```
    pub fn extrema(&self, property: SpinodalProperty) -> Result<BaselineSpinodalExtrema> {
        extrema_on_pure_curve(
            self.model,
            &self.composition,
            &self.curve,
            property,
            Some(&self.properties[property.code() - 1]),
        )
    }
    /// Query the original pure branch/reporting conventions using a cached grid.
    pub fn locate(
        &self,
        property: SpinodalProperty,
        target: f64,
        phase: Option<Phase>,
    ) -> Result<Option<SpinodalPoint>> {
        if !target.is_finite() || phase == Some(Phase::MinimumGibbs) {
            return Err(Error::InvalidInput("spinodal intersection specification"));
        }
        locate_on_pure_curve(
            self.model,
            &self.composition,
            &self.curve,
            property,
            target,
            phase,
            Some(&self.properties[property.code() - 1]),
        )
    }
}

fn extrema_on_pure_curve(
    model: &Cubic,
    z: &[f64],
    curve: &SpinodalCurve,
    property: SpinodalProperty,
    cached: Option<&[f64]>,
) -> Result<BaselineSpinodalExtrema> {
    let mut candidates = [(curve.points[0], 0.); 2];
    for (branch, liquid) in [true, false].into_iter().enumerate() {
        let index = |i| if liquid { i } else { 199 - i };
        let point = |i| curve.points[index(i)];
        let mut values = [0.; 100];
        for (i, value) in values.iter_mut().enumerate() {
            *value = cached.map_or_else(
                || property.evaluate(model, z, point(i)),
                |v| Ok(v[index(i)]),
            )?;
        }
        let turn = (2..100)
            .filter(|&i| (values[i] - values[i - 1]) * (values[i - 1] - values[i - 2]) < 0.)
            .map(|i| i - 1)
            .next_back();
        let candidate = if let Some(j) = turn {
            let p = refine_property(
                model,
                z,
                point(j - 1),
                point((j + 1).min(99)),
                property,
                0.,
                Refinement::SourceGradient,
            )?;
            (p, property.evaluate(model, z, p)?)
        } else {
            (point(0), values[0])
        };
        candidates[branch] = candidate;
    }
    let ordering_ambiguous = property == SpinodalProperty::Pressure
        && (candidates[0].1 - candidates[1].1).abs()
            <= 2. + 1e-7 * candidates[0].1.abs().max(candidates[1].1.abs());
    let ((mut minimum, _), (mut maximum, _)) = if candidates[1].1 > candidates[0].1 {
        (candidates[0], candidates[1])
    } else {
        (candidates[1], candidates[0])
    };
    // Source typo: Ps_min = pressure(Ts_max,vs_min,z).
    let maximum_temperature = model.prepare(maximum.temperature, z)?;
    minimum.pressure = maximum_temperature.pressure(minimum.volume)?.value;
    maximum.pressure = maximum_temperature.pressure(maximum.volume)?.value;
    Ok(BaselineSpinodalExtrema {
        minimum,
        maximum,
        ordering_ambiguous,
    })
}

/// The fork's mixture selection: map from P0 to model Tmin and refine the
/// last strict property sign change. Its uncorrected output-volume convention
/// is preserved, as for the pure compatibility API. No crossing returns None.
pub fn baseline_mixture_spinodal_intersection(
    model: &Cubic,
    z: &[f64],
    property: SpinodalProperty,
    target: f64,
    initial_pressure: f64,
    liquid_temperature: Option<f64>,
) -> Result<Option<SpinodalPoint>> {
    if !target.is_finite() {
        return Err(Error::InvalidInput("finite spinodal target required"));
    }
    if z.iter().filter(|&&x| x > 0.0).count() < 2 {
        return Err(Error::InvalidInput("mixture spinodal composition required"));
    }
    let curve = crate::mixture_spinodal::baseline_trace_spinodal(
        model,
        z,
        crate::mixture_spinodal::SpinodalTraceOptions {
            initial_pressure,
            initial_liquid_temperature: liquid_temperature,
            minimum_vapor_temperature: Some(model.temperature_bounds().0),
            ..Default::default()
        },
    )?;
    if curve.end == SpinodalTraceEnd::PointLimit {
        return Err(Error::NotConverged("incomplete spinodal property curve"));
    }
    let values = curve
        .points
        .iter()
        .map(|&p| property.evaluate(model, z, p))
        .collect::<Result<Vec<_>>>()?;
    let last = values
        .windows(2)
        .rposition(|w| (w[0] - target) * (w[1] - target) < 0.0);
    last.map(|i| {
        refine_property(
            model,
            z,
            curve.points[i],
            curve.points[i + 1],
            property,
            target,
            Refinement::SourceValue,
        )
    })
    .transpose()
}

fn locate_on_pure_curve(
    model: &Cubic,
    z: &[f64],
    curve: &SpinodalCurve,
    property: SpinodalProperty,
    target: f64,
    phase: Option<Phase>,
    cached: Option<&Vec<f64>>,
) -> Result<Option<SpinodalPoint>> {
    let mut found = None;
    for branch in [Phase::Liquid, Phase::Vapor] {
        if phase.is_some_and(|p| p != branch) {
            continue;
        }
        let points: Vec<_> = if branch == Phase::Liquid {
            curve.points[..100].to_vec()
        } else {
            curve.points[100..].iter().rev().copied().collect()
        };
        let mut values = if let Some(cache) = cached {
            if branch == Phase::Liquid {
                cache[..100].to_vec()
            } else {
                cache[100..].iter().rev().copied().collect()
            }
        } else {
            points
                .iter()
                .map(|&p| property.evaluate(model, z, p))
                .collect::<Result<Vec<_>>>()?
        };
        // The source's uncorrected reported volume depends on which adjacent
        // grid interval owns an exact target. Stabilize equality against the
        // independently solved grid's floating-point/iteration error.
        for v in &mut values {
            if (*v - target).abs() < 2e-8 * target.abs().max(1e-12) {
                *v = target;
            }
        }
        let mut turn = None;
        for i in 2..100 {
            if (values[i] - values[i - 1]) * (values[i - 1] - values[i - 2]) < 0.0 {
                turn = Some(i - 1);
            }
        }
        let (mut previous, start) = if let Some(j) = turn {
            let point = refine_property(
                model,
                z,
                points[j - 1],
                points[(j + 1).min(99)],
                property,
                0.0,
                Refinement::SourceGradient,
            )?;
            values[j] = property.evaluate(model, z, point)?;
            (values[j], j + 1)
        } else {
            (values[0], 1)
        };
        if (previous - target) * (values[99] - target) <= 0.0 {
            for i in start..100 {
                if (previous - target) * (values[i] - target) <= 0.0 {
                    found = Some(refine_property(
                        model,
                        z,
                        points[i - 1],
                        points[i],
                        property,
                        target,
                        Refinement::SourceValue,
                    )?);
                    break;
                }
                previous = values[i];
            }
        }
    }
    Ok(found)
}
