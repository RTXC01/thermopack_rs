//! Radial distribution functions of the association models.
//!
//! The original `src/saft_rdf.f90` routes a request for the radial
//! distribution function of a pair of components to the model of the active
//! thermodynamic model. The native module provides the branch of the cubic
//! association model, which needs no temperature dependence, the branch of the
//! original PC-SAFT model (Boublik-Mansoori-Carnahan-Starling-Leland), and the
//! branch of the simplified PC-SAFT model (Carnahan-Starling average packing
//! fraction).

use crate::Result;
use crate::bmcsl::{self, HsDiameters};
use crate::cpa::{Cpa, PreparedCpa};
use crate::pc_saft::{N_AVOGADRO, PcSaft, PreparedPcSaft};

/// The radial distribution function and its derivatives.
///
/// The derivatives follow the original: the temperature derivatives of the
/// hard-sphere contact function of the cubic association model vanish.
#[derive(Debug, Clone, Default)]
pub struct RdfDerivatives {
    /// The value of the radial distribution function.
    pub g: f64,
    /// The derivative with respect to the temperature.
    pub temperature: f64,
    /// The derivative with respect to the volume.
    pub volume: f64,
    /// The derivatives with respect to the number of moles.
    pub moles: Vec<f64>,
    /// The second derivative with respect to the temperature.
    pub temperature_second: f64,
    /// The mixed derivative with respect to the temperature and the volume.
    pub temperature_volume: f64,
    /// The mixed derivatives with respect to the temperature and the number
    /// of moles.
    pub temperature_moles: Vec<f64>,
    /// The second derivative with respect to the volume.
    pub volume_second: f64,
    /// The mixed derivatives with respect to the volume and the number of
    /// moles.
    pub volume_moles: Vec<f64>,
    /// The second derivatives with respect to the number of moles.
    pub moles_second: Vec<f64>,
}

/// The hard-sphere contact radial distribution function of the cubic
/// association model.
///
/// The volume is in cubic metre and the number of moles is the amount of each
/// component. The function reproduces the original `g_rdf_cpa`, including the
/// packing fraction of the mixture and its derivatives.
pub fn cpa_rdf(
    model: &Cpa,
    temperature: f64,
    volume: f64,
    moles: &[f64],
) -> Result<RdfDerivatives> {
    crate::error::positive(temperature, "the temperature")?;
    crate::error::positive(volume, "the volume")?;
    let count = model.num_components();
    if moles.len() != count {
        return Err(crate::Error::InvalidInput(
            "the number of moles must match the components of the model",
        ));
    }
    let total: f64 = moles.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(crate::Error::InvalidInput(
            "the number of moles must be strictly positive",
        ));
    }
    let composition: Vec<f64> = moles.iter().map(|amount| amount / total).collect();
    let prepared = PreparedCpa::new(model, temperature, &composition)?;
    // The mixture covolume and the covolume of every component.
    let mixture_covolume = prepared.b_mix;
    let covolumes: Vec<f64> = model
        .components
        .iter()
        .map(|component| component.b)
        .collect();
    // The packing fraction of the mixture and its derivatives.
    let total_covolume = total * mixture_covolume;
    let inverse_volume = 1.0 / volume;
    let inverse_volume_squared = inverse_volume * inverse_volume;
    let inverse_volume_cubed = inverse_volume_squared * inverse_volume;
    let eta_volume = -total_covolume * inverse_volume_squared / 4.0;
    let eta_moles: Vec<f64> = covolumes
        .iter()
        .map(|covolume| covolume * inverse_volume / 4.0)
        .collect();
    let eta_volume_second = total_covolume * inverse_volume_cubed / 2.0;
    let eta_volume_moles: Vec<f64> = covolumes
        .iter()
        .map(|covolume| -covolume * inverse_volume_squared / 4.0)
        .collect();
    // The radial distribution function and the derivatives with respect to the
    // packing fraction.
    let density = total / volume;
    let (g, g_eta, g_eta_second) = prepared.g_rdf(density);
    let mut derivatives = RdfDerivatives {
        g,
        temperature: 0.0,
        volume: g_eta * eta_volume,
        moles: eta_moles.iter().map(|eta_n| g_eta * eta_n).collect(),
        temperature_second: 0.0,
        temperature_volume: 0.0,
        temperature_moles: vec![0.0; count],
        volume_second: g_eta_second * eta_volume * eta_volume + g_eta * eta_volume_second,
        volume_moles: Vec::with_capacity(count),
        moles_second: vec![0.0; count * count],
    };
    for i in 0..count {
        derivatives
            .volume_moles
            .push(g_eta_second * eta_volume * eta_moles[i] + g_eta * eta_volume_moles[i]);
        for j in 0..count {
            derivatives.moles_second[i * count + j] = g_eta_second * eta_moles[i] * eta_moles[j];
        }
    }
    Ok(derivatives)
}

/// Evaluates the BMCSL contact RDF at unit total moles `(v_unit = V/N, x = n/N)`.
pub(crate) fn pcsaft_rdf_unit(
    prepared: &PreparedPcSaft<'_>,
    v_unit: f64,
    first: usize,
    second: usize,
) -> RdfDerivatives {
    let count = prepared.model.num_components();
    let diameters = HsDiameters {
        d: prepared.d.clone(),
        d_t: prepared.d_t.clone(),
        d_tt: prepared.d_tt.clone(),
    };
    let segments: Vec<f64> = prepared
        .model
        .components()
        .iter()
        .map(|component| component.m)
        .collect();
    let packing = bmcsl::packing_fractions(v_unit, &prepared.mole_fractions, &diameters, &segments);
    let correlation = bmcsl::pair_correlation(&packing, &diameters, first, second);

    let mut moles_second = Vec::with_capacity(count * count);
    for row in &correlation.d_nn[..count] {
        moles_second.extend_from_slice(&row[..count]);
    }
    RdfDerivatives {
        g: correlation.value,
        temperature: correlation.d_t,
        volume: correlation.d_v,
        moles: correlation.d_n,
        temperature_second: correlation.d_tt,
        temperature_volume: correlation.d_tv,
        temperature_moles: correlation.d_tn,
        volume_second: correlation.d_vv,
        volume_moles: correlation.d_vn,
        moles_second,
    }
}

/// Evaluates the simplified PC-SAFT (Carnahan-Starling) contact RDF at unit
/// total moles `(v_unit = V/N, x = n/N)`.
pub(crate) fn spcsaft_rdf_unit(prepared: &PreparedPcSaft<'_>, v_unit: f64) -> RdfDerivatives {
    let count = prepared.model.num_components();
    let inv_v = 1.0 / v_unit;
    let pref = (std::f64::consts::PI * N_AVOGADRO / 6.0) * inv_v;
    let mut eta = 0.0;
    let mut eta_t = 0.0;
    let mut eta_tt = 0.0;
    let mut eta_x = Vec::with_capacity(count);
    let mut eta_tx = Vec::with_capacity(count);
    let mut eta_vx = Vec::with_capacity(count);

    for i in 0..count {
        let m = prepared.model.components()[i].m;
        let d = prepared.d[i];
        let dt = prepared.d_t[i];
        let dtt = prepared.d_tt[i];
        let d2 = d * d;
        let d3 = d2 * d;
        let xi = prepared.mole_fractions[i];
        let ex = pref * m * d3;
        let etx = pref * m * (3.0 * d2 * dt);
        let ettx = pref * m * (6.0 * d * dt * dt + 3.0 * d2 * dtt);
        eta += xi * ex;
        eta_t += xi * etx;
        eta_tt += xi * ettx;
        eta_x.push(ex);
        eta_tx.push(etx);
        eta_vx.push(-ex * inv_v);
    }

    let eta_v = -eta * inv_v;
    let eta_vv = 2.0 * eta * inv_v * inv_v;
    let eta_tv = -eta_t * inv_v;

    let om = 1.0 - eta;
    let om2 = om * om;
    let om3 = om2 * om;
    let om4 = om2 * om2;
    let om5 = om4 * om;
    let g = (1.0 - 0.5 * eta) / om3;
    let g_eta = (5.0 - 2.0 * eta) / (2.0 * om4);
    let g_etaeta = 3.0 * (3.0 - eta) / om5;

    let mut moles = Vec::with_capacity(count);
    let mut temperature_moles = Vec::with_capacity(count);
    let mut volume_moles = Vec::with_capacity(count);
    let mut moles_second = vec![0.0; count * count];

    for i in 0..count {
        moles.push(g_eta * eta_x[i]);
        temperature_moles.push(g_etaeta * eta_t * eta_x[i] + g_eta * eta_tx[i]);
        volume_moles.push(g_etaeta * eta_v * eta_x[i] + g_eta * eta_vx[i]);
        for j in 0..count {
            moles_second[i * count + j] = g_etaeta * eta_x[i] * eta_x[j];
        }
    }

    RdfDerivatives {
        g,
        temperature: g_eta * eta_t,
        volume: g_eta * eta_v,
        moles,
        temperature_second: g_etaeta * eta_t * eta_t + g_eta * eta_tt,
        temperature_volume: g_etaeta * eta_t * eta_v + g_eta * eta_tv,
        temperature_moles,
        volume_second: g_etaeta * eta_v * eta_v + g_eta * eta_vv,
        volume_moles,
        moles_second,
    }
}

fn scale_rdf_from_unit(mut unit: RdfDerivatives, total: f64) -> Result<RdfDerivatives> {
    let inv_n = 1.0 / total;
    let inv_n2 = inv_n * inv_n;
    unit.volume *= inv_n;
    unit.temperature_volume *= inv_n;
    unit.volume_second *= inv_n2;
    for v in &mut unit.moles {
        *v *= inv_n;
    }
    for v in &mut unit.temperature_moles {
        *v *= inv_n;
    }
    for v in &mut unit.volume_moles {
        *v *= inv_n2;
    }
    for v in &mut unit.moles_second {
        *v *= inv_n2;
    }
    let all_finite = unit.g.is_finite()
        && unit.temperature.is_finite()
        && unit.volume.is_finite()
        && unit.temperature_second.is_finite()
        && unit.temperature_volume.is_finite()
        && unit.volume_second.is_finite()
        && unit.moles.iter().all(|v| v.is_finite())
        && unit.temperature_moles.iter().all(|v| v.is_finite())
        && unit.volume_moles.iter().all(|v| v.is_finite())
        && unit.moles_second.iter().all(|v| v.is_finite());
    if !all_finite {
        return Err(crate::Error::InvalidInput(
            "non-finite PC-SAFT radial distribution function or derivative",
        ));
    }
    Ok(unit)
}

/// The hard-sphere contact radial distribution function of the original PC-SAFT
/// model.
///
/// The volume is in cubic metre and the number of moles is the amount of each
/// component. The function reproduces the original `g_pc_saft_TVn`, which
/// applies the chain rule of the hard-sphere model of Boublik to the packing
/// fractions of the mixture.
pub fn pcsaft_rdf(
    model: &PcSaft,
    temperature: f64,
    volume: f64,
    moles: &[f64],
    first: usize,
    second: usize,
) -> Result<RdfDerivatives> {
    crate::error::positive(temperature, "the temperature")?;
    crate::error::positive(volume, "the volume")?;
    let count = model.num_components();
    if moles.len() != count {
        return Err(crate::Error::InvalidInput(
            "the number of moles must match the components of the model",
        ));
    }
    if first >= count || second >= count {
        return Err(crate::Error::InvalidInput(
            "the components of the pair must match the components of the model",
        ));
    }
    if moles.iter().any(|&v| !v.is_finite() || v < 0.0) {
        return Err(crate::Error::InvalidInput(
            "mole numbers must be finite and non-negative",
        ));
    }
    let total: f64 = moles.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(crate::Error::InvalidInput(
            "the number of moles must be strictly positive",
        ));
    }
    let v_unit = volume / total;
    if !v_unit.is_finite() || v_unit <= 0.0 {
        return Err(crate::Error::InvalidInput(
            "molar volume must be finite and strictly positive",
        ));
    }
    let composition: Vec<f64> = moles.iter().map(|amount| amount / total).collect();
    let prepared = PreparedPcSaft::new(model, temperature, &composition)?;
    let unit = pcsaft_rdf_unit(&prepared, v_unit, first, second);
    scale_rdf_from_unit(unit, total)
}

/// The hard-sphere contact radial distribution function of the simplified
/// PC-SAFT model (`sPC-SAFT` and `sPCP-SAFT`).
///
/// Reproduces the original `g_spc_saft_TVn`, using the Carnahan-Starling
/// average packing fraction. The result is independent of component pair
/// indices.
pub fn spcsaft_rdf(
    model: &PcSaft,
    temperature: f64,
    volume: f64,
    moles: &[f64],
) -> Result<RdfDerivatives> {
    crate::error::positive(temperature, "the temperature")?;
    crate::error::positive(volume, "the volume")?;
    let count = model.num_components();
    if moles.len() != count {
        return Err(crate::Error::InvalidInput(
            "the number of moles must match the components of the model",
        ));
    }
    if moles.iter().any(|&v| !v.is_finite() || v < 0.0) {
        return Err(crate::Error::InvalidInput(
            "mole numbers must be finite and non-negative",
        ));
    }
    let total: f64 = moles.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(crate::Error::InvalidInput(
            "the number of moles must be strictly positive",
        ));
    }
    let v_unit = volume / total;
    if !v_unit.is_finite() || v_unit <= 0.0 {
        return Err(crate::Error::InvalidInput(
            "molar volume must be finite and strictly positive",
        ));
    }
    let composition: Vec<f64> = moles.iter().map(|amount| amount / total).collect();
    let prepared = PreparedPcSaft::new(model, temperature, &composition)?;
    let unit = spcsaft_rdf_unit(&prepared, v_unit);
    scale_rdf_from_unit(unit, total)
}

/// Dispatches to [`pcsaft_rdf`] for original PC-SAFT / PCP-SAFT models or
/// [`spcsaft_rdf`] for simplified sPC-SAFT / sPCP-SAFT models.
pub fn pcsaft_rdf_for_model(
    model: &PcSaft,
    temperature: f64,
    volume: f64,
    moles: &[f64],
    first: usize,
    second: usize,
) -> Result<RdfDerivatives> {
    if model.kind().is_simplified() {
        let count = model.num_components();
        if first >= count || second >= count {
            return Err(crate::Error::InvalidInput(
                "the components of the pair must match the components of the model",
            ));
        }
        spcsaft_rdf(model, temperature, volume, moles)
    } else {
        pcsaft_rdf(model, temperature, volume, moles, first, second)
    }
}
