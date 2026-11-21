//! Radial distribution functions of the association models.
//!
//! The original `src/saft_rdf.f90` routes a request for the radial
//! distribution function of a pair of components to the model of the active
//! thermodynamic model. The native module provides the branch of the cubic
//! association model, which needs no temperature dependence, and the branch of
//! the PC-SAFT model, which follows the hard-sphere model of Boublik.

use crate::Result;
use crate::bmcsl::{self, HsDiameters};
use crate::cpa::{Cpa, PreparedCpa};
use crate::pc_saft::{PcSaft, PreparedPcSaft};

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

/// The hard-sphere contact radial distribution function of the PC-SAFT model.
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
    let total: f64 = moles.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return Err(crate::Error::InvalidInput(
            "the number of moles must be strictly positive",
        ));
    }
    let composition: Vec<f64> = moles.iter().map(|amount| amount / total).collect();
    let prepared = PreparedPcSaft::new(model, temperature, &composition)?;
    let diameters = HsDiameters {
        d: prepared.d.clone(),
        d_t: prepared.d_t.clone(),
        d_tt: prepared.d_tt.clone(),
    };
    let segments: Vec<f64> = model
        .components()
        .iter()
        .map(|component| component.m)
        .collect();
    let packing = bmcsl::packing_fractions(volume, moles, &diameters, &segments);
    let correlation = bmcsl::pair_correlation(&packing, &diameters, first, second);

    let mut derivatives = RdfDerivatives {
        g: correlation.value,
        temperature: correlation.d_t,
        volume: correlation.d_v,
        moles: correlation.d_n.clone(),
        temperature_second: correlation.d_tt,
        temperature_volume: correlation.d_tv,
        temperature_moles: correlation.d_tn.clone(),
        volume_second: correlation.d_vv,
        volume_moles: correlation.d_vn.clone(),
        moles_second: Vec::with_capacity(count * count),
    };
    for row in &correlation.d_nn[..count] {
        derivatives.moles_second.extend_from_slice(&row[..count]);
    }
    Ok(derivatives)
}
