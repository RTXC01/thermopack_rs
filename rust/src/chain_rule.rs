//! The chain rule of the packing fraction and the reduced distance.
//! Ported from `src/saftvrmie_utils.f90`.
//! Copyright thermotools and Contributors. Apache-2.0.
//!
//! The models of the source express their terms in the packing fraction and in
//! the reduced distance of the hard-sphere centres. The chain rule turns the
//! derivatives with respect to those two variables into the derivatives with
//! respect to the temperature, the volume and the number of moles.

/// The derivatives of a packing fraction with respect to the state variables.
#[derive(Debug, Clone)]
pub struct PackingFractionDerivatives {
    /// The packing fraction.
    pub value: f64,
    /// The derivative with respect to the temperature.
    pub temperature: f64,
    /// The second derivative with respect to the temperature.
    pub temperature_second: f64,
    /// The derivative with respect to the volume.
    pub volume: f64,
    /// The second derivative with respect to the volume.
    pub volume_second: f64,
    /// The third derivative with respect to the volume.
    pub volume_third: f64,
    /// The mixed derivative with respect to the temperature and the volume.
    pub temperature_volume: f64,
    /// The second derivative with respect to the volume and the temperature.
    pub volume_temperature_second: f64,
    /// The second derivative with respect to the volume and the temperature.
    pub volume_temperature_third: f64,
    /// The derivatives with respect to the number of moles.
    pub moles: Vec<f64>,
    /// The second derivatives with respect to the number of moles.
    pub moles_second: Vec<f64>,
    /// The mixed derivatives with respect to the temperature and the number of
    /// moles.
    pub temperature_moles: Vec<f64>,
    /// The mixed derivatives with respect to the volume and the number of
    /// moles.
    pub volume_moles: Vec<f64>,
    /// The second derivatives with respect to the volume and the number of
    /// moles.
    pub volume_moles_second: Vec<f64>,
    /// The mixed derivatives with respect to the volume, the temperature and
    /// the number of moles.
    pub volume_temperature_moles: Vec<f64>,
}

/// The value of a term and its derivatives with respect to the packing
/// fraction and the reduced distance.
#[derive(Debug, Clone, Copy)]
pub struct TermDerivatives {
    /// The value of the term.
    pub value: f64,
    /// The derivative with respect to the packing fraction.
    pub packing_fraction: f64,
    /// The derivative with respect to the reduced distance.
    pub reduced_distance: f64,
    /// The second derivative with respect to the packing fraction.
    pub packing_fraction_second: f64,
    /// The second derivative with respect to the reduced distance.
    pub reduced_distance_second: f64,
    /// The mixed derivative with respect to the packing fraction and the
    /// reduced distance.
    pub packing_fraction_reduced_distance: f64,
    /// The third derivative with respect to the packing fraction.
    pub packing_fraction_third: f64,
    /// The second derivative with respect to the packing fraction and the
    /// reduced distance.
    pub packing_fraction_second_reduced_distance: f64,
    /// The third derivative with respect to the packing fraction and the
    /// reduced distance.
    pub packing_fraction_reduced_distance_second: f64,
}

/// The derivatives of a term with respect to the state variables.
#[derive(Debug, Clone)]
pub struct StateDerivatives {
    /// The value of the term.
    pub value: f64,
    /// The derivative with respect to the temperature.
    pub temperature: f64,
    /// The derivative with respect to the volume.
    pub volume: f64,
    /// The derivatives with respect to the number of moles.
    pub moles: Vec<f64>,
    /// The second derivative with respect to the temperature.
    pub temperature_second: f64,
    /// The second derivative with respect to the volume.
    pub volume_second: f64,
    /// The mixed derivative with respect to the temperature and the volume.
    pub temperature_volume: f64,
    /// The mixed derivatives with respect to the temperature and the number of
    /// moles.
    pub temperature_moles: Vec<f64>,
    /// The mixed derivatives with respect to the volume and the number of
    /// moles.
    pub volume_moles: Vec<f64>,
    /// The second derivatives with respect to the number of moles.
    pub moles_second: Vec<f64>,
    /// The third derivative with respect to the volume.
    pub volume_third: f64,
    /// The second derivative with respect to the volume and the temperature.
    pub volume_temperature_second: f64,
    /// The second derivative with respect to the volume and the temperature.
    pub volume_temperature_third: f64,
    /// The second derivatives with respect to the volume and the number of
    /// moles.
    pub volume_moles_second: Vec<f64>,
    /// The mixed derivatives with respect to the volume, the temperature and
    /// the number of moles.
    pub volume_temperature_moles: Vec<f64>,
}

/// The chain rule of the source.
///
/// The term carries the derivatives with respect to the packing fraction and
/// the reduced distance, the packing fraction carries its derivatives with
/// respect to the state variables, and the reduced distance carries its first
/// and second temperature derivatives. The reduced distance itself does not
/// enter the chain rule, as in the source. The level selects the order of the
/// derivatives, the first level gives the first derivatives, the second level
/// the second derivatives and the third level the third derivatives.
pub fn convert_zeta_x_to_tvn(
    term: &TermDerivatives,
    _reduced_distance: f64,
    reduced_distance_t: f64,
    reduced_distance_tt: f64,
    zeta: &PackingFractionDerivatives,
    level: usize,
) -> StateDerivatives {
    let count = zeta.moles.len();
    let mut result = StateDerivatives {
        value: term.value,
        temperature: 0.0,
        volume: 0.0,
        moles: vec![0.0; count],
        temperature_second: 0.0,
        volume_second: 0.0,
        temperature_volume: 0.0,
        temperature_moles: vec![0.0; count],
        volume_moles: vec![0.0; count],
        moles_second: vec![0.0; count * count],
        volume_third: 0.0,
        volume_temperature_second: 0.0,
        volume_temperature_third: 0.0,
        volume_moles_second: vec![0.0; count * count],
        volume_temperature_moles: vec![0.0; count],
    };
    result.temperature =
        term.packing_fraction * zeta.temperature + term.reduced_distance * reduced_distance_t;
    result.volume = term.packing_fraction * zeta.volume;
    for (index, value) in zeta.moles.iter().enumerate() {
        result.moles[index] = term.packing_fraction * value;
    }
    if level > 1 {
        result.temperature_second = term.packing_fraction * zeta.temperature_second
            + term.reduced_distance * reduced_distance_tt
            + term.packing_fraction_second * zeta.temperature * zeta.temperature
            + term.reduced_distance_second * reduced_distance_t * reduced_distance_t
            + 2.0 * term.packing_fraction_reduced_distance * reduced_distance_t * zeta.temperature;
        result.volume_second = term.packing_fraction * zeta.volume_second
            + term.packing_fraction_second * zeta.volume * zeta.volume;
        result.temperature_volume = term.packing_fraction_second * zeta.temperature * zeta.volume
            + term.packing_fraction_reduced_distance * reduced_distance_t * zeta.volume
            + term.packing_fraction * zeta.temperature_volume;
        for (index, value) in zeta.moles.iter().enumerate() {
            result.temperature_moles[index] =
                term.packing_fraction_second * zeta.temperature * value
                    + term.packing_fraction_reduced_distance * reduced_distance_t * value
                    + term.packing_fraction * zeta.temperature_moles[index];
            result.volume_moles[index] = term.packing_fraction_second * zeta.volume * value
                + term.packing_fraction * zeta.volume_moles[index];
        }
        for index in 0..count {
            for other in 0..count {
                result.moles_second[index * count + other] =
                    term.packing_fraction_second * zeta.moles[index] * zeta.moles[other]
                        + term.packing_fraction * zeta.moles_second[index * count + other];
            }
        }
    }
    if level > 2 {
        result.volume_third = term.packing_fraction_third * zeta.volume * zeta.volume * zeta.volume
            + 3.0 * term.packing_fraction_second * zeta.volume * zeta.volume_second
            + term.packing_fraction * zeta.volume_third;
        result.volume_temperature_second =
            term.packing_fraction_third * zeta.volume * zeta.volume * zeta.temperature
                + term.packing_fraction_second * zeta.volume_second * zeta.temperature
                + term.packing_fraction_reduced_distance * reduced_distance_t * zeta.volume_second
                + term.packing_fraction_second_reduced_distance
                    * reduced_distance_t
                    * zeta.volume
                    * zeta.volume
                + term.packing_fraction * zeta.volume_temperature_second
                + 2.0 * term.packing_fraction_second * zeta.volume * zeta.temperature_volume;
        result.volume_temperature_third =
            term.packing_fraction_second * zeta.volume * zeta.temperature_second
                + term.packing_fraction_reduced_distance * zeta.volume * reduced_distance_tt
                + term.packing_fraction_third * zeta.volume * zeta.temperature * zeta.temperature
                + term.packing_fraction_reduced_distance_second
                    * zeta.volume
                    * reduced_distance_t
                    * reduced_distance_t
                + term.packing_fraction * zeta.volume_temperature_third
                + 2.0 * term.packing_fraction_second * zeta.temperature * zeta.temperature_volume
                + 2.0
                    * term.packing_fraction_second_reduced_distance
                    * zeta.volume
                    * reduced_distance_t
                    * zeta.temperature
                + 2.0
                    * term.packing_fraction_reduced_distance
                    * zeta.temperature_volume
                    * reduced_distance_t;
        for (index, value) in zeta.moles.iter().enumerate() {
            result.volume_temperature_moles[index] = term.packing_fraction_third
                * zeta.volume
                * zeta.temperature
                * value
                + term.packing_fraction_second * zeta.temperature * zeta.volume_moles[index]
                + term.packing_fraction_second * zeta.temperature_volume * value
                + term.packing_fraction_reduced_distance * reduced_distance_t * zeta.volume * value
                + term.packing_fraction_reduced_distance
                    * zeta.volume_moles[index]
                    * reduced_distance_t
                + term.packing_fraction_second * zeta.volume * zeta.temperature_moles[index]
                + term.packing_fraction * zeta.volume_temperature_moles[index];
        }
        for index in 0..count {
            for other in 0..count {
                result.volume_moles_second[index * count + other] = term.packing_fraction_third
                    * zeta.volume
                    * zeta.moles[index]
                    * zeta.moles[other]
                    + 2.0
                        * term.packing_fraction_second
                        * zeta.volume_moles[index]
                        * zeta.moles[other]
                    + 2.0
                        * term.packing_fraction_second
                        * zeta.moles[index]
                        * zeta.volume_moles[other]
                    + term.packing_fraction_reduced_distance
                        * reduced_distance_t
                        * zeta.volume
                        * zeta.moles[index]
                        * zeta.moles[other]
                    + term.packing_fraction * zeta.volume_moles_second[index * count + other];
            }
        }
    }
    result
}
