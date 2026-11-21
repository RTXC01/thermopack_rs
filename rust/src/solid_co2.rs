//! Thermodynamic potentials of solid carbon dioxide.
//! Ported from `src/co2_gibbs.f90`, which follows
//! "Equation of State for Solid Carbon Dioxide Based on the Gibbs Free Energy",
//! Jäger and Span, J. Chem. Eng. Data 57 (2012) 590-597.
//! Copyright thermotools and Contributors. Apache-2.0.

/// The gas constant of the original model [J/mol/K].
///
/// The source uses its own value, which differs slightly from the value of
/// the equation of state library.
const R: f64 = 8.314472;
/// Reference temperature of the model [K].
const TREF: f64 = 150.0;
/// Reference pressure of the model [Pa].
const PREF: f64 = 101325.0;
/// Melting enthalpy at the triple point [J/mol].
const DH_MELT_TR: f64 = 8875.0;
/// Mole weight of carbon dioxide [g/mol].
const MWCO2: f64 = 44.01;
/// Exponent of the reduced pressure term.
const N: f64 = 7.0;
/// One minus the reciprocal of the exponent.
const NN: f64 = (N - 1.0) / N;

const G0: f64 = -28.605849252957245;
const G1: f64 = 4.5088732 + 14.500287452234922;
const G2: f64 = -2.0109135;
const G3: f64 = -2.7976237;
const G4: f64 = 2.6427834e-1;
const G5: f64 = 3.8259935;
const G6: f64 = 3.1711996e-1;
const G7: f64 = 2.2087195e-3;
const G8: f64 = -1.1289668;
const G9: f64 = 9.2923982e-3;
const G10: f64 = 3.3914617e3;

const G0A: f64 = 3.9993365e-2;
const G1A: f64 = 2.3945101e-3;
const G2A: f64 = 3.2839467e-1;
const G3A: f64 = 5.7918471e-2;
const G4A: f64 = G1A;
const G5A: f64 = -2.6531689e-3;
const G6A: f64 = 1.6419734e-1;
const G7A: f64 = 1.7594802e-1;
const G8A: f64 = 2.6531689e-3;
const G0K: f64 = 2.2690751e-1;
const G1K: f64 = -7.5019750e-2;
const G2K: f64 = 2.6442913e-1;

/// The Gibbs energy model of solid carbon dioxide.
///
/// The reference levels of the model are shifted so that the solid and the
/// fluid phase of the equation of state share the same reference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolidCo2 {
    /// Reference level of the Gibbs energy.
    g0: f64,
    /// Reference level of the entropy.
    g1: f64,
}

impl Default for SolidCo2 {
    fn default() -> Self {
        Self { g0: G0, g1: G1 }
    }
}

impl SolidCo2 {
    /// The model with the reference levels of the source.
    pub fn new() -> Self {
        Self::default()
    }

    /// Shift the reference levels to the triple point of the equation of state.
    ///
    /// The original `sco2init` takes the liquid entropy and Gibbs energy of
    /// the fluid phase at the triple point together with the triple point
    /// temperature and pressure.
    pub fn initialize(
        &mut self,
        liquid_entropy: f64,
        liquid_gibbs: f64,
        temperature: f64,
        pressure: f64,
    ) {
        let solid_entropy = liquid_entropy - DH_MELT_TR / temperature;
        let current_entropy = self.entropy(temperature, pressure);
        self.g1 -= (solid_entropy - current_entropy) / R;
        let current_gibbs = self.gibbs(temperature, pressure);
        self.g0 += (liquid_gibbs - current_gibbs) / (R * TREF);
    }

    /// The reduced temperature of the model.
    #[inline]
    fn reduced_temperature(temperature: f64) -> f64 {
        temperature / TREF
    }

    /// The reduced pressure of the model.
    #[inline]
    fn reduced_pressure(pressure: f64) -> f64 {
        pressure / PREF
    }

    /// The bulk modulus term of the model.
    #[inline]
    fn k(v: f64) -> f64 {
        G0K * v * v + G1K * v + G2K
    }

    /// The first temperature derivative of the bulk modulus term.
    #[inline]
    fn dk_dv(v: f64) -> f64 {
        2.0 * G0K * v + G1K
    }

    /// The second temperature derivative of the bulk modulus term.
    #[inline]
    fn d2k_dv2(_v: f64) -> f64 {
        2.0 * G0K
    }

    /// The exponential term of the model.
    #[inline]
    fn f_alpha(v: f64) -> f64 {
        G0A * (v * v - 1.0)
            + G1A * ((v * v - G2A * v + G3A) / (1.0 - G2A + G3A)).ln()
            + G4A * ((v * v + G2A * v + G3A) / (1.0 + G2A + G3A)).ln()
            + G5A * (((v - G6A) / G7A).atan() - ((1.0 - G6A) / G7A).atan())
            + G8A * (((v + G6A) / G7A).atan() - ((1.0 + G6A) / G7A).atan())
    }

    /// The first temperature derivative of the exponential term.
    #[inline]
    fn d_f_alpha_dv(v: f64) -> f64 {
        let u1 = (v - G6A) / G7A;
        let u2 = (v + G6A) / G7A;
        2.0 * G0A * v
            + G1A * (2.0 * v - G2A) / (v * v - G2A * v + G3A)
            + G4A * (2.0 * v + G2A) / (v * v + G2A * v + G3A)
            + (G5A / G7A) / (1.0 + u1 * u1)
            + (G8A / G7A) / (1.0 + u2 * u2)
    }

    /// The second temperature derivative of the exponential term.
    #[inline]
    fn d2_f_alpha_dv2(v: f64) -> f64 {
        let q1 = v * v - G2A * v + G3A;
        let d1 = 2.0 * v - G2A;
        let q2 = v * v + G2A * v + G3A;
        let d2 = 2.0 * v + G2A;
        let u1 = (v - G6A) / G7A;
        let p1 = 1.0 + u1 * u1;
        let u2 = (v + G6A) / G7A;
        let p2 = 1.0 + u2 * u2;
        let g7a2 = G7A * G7A;
        2.0 * G0A + G1A * (2.0 * q1 - d1 * d1) / (q1 * q1) + G4A * (2.0 * q2 - d2 * d2) / (q2 * q2)
            - (G5A / g7a2) / (p1 * p1) * 2.0 * (v - G6A) / G7A
            - (G8A / g7a2) / (p2 * p2) * 2.0 * (v + G6A) / G7A
    }

    /// The Gibbs energy of solid carbon dioxide [J/mol].
    pub fn gibbs(&self, temperature: f64, pressure: f64) -> f64 {
        let v = Self::reduced_temperature(temperature);
        let dv = v - 1.0;
        let pi = Self::reduced_pressure(pressure);
        let dpi = pi - 1.0;
        let k = Self::k(v);
        let f = Self::f_alpha(v);
        let g = self.g0
            + self.g1 * dv
            + G2 * dv * dv
            + G3 * (((v * v + G4 * G4) / (1.0 + G4 * G4)).ln()
                - (2.0 * v / G4) * ((v / G4).atan() - (1.0 / G4).atan()))
            + G5 * (((v * v + G6 * G6) / (1.0 + G6 * G6)).ln()
                - (2.0 * v / G6) * ((v / G6).atan() - (1.0 / G6).atan()))
            + G7 * dpi * (f.exp() + k * G8)
            + G9 * k * ((pi + G10).powf(NN) - (1.0 + G10).powf(NN));
        g * R * TREF
    }

    /// The pressure derivative of the Gibbs energy [m3/mol].
    pub fn dg_dp(&self, temperature: f64, pressure: f64) -> f64 {
        let v = Self::reduced_temperature(temperature);
        let pi = Self::reduced_pressure(pressure);
        let k = Self::k(v);
        let f = Self::f_alpha(v);
        (R * TREF / PREF) * (G7 * (f.exp() + k * G8) + G9 * k * NN * (pi + G10).powf(-1.0 / N))
    }

    /// The second pressure derivative of the Gibbs energy [m3/mol/Pa].
    pub fn d2g_dp2(&self, temperature: f64, pressure: f64) -> f64 {
        let v = Self::reduced_temperature(temperature);
        let pi = Self::reduced_pressure(pressure);
        let k = Self::k(v);
        -R * TREF / (PREF * PREF) * G9 * k * NN / N * (pi + G10).powf(-(N + 1.0) / N)
    }

    /// The temperature derivative of the Gibbs energy [J/mol/K].
    pub fn dg_dt(&self, temperature: f64, pressure: f64) -> f64 {
        let v = Self::reduced_temperature(temperature);
        let dv = v - 1.0;
        let pi = Self::reduced_pressure(pressure);
        let dpi = pi - 1.0;
        let ef = Self::f_alpha(v).exp();
        let dfdv = Self::d_f_alpha_dv(v);
        let dkdv = Self::dk_dv(v);
        R * (self.g1 + G2 * 2.0 * dv
            - 2.0 * G3 / G4 * ((v / G4).atan() - (1.0 / G4).atan())
            - 2.0 * G5 / G6 * ((v / G6).atan() - (1.0 / G6).atan())
            + G7 * dpi * (ef * dfdv + G8 * dkdv)
            + G9 * dkdv * ((pi + G10).powf(NN) - (1.0 + G10).powf(NN)))
    }

    /// The second temperature derivative of the Gibbs energy [J/mol/K2].
    pub fn d2g_dt2(&self, temperature: f64, pressure: f64) -> f64 {
        let v = Self::reduced_temperature(temperature);
        let pi = Self::reduced_pressure(pressure);
        let dpi = pi - 1.0;
        let ef = Self::f_alpha(v).exp();
        let dfdv = Self::d_f_alpha_dv(v);
        let d2fdv2 = Self::d2_f_alpha_dv2(v);
        let d2kdv2 = Self::d2k_dv2(v);
        R / TREF
            * (2.0 * G2 - 2.0 * G3 / (v * v + G4 * G4) - 2.0 * G5 / (v * v + G6 * G6)
                + G7 * dpi * (ef * d2fdv2 + ef * dfdv * dfdv + G8 * d2kdv2)
                + G9 * d2kdv2 * ((pi + G10).powf(NN) - (1.0 + G10).powf(NN)))
    }

    /// The mixed temperature and pressure derivative of the Gibbs energy
    /// [J/mol/K/Pa].
    pub fn d2g_dtdp(&self, temperature: f64, pressure: f64) -> f64 {
        let v = Self::reduced_temperature(temperature);
        let pi = Self::reduced_pressure(pressure);
        let ef = Self::f_alpha(v).exp();
        let dfdv = Self::d_f_alpha_dv(v);
        let dkdv = Self::dk_dv(v);
        R / PREF * (G7 * (ef * dfdv + G8 * dkdv) + G9 * dkdv * NN * (pi + G10).powf(-1.0 / N))
    }

    /// The specific volume of solid carbon dioxide [m3/mol].
    pub fn specific_volume(&self, temperature: f64, pressure: f64) -> f64 {
        self.dg_dp(temperature, pressure)
    }

    /// The specific entropy of solid carbon dioxide [J/mol/K].
    pub fn entropy(&self, temperature: f64, pressure: f64) -> f64 {
        -self.dg_dt(temperature, pressure)
    }

    /// The specific enthalpy of solid carbon dioxide [J/mol].
    pub fn enthalpy(&self, temperature: f64, pressure: f64) -> f64 {
        self.gibbs(temperature, pressure) + temperature * self.entropy(temperature, pressure)
    }

    /// The specific internal energy of solid carbon dioxide [J/mol].
    pub fn internal_energy(&self, temperature: f64, pressure: f64) -> f64 {
        self.enthalpy(temperature, pressure)
            - pressure * self.specific_volume(temperature, pressure)
    }

    /// The specific Helmholtz energy of solid carbon dioxide [J/mol].
    pub fn helmholtz(&self, temperature: f64, pressure: f64) -> f64 {
        self.gibbs(temperature, pressure) - pressure * self.specific_volume(temperature, pressure)
    }

    /// The specific heat capacity of solid carbon dioxide [J/mol/K].
    pub fn heat_capacity(&self, temperature: f64, pressure: f64) -> f64 {
        -temperature * self.d2g_dt2(temperature, pressure)
    }

    /// The speed of sound in solid carbon dioxide [m/s].
    pub fn speed_of_sound(&self, temperature: f64, pressure: f64) -> f64 {
        let volume = self.specific_volume(temperature, pressure);
        let dsdt = -self.d2g_dt2(temperature, pressure);
        let dsdp = -self.d2g_dtdp(temperature, pressure);
        let dtdp = -dsdp / dsdt;
        let dvdp = self.d2g_dp2(temperature, pressure);
        let dvdt = self.d2g_dtdp(temperature, pressure);
        let denominator = dvdp + dvdt * dtdp;
        volume / (-denominator / MWCO2 / 1.0e3).sqrt()
    }
}
