//! Source-derived fluid database. No runtime parsing or external data files.
use crate::{Error, Result};

#[derive(Debug, Clone, Copy)]
pub struct Component {
    pub ident: &'static str,
    pub name: &'static str,
    pub aliases: &'static [&'static str],
    /// Molecular mass [g/mol].
    pub molar_mass: f64,
    /// Critical temperature (K).
    pub critical_temperature: f64,
    /// Critical pressure (Pa).
    pub critical_pressure: f64,
    pub acentric_factor: f64,
    pub(crate) cp_type: u8,
    pub(crate) cp: [f64; 21],
    pub(crate) reference_enthalpy: f64,
    pub(crate) reference_entropy: f64,
    pub(crate) reference_pressure: f64,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct Interaction {
    pub eos: &'static str,
    pub first: &'static str,
    pub second: &'static str,
    pub value: f64,
    pub covolume: bool,
}
include!("data.rs");

pub fn components() -> &'static [Component] {
    COMPONENTS
}

pub fn component(name: &str) -> Result<&'static Component> {
    let name = name.trim();
    COMPONENTS
        .iter()
        .find(|c| {
            c.ident.eq_ignore_ascii_case(name)
                || c.name.eq_ignore_ascii_case(name)
                || c.aliases.iter().any(|a| a.eq_ignore_ascii_case(name))
        })
        .ok_or_else(|| Error::UnknownComponent(name.to_owned()))
}

pub(crate) fn interaction(eos: &str, first: &str, second: &str, covolume: bool) -> f64 {
    INTERACTIONS
        .iter()
        .rev()
        .find(|k| {
            k.eos == eos
                && k.covolume == covolume
                && ((k.first == first && k.second == second)
                    || (k.first == second && k.second == first))
        })
        .map_or(0.0, |k| k.value)
}
