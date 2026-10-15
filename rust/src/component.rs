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
    pub rackett_factor: f64,
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
    pub reference: &'static str,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct TranslationRecord {
    pub eos: &'static str,
    pub component: &'static str,
    pub reference: &'static str,
    pub correlation: u8,
    pub coefficients: [f64; 6],
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct IdealRecord {
    pub component: &'static str,
    pub reference: &'static str,
    pub correlation: u8,
    pub coefficients: [f64; 21],
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct AlphaRecord {
    pub eos: &'static str,
    pub component: &'static str,
    pub reference: &'static str,
    pub mathias_copeman: bool,
    pub parameters: [f64; 3],
}
include!("data.rs");

pub fn components() -> &'static [Component] {
    COMPONENTS
}

pub(crate) fn alpha_records() -> &'static [AlphaRecord] {
    ALPHA_RECORDS
}

pub(crate) fn translation_record(
    eos: &str,
    component: &str,
    reference: &str,
) -> Option<&'static TranslationRecord> {
    select_record(TRANSLATION_RECORDS, reference, |r| {
        (r.eos.eq_ignore_ascii_case(eos) && r.component == component).then_some(r.reference)
    })
}

pub(crate) fn ideal_component(component: &Component) -> Component {
    let mut result = *component;
    // init_cubic passes paramref_loc(len=100) to isRef, whose database ref
    // has length 40. Fortran INDEX cannot match the longer substring, so
    // getCpDBindex always selects its last case-sensitive Default fallback.
    if let Some(record) = IDEAL_RECORDS
        .iter()
        .rev()
        .find(|r| r.component == component.ident && r.reference.contains("Default"))
    {
        result.cp_type = record.correlation;
        result.cp = record.coefficients;
    }
    result
}

// First exact non-DEFAULT token at the best caller priority; otherwise last
// default. The substring default check follows parameters.f90/stringmod.f90.
pub(crate) fn select_record<'a, T>(
    records: &'a [T],
    reference: &str,
    classify: impl Fn(&T) -> Option<&str>,
) -> Option<&'a T> {
    let (mut best, mut fallback, mut rank) = (None, None, usize::MAX);
    for record in records {
        let Some(label) = classify(record) else {
            continue;
        };
        let matching = (!reference.trim().is_empty() && !label.trim().is_empty())
            .then(|| {
                reference
                    .trim_end()
                    .split_terminator('/')
                    .map(str::trim)
                    .filter(|s| !s.eq_ignore_ascii_case("DEFAULT"))
                    .position(|s| {
                        label
                            .trim_end()
                            .split_terminator('/')
                            .any(|r| s.eq_ignore_ascii_case(r.trim()))
                    })
            })
            .flatten();
        if let Some(position) = matching {
            if position < rank {
                best = Some(record);
                rank = position;
            }
        } else if label.to_ascii_uppercase().contains("DEFAULT") {
            fallback = Some(record);
        }
    }
    best.or(fallback)
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

pub(crate) fn interaction_with_reference(
    eos: &str,
    first: &str,
    second: &str,
    covolume: bool,
    reference: &str,
) -> f64 {
    select_record(INTERACTIONS, reference, |k| {
        (k.eos == eos
            && k.covolume == covolume
            && ((k.first == first && k.second == second)
                || (k.first == second && k.second == first)))
            .then_some(k.reference)
    })
    .map_or(0., |k| k.value)
}
