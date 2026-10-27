//! Qualification-only coordinate order; source expectations are independent.
use thermopack_rs::{MbwrSourceTpProperties, MbwrTpProperties, MbwrTpQuantity};
fn quantity(q: MbwrTpQuantity, out: &mut Vec<f64>) {
    out.push(q.value);
    if let Some(d) = q.derivatives {
        out.extend([1., d.dt, d.dp, d.dn]);
    } else {
        out.extend([0.; 4]);
    }
}
pub fn physical(q: &MbwrTpProperties) -> Vec<f64> {
    let mut out = vec![q.pressure];
    for p in [
        q.volume,
        q.internal_energy,
        q.enthalpy,
        q.entropy,
        q.helmholtz_energy,
        q.gibbs_energy,
        q.chemical_potential,
        q.compressibility,
        q.ln_fugacity_coefficient,
    ] {
        quantity(p, &mut out);
    }
    out.push(q.heat_capacity_v);
    for p in [q.heat_capacity_p, q.speed_of_sound] {
        out.extend([f64::from(p.is_some()), p.unwrap_or(0.)]);
    }
    assert_eq!(out.len(), 51);
    out
}
pub fn source(q: &MbwrSourceTpProperties) -> Vec<f64> {
    let mut out = vec![q.pressure, q.volume];
    for p in [
        q.enthalpy,
        q.entropy,
        q.compressibility,
        q.ln_fugacity_coefficient,
    ] {
        quantity(p, &mut out);
    }
    assert_eq!(out.len(), 22);
    out
}
