//! Qualification-only fixed coordinate order, independent source expectations.
use thermopack_rs::{MbwrTvProperties, MbwrTvQuantity, PureHelmholtz};
fn jet(f: PureHelmholtz) -> [f64; 10] {
    [
        f.value, f.dt, f.dv, f.dn, f.dtt, f.dtv, f.dtn, f.dvv, f.dvn, f.dnn,
    ]
}
fn quantity(q: MbwrTvQuantity) -> [f64; 4] {
    [q.value, q.dt, q.dv, q.dn]
}
pub fn coordinates(q: &MbwrTvProperties) -> Vec<f64> {
    let mut out = jet(q.reduced_helmholtz).to_vec();
    let p = q.pressure;
    out.extend([p.value, p.dt, p.dv, p.dn]);
    for prop in [
        q.internal_energy,
        q.enthalpy,
        q.entropy,
        q.helmholtz_energy,
        q.gibbs_energy,
        q.chemical_potential,
    ] {
        out.extend(quantity(prop));
    }
    out.push(q.heat_capacity_v);
    for value in [
        q.heat_capacity_p,
        q.ln_fugacity_coefficient,
        q.speed_of_sound,
    ] {
        out.extend([f64::from(value.is_some()), value.unwrap_or(0.)]);
    }
    out
}
