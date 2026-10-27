//! Source-probe protocol shared by the parity executable and fixture replay.
use thermopack_rs::MultiparameterIdealMixtureState;

pub fn coordinates(q: &MultiparameterIdealMixtureState) -> Vec<f64> {
    let mut out = vec![
        q.compressibility_factor,
        q.compressibility_factor_dt,
        q.compressibility_factor_dp,
    ];
    out.extend(&q.compressibility_factor_composition_coefficients);
    out.extend([
        q.residual_enthalpy,
        q.residual_enthalpy_dt,
        q.residual_enthalpy_dp,
    ]);
    out.extend(&q.residual_enthalpy_composition_coefficients);
    out.extend([
        q.residual_entropy,
        q.residual_entropy_dt,
        q.residual_entropy_dp,
    ]);
    out.extend(&q.residual_entropy_composition_coefficients);
    out.extend([
        q.residual_gibbs_energy,
        q.residual_gibbs_energy_dt,
        q.baseline_gibbs_dp,
        q.residual_gibbs_energy_dp,
    ]);
    out.extend(&q.residual_gibbs_composition_coefficients);
    out.extend(
        q.fugacity_coefficients
            .as_ref()
            .expect("source protocol requires representable phi"),
    );
    out.extend(&q.ln_fugacity_coefficients_dt);
    out.extend(&q.ln_fugacity_coefficients_dp);
    out.extend(
        q.ln_fugacity_coefficients_composition_derivatives
            .iter()
            .flatten(),
    );
    out.extend(&q.component_volumes);
    out
}
