// Line-oriented reference protocol: KIND|COMP,COMP|T|P|x,x|L or V
use std::io::{self, BufRead};
use thermopack_rs::{Cubic, CubicKind, Phase};

fn translated_model(
    kind: CubicKind,
    components: &str,
    config: &[&str],
) -> thermopack_rs::Result<Cubic> {
    let mut model = if config.contains(&"qhelper") {
        Cubic::quantum(&components.split(',').collect::<Vec<_>>())?
    } else {
        Cubic::with_options(
            kind,
            &components.split(',').collect::<Vec<_>>(),
            thermopack_rs::CubicOptions {
                alpha: Some(thermopack_rs::AlphaCorrelation::from_name(config[1])?),
                parameter_reference: config[2],
                volume_translation: config[3] == "1",
                ..Default::default()
            },
        )?
    };
    for operation in &config[4..] {
        let data: Vec<_> = operation.split(':').collect();
        match data[0] {
            "init" => model.init_peneloux_volume_translation(data[1])?,
            "disable" => model.disable_volume_translation(),
            "bounds" => {} // Applied to native caloric options below.
            "qhelper" => {}
            "rr" => {}      // RR iteration inputs, handled after model construction.
            "quality" => {} // Explicit pure two-phase acoustic state.
            "roots" => {}   // Optional fluid acoustic phase-root overrides.
            "amount" => {}  // Total moles in an extensive TV acoustic input.
            "guess" => {}   // UV/SV temperature and pressure guesses.
            "kij" => model.set_kij(
                data[1].parse().unwrap(),
                data[2].parse().unwrap(),
                data[3].parse().unwrap(),
            )?,
            "lij" => model.set_lij(
                data[1].parse().unwrap(),
                data[2].parse().unwrap(),
                data[3].parse().unwrap(),
            )?,
            "beta" => {
                let params = data[3]
                    .split(',')
                    .map(str::parse)
                    .collect::<Result<Vec<f64>, _>>()
                    .map_err(|_| thermopack_rs::Error::InvalidInput("beta parameters"))?;
                model.set_beta(
                    data[1].parse().unwrap(),
                    thermopack_rs::Beta::new(
                        thermopack_rs::BetaCorrelation::from_name(data[2])?,
                        params.try_into().map_err(|_| {
                            thermopack_rs::Error::InvalidInput("two beta parameters")
                        })?,
                    )?,
                )?;
            }
            "set" => {
                let coefficients = data[3]
                    .split(',')
                    .map(str::parse)
                    .collect::<Result<Vec<f64>, _>>()
                    .map_err(|_| thermopack_rs::Error::InvalidInput("translation coefficients"))?;
                model.set_volume_translation(
                    data[1]
                        .parse()
                        .map_err(|_| thermopack_rs::Error::InvalidInput("component index"))?,
                    thermopack_rs::VolumeTranslation::new(
                        thermopack_rs::TranslationKind::from_code(data[2].parse().unwrap())?,
                        coefficients.try_into().map_err(|_| {
                            thermopack_rs::Error::InvalidInput("six coefficients required")
                        })?,
                    )?,
                )?;
            }
            _ => return Err(thermopack_rs::Error::InvalidInput("translation operation")),
        }
    }
    Ok(model)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for line in io::stdin().lock().lines() {
        let line = line?;
        let fields: Vec<_> = line.split('|').collect();
        if fields[0] != "SOLIDCO2" && fields[0] != "HYPERDUAL" && fields.len() != 6 {
            return Err("expected six fields".into());
        }
        if fields[0] == "MEOS" {
            let outcome = (|| -> thermopack_rs::Result<Vec<f64>> {
                use thermopack_rs::{LennardJonesKind, Multiparameter};
                let config: Vec<_> = fields[1].split(':').collect();
                let model = match config[0] {
                    "LJ_MEOS" => Multiparameter::lennard_jones(config[1], LennardJonesKind::Full)?,
                    "LJTS_MEOS" => Multiparameter::lennard_jones(
                        config[1],
                        LennardJonesKind::TruncatedShifted,
                    )?,
                    "LJTS_REDUCED" => Multiparameter::reduced_lennard_jones_truncated_shifted(),
                    _ => Multiparameter::new(config[1])?,
                };
                if fields[5] == "DENSITY" {
                    let t = fields[2].parse().unwrap();
                    let p = fields[3].parse().unwrap();
                    let phase = fields[4].parse().unwrap();
                    let prep = model.prepare(t)?;
                    let r = prep.baseline_density_with_phase_code(p, phase)?;
                    let strict = match phase {
                        1 => Some(Phase::Liquid),
                        2 => Some(Phase::Vapor),
                        3 => Some(Phase::MinimumGibbs),
                        _ => None,
                    }
                    .and_then(|phase| prep.volume(p, phase).ok());
                    return Ok(vec![
                        r.density,
                        r.phase as f64,
                        r.source_error_code as f64,
                        r.pressure,
                        r.pressure_density_derivative,
                        r.iterations as f64,
                        r.phase_switches as f64,
                        f64::from(r.pressure_converged),
                        f64::from(r.mechanically_stable),
                        r.relative_pressure_error,
                        strict.unwrap_or(f64::NAN),
                        f64::from(r.relaxed_convergence),
                    ]);
                }
                if fields[5] == "META" {
                    let r = model.record();
                    let (mt, mp) = model.lennard_jones_limits().unwrap_or((0., 0.));
                    let (sigma, eps) = model
                        .lennard_jones_parameters()
                        .map_or((0., 0.), |p| (p.sigma, p.epsilon_over_k));
                    let (np, ne, ng) = if config[0] == "NIST_MEOS" {
                        (0, 0, 0)
                    } else {
                        (
                            r.pol_terms.len(),
                            r.pol_terms.len() + r.exp_terms.len(),
                            r.pol_terms.len() + r.exp_terms.len() + r.expexp_terms.len(),
                        )
                    };
                    return Ok(vec![
                        r.tc,
                        r.pc,
                        r.rc,
                        r.acf,
                        model.baseline_molar_mass(),
                        r.t_triple,
                        r.p_triple,
                        r.rholiq_triple,
                        r.rhovap_triple,
                        r.rgas_fit,
                        model.gas_constant(),
                        mt,
                        mp,
                        sigma,
                        eps,
                        np as f64,
                        ne as f64,
                        ng as f64,
                    ]);
                }
                let t = fields[2].parse().unwrap();
                let v: f64 = fields[3].parse().unwrap();
                let n: f64 = fields[4].parse().unwrap();
                let prep = model.prepare(t)?;
                let s = prep.state(v)?;
                let i = s.ideal;
                let r = s.residual;
                let mut values = vec![
                    i.alpha,
                    i.d_tau,
                    i.d_tautau,
                    i.d_delta,
                    i.d_deltadelta,
                    i.d_deltatau,
                    r.alpha,
                    r.d_tau,
                    r.d_tautau,
                    r.d_delta,
                    r.d_deltadelta,
                    r.d_deltatau,
                    s.pressure.value,
                    s.pressure.dt,
                    s.pressure.dv,
                    s.internal_energy,
                    s.entropy,
                    s.enthalpy,
                    s.helmholtz_energy,
                    s.gibbs_energy,
                    s.heat_capacity_v,
                    s.heat_capacity_p,
                    s.ln_fugacity.unwrap_or(f64::NAN),
                    s.baseline_speed_of_sound.unwrap_or(f64::NAN),
                ];
                for a in [
                    prep.residual_helmholtz(v * n, n)?,
                    prep.ideal_helmholtz(v * n, n)?,
                ] {
                    values.extend([
                        a.value, a.dt, a.dv, a.dn, a.dtt, a.dtv, a.dtn, a.dvv, a.dvn, a.dnn,
                    ]);
                }
                values.extend([
                    prep.sat_delta_estimate(Phase::Liquid).unwrap_or(f64::NAN),
                    prep.sat_delta_estimate(Phase::Vapor).unwrap_or(f64::NAN),
                ]);
                Ok(values)
            })();
            match outcome {
                Ok(values) => println!(
                    "{}",
                    values
                        .iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
                Err(e) => println!("ERROR {e}"),
            }
            continue;
        }
        if fields[0] == "FMT" {
            let outcome = (|| -> thermopack_rs::Result<Vec<f64>> {
                let model = thermopack_rs::fmt::FmtModel::from_code(fields[1].parse().unwrap())?;
                let n: Vec<f64> = fields[2].split(',').map(|v| v.parse().unwrap()).collect();
                let s =
                    model.baseline_energy_density(n.try_into().map_err(|_| {
                        thermopack_rs::Error::InvalidInput("six weighted densities")
                    })?)?;
                let flags: u8 = fields[3].parse().unwrap();
                let mut values = vec![if flags & 1 != 0 { s.value } else { 0. }];
                values.extend(s.gradient.map(|v| if flags & 2 != 0 { v } else { 0. }));
                values.extend(
                    s.hessian
                        .into_iter()
                        .flatten()
                        .map(|v| if flags & 4 != 0 { v } else { 0. }),
                );
                Ok(values)
            })();
            match outcome {
                Ok(values) => println!(
                    "{}",
                    values
                        .iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
                Err(e) => println!("ERROR {e}"),
            }
            continue;
        }
        if fields[0] == "SOLIDH2O" {
            let outcome = (|| -> thermopack_rs::Result<Vec<f64>> {
                let mut model = thermopack_rs::SolidWater::new();
                if fields[1] != "raw" {
                    for shift in fields[1].split(':') {
                        let a: Vec<f64> = shift.split(',').map(|v| v.parse().unwrap()).collect();
                        model.initialize(a[0], a[1])?;
                    }
                }
                let s = model.state(fields[2].parse().unwrap(), fields[3].parse().unwrap())?;
                let g = s.gibbs;
                Ok(vec![
                    g.value,
                    g.dp,
                    g.dpp,
                    g.dt,
                    g.dtt,
                    g.dtp,
                    s.volume,
                    s.entropy,
                    s.enthalpy,
                    s.internal_energy,
                    s.helmholtz_energy,
                    s.baseline_heat_capacity_p,
                ])
            })();
            match outcome {
                Ok(values) => println!(
                    "{}",
                    values
                        .iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
                Err(e) => println!("ERROR {e}"),
            }
            continue;
        }
        if fields[0] == "SOLIDCO2" {
            let mut model = thermopack_rs::solid_co2::SolidCo2::new();
            model.initialize(
                fields[1].parse().unwrap(),
                fields[2].parse().unwrap(),
                fields[3].parse().unwrap(),
                fields[4].parse().unwrap(),
            );
            let temperature: f64 = fields[5].parse().unwrap();
            let pressure: f64 = fields[6].parse().unwrap();
            let values = [
                model.gibbs(temperature, pressure),
                model.dg_dp(temperature, pressure),
                model.d2g_dp2(temperature, pressure),
                model.dg_dt(temperature, pressure),
                model.d2g_dt2(temperature, pressure),
                model.d2g_dtdp(temperature, pressure),
                model.specific_volume(temperature, pressure),
                model.entropy(temperature, pressure),
                model.enthalpy(temperature, pressure),
                model.internal_energy(temperature, pressure),
                model.helmholtz(temperature, pressure),
                model.heat_capacity(temperature, pressure),
                model.speed_of_sound(temperature, pressure),
            ];
            println!(
                "{}",
                values
                    .iter()
                    .map(|v| format!("{v:.17e}"))
                    .collect::<Vec<_>>()
                    .join(",")
            );
            continue;
        }
        if fields[0] == "SOLIDSAT" {
            let component = fields[1];
            let temperature: f64 = fields[2].parse().unwrap();
            let pressure: f64 = fields[3].parse().unwrap();
            let phase = match fields[4].parse::<i32>().unwrap() {
                1 => Phase::Liquid,
                _ => Phase::Vapor,
            };
            let specification: Vec<&str> = fields[5].split(':').collect();
            let kind = specification[0];
            let value: f64 = specification[1].parse().unwrap();
            let model =
                thermopack_rs::Cubic::new(thermopack_rs::CubicKind::PengRobinson, &[component])?;
            let solid = thermopack_rs::solid_saturation::solid_co2_with_fluid_reference(&model)?;
            let ideal =
                thermopack_rs::IdealGas::new(thermopack_rs::component::component(component)?)?;
            let (specification, initial) = match kind {
                "T" => (
                    thermopack_rs::solid_saturation::SolidSaturationSpecification::FixedTemperature(
                        value,
                    ),
                    (value, pressure),
                ),
                _ => (
                    thermopack_rs::solid_saturation::SolidSaturationSpecification::FixedPressure(
                        value,
                    ),
                    (temperature, value),
                ),
            };
            let solution = thermopack_rs::solid_saturation::solid_co2_saturation(
                &model,
                &solid,
                &ideal,
                specification,
                initial,
                phase,
            )?;
            println!(
                "{:.17e},{:.17e},{:.17e}",
                solution.temperature, solution.pressure, solution.residual
            );
            continue;
        }
        if fields[0] == "HYPERDUAL" {
            let inputs: Vec<f64> = fields[1]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let mut number =
                thermopack_rs::hyperdual::HyperDual::with_order(inputs[0], inputs[8] as u8);
            number.f1 = inputs[1];
            number.f2 = inputs[2];
            number.f3 = inputs[3];
            number.f12 = inputs[4];
            number.f13 = inputs[5];
            number.f23 = inputs[6];
            number.f123 = inputs[7];
            let functions: Vec<thermopack_rs::hyperdual::HyperDual> = vec![
                number.exp(),
                number.ln(),
                number.sqrt(),
                number.sin(),
                number.cos(),
                number.tan(),
                number.sinh(),
                number.cosh(),
                number.tanh(),
                number.asin(),
                number.acos(),
                number.atan(),
            ];
            let mut values = Vec::with_capacity(96);
            for value in functions {
                values.extend([
                    value.f0, value.f1, value.f2, value.f3, value.f12, value.f13, value.f23,
                    value.f123,
                ]);
            }
            println!(
                "{}",
                values
                    .iter()
                    .map(|v| format!("{v:.17e}"))
                    .collect::<Vec<_>>()
                    .join(",")
            );
            continue;
        }
        if fields[0] == "MUT" {
            let components: Vec<&str> = fields[1].split(',').collect();
            let temperature: f64 = fields[2].parse().unwrap();
            let densities: Vec<f64> = fields[4]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let config: Vec<&str> = fields[5].split(':').collect();
            let model =
                thermopack_rs::Cubic::new(thermopack_rs::CubicKind::PengRobinson, &components)?;
            let target: Vec<f64> = if config.len() > 1 {
                config[1]
                    .split(',')
                    .map(|value| value.parse().unwrap())
                    .collect()
            } else {
                Vec::new()
            };
            let kind = config[0];
            let outcome = match kind {
                "MU" => {
                    thermopack_rs::mut_flash::solve_mu_t(&model, temperature, &target, &densities)
                        .map(|solution| solution.density)
                }
                "LNF" => {
                    thermopack_rs::mut_flash::solve_lnf_t(&model, temperature, &target, &densities)
                        .map(|solution| solution.density)
                }
                "MUTARGET" => thermopack_rs::mut_flash::log_fugacity_and_potential(
                    &model,
                    temperature,
                    &densities,
                )
                .map(|(_, potential)| potential),
                "LNFTARGET" => thermopack_rs::mut_flash::log_fugacity_and_potential(
                    &model,
                    temperature,
                    &densities,
                )
                .map(|(fugacity, _)| fugacity),
                _ => Err(thermopack_rs::Error::InvalidInput("unknown mu-T mode")),
            };
            let values = match outcome {
                Ok(values) => values,
                Err(error) => {
                    println!("ERROR: {error}");
                    continue;
                }
            };
            println!(
                "{}",
                values
                    .iter()
                    .map(|v| format!("{v:.17e}"))
                    .collect::<Vec<_>>()
                    .join(",")
            );
            continue;
        }
        if fields[0] == "SOLID" {
            let component = fields[1];
            let kind = match fields[5] {
                "MELT" => 0,
                "SUBL" => 1,
                other => return Err(format!("unknown solid line kind: {other}").into()),
            };
            let points: usize = fields[3].parse().unwrap();
            let outcome = if kind == 0 {
                thermopack_rs::solid::melting_line(component, fields[2].parse().unwrap(), points)
            } else {
                thermopack_rs::solid::sublimation_line(
                    component,
                    fields[2].parse().unwrap(),
                    points,
                )
            };
            let line = match outcome {
                Ok(line) => line,
                Err(error) => {
                    println!("ERROR: {error}");
                    continue;
                }
            };
            let mut values = Vec::with_capacity(2 * points);
            values.extend(line.temperature.iter().copied());
            values.extend(line.pressure.iter().copied());
            println!(
                "{}",
                values
                    .iter()
                    .map(|v| format!("{v:.17e}"))
                    .collect::<Vec<_>>()
                    .join(",")
            );
            continue;
        }
        if fields[0] == "UNIFAC" {
            let components: Vec<&str> = fields[1].split(',').collect();
            let composition: Vec<f64> = fields[4]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let variant = match fields[5] {
                "UNIFAC" => thermopack_rs::UnifacVariant::Unifac,
                "UMR" => thermopack_rs::UnifacVariant::Umr,
                "VTPR" => thermopack_rs::UnifacVariant::Vtpr,
                other => return Err(format!("unknown unifac variant: {other}").into()),
            };
            let model = thermopack_rs::Unifac::new(&components, variant)?;
            let result = model.excess_gibbs(fields[2].parse().unwrap(), &composition)?;
            let mut values = vec![result.gibbs, result.dg_dt, result.d2g_dt2];
            values.extend(result.dg_dn.iter().copied());
            values.extend(result.d2g_dn_dt.iter().copied());
            values.extend(result.d2g_dn_dn.iter().copied());
            println!(
                "{}",
                values
                    .iter()
                    .map(|v| format!("{v:.17e}"))
                    .collect::<Vec<_>>()
                    .join(",")
            );
            continue;
        }
        if fields[0] == "MBWR" {
            let model = thermopack_rs::Mbwr::from_record(fields[1])?;
            let result = (|| -> thermopack_rs::Result<Vec<f64>> {
                if fields[5] == "MBWR_INITIALIZATION" {
                    let s = model.source_initialization()?;
                    let a = s.parameters;
                    return Ok(vec![
                        a.critical_temperature,
                        a.critical_pressure,
                        a.critical_density,
                        s.compressibility,
                        a.acentric_factor,
                        s.srk_covolume,
                        s.srk_alpha_coefficient,
                        s.srk_attraction,
                    ]);
                }
                if fields[5].starts_with("MBWR_CRITICAL_GRID@") {
                    let a: Vec<f64> = fields[5]
                        .split('@')
                        .nth(1)
                        .unwrap()
                        .split(',')
                        .map(|s| s.parse().unwrap())
                        .collect();
                    let s = model.baseline_critical_grid(
                        thermopack_rs::MbwrDensityParameters::new(a[0], a[1], a[2], 0.)?,
                    )?;
                    return Ok(vec![s.temperature, s.pressure, 1. / s.volume]);
                }
                if fields[5].starts_with("MBWR_DENSITY@") {
                    let a: Vec<f64> = fields[5]
                        .split('@')
                        .nth(1)
                        .unwrap()
                        .split(',')
                        .map(|s| s.parse().unwrap())
                        .collect();
                    let parameters =
                        thermopack_rs::MbwrDensityParameters::new(a[0], a[1], a[2], a[3])?;
                    let phase = if fields[4] == "L" {
                        Phase::Liquid
                    } else {
                        Phase::Vapor
                    };
                    let result = model.baseline_density(
                        fields[2].parse().unwrap(),
                        fields[3].parse().unwrap(),
                        phase,
                        parameters,
                    )?;
                    return Ok(vec![
                        result.density,
                        if result.phase == Phase::Liquid {
                            1.
                        } else {
                            2.
                        },
                        result.pressure,
                        result.pressure_density_derivative,
                        result.path as u8 as f64,
                    ]);
                }
                if fields[5] == "MBWR_METADATA" {
                    let mut values = vec![0.; 8];
                    if let Some(c) = model.source_critical_point() {
                        values[..5].copy_from_slice(&[
                            1.,
                            c.temperature,
                            c.pressure,
                            1. / c.volume,
                            c.compressibility,
                        ]);
                    }
                    if let Some((t, p)) = model.source_triple_point() {
                        values[5..].copy_from_slice(&[1., t, p]);
                    }
                    return Ok(values);
                }
                if fields[5] == "MBWR_CRITICAL" {
                    let c = model.critical_point_from(
                        fields[2].parse().unwrap(),
                        fields[3].parse().unwrap(),
                    )?;
                    return Ok(vec![
                        c.temperature,
                        c.pressure,
                        1. / c.volume,
                        c.compressibility,
                    ]);
                }
                if fields[5] == "MBWR_ROOTS" {
                    let t = fields[2].parse().unwrap();
                    let p = fields[3].parse().unwrap();
                    let maximum = fields[4].parse().unwrap();
                    return model.prepare_equilibrium(t, maximum)?.densities(p);
                }
                if fields[5].starts_with("MBWR_TP_PROPERTIES@") {
                    let p = fields[5].split('@').nth(1).unwrap().parse().unwrap();
                    let t = fields[2].parse().unwrap();
                    let v = fields[3].parse().unwrap();
                    let n = fields[4].parse().unwrap();
                    let state = model.baseline_residual_tp(t, p, v, n)?;
                    let mut values = vec![state.pressure, state.gibbs_energy];
                    for q in [
                        state.enthalpy,
                        state.entropy,
                        state.compressibility,
                        state.ln_fugacity_coefficient,
                    ] {
                        let d = q.derivatives.ok_or(thermopack_rs::Error::NotConverged(
                            "singular residual TP response",
                        ))?;
                        values.extend([q.value, d.dt, d.dp, d.dn]);
                    }
                    return Ok(values);
                }
                if fields[5] == "MBWR_INTEGRATED" {
                    return Ok(model
                        .baseline_helmholtz_terms()
                        .iter()
                        .map(|t| t.coefficient)
                        .collect());
                }
                if fields[5] == "MBWR_PARAMETERS" {
                    let (c, gamma) = model.parameters();
                    let mut values = vec![model.kind().coefficient_count() as f64];
                    values.extend(c);
                    values.push(gamma);
                    return Ok(values);
                }
                let t = fields[2].parse().unwrap();
                let v = fields[3].parse().unwrap();
                let n = fields[4].parse().unwrap();
                let prepared = model.prepare(t)?;
                let baseline = fields[5] == "MBWR_BASELINE_CORE";
                let a = if baseline {
                    prepared.baseline_residual_density(n / v)?
                } else {
                    prepared.residual_density(n / v)?
                };
                let f = if baseline {
                    prepared.baseline_residual_helmholtz(v, n)?
                } else {
                    prepared.residual_helmholtz(v, n)?
                };
                let p = if baseline {
                    prepared.baseline_pressure_tv(v, n)?
                } else {
                    prepared.pressure_tv(v, n)?
                };
                let d = prepared.pressure_density(n / v)?;
                Ok(vec![
                    a.value, a.dt, a.drho, a.dtt, a.dtrho, a.drhorho, f.value, f.dt, f.dv, f.dn,
                    f.dtt, f.dtv, f.dtn, f.dvv, f.dvn, f.dnn, p.value, p.dt, p.dv, p.dn, d.value,
                    d.dt, d.drho, d.dtt, d.dtrho, d.drhorho,
                ])
            })();
            match result {
                Ok(v) => println!(
                    "{}",
                    v.iter()
                        .map(|x| format!("{x:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
                Err(e) => println!("ERROR {e}"),
            }
            continue;
        }
        if fields[0] == "PETS" {
            let parameters: Vec<f64> = fields[1]
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?;
            let model = thermopack_rs::Pets::with_parameters(parameters[0], parameters[1])?;
            let result = (|| -> thermopack_rs::Result<Vec<f64>> {
                let config: Vec<_> = fields[5].split('@').collect();
                if matches!(config[0], "PETS_UV" | "PETS_SV") {
                    let target = fields[2].parse().unwrap();
                    let volume = fields[3].parse().unwrap();
                    let options = thermopack_rs::VolumeFlashOptions {
                        temperature_guess: if config[1] == "none" {
                            None
                        } else {
                            Some(config[1].parse().unwrap())
                        },
                        pressure_guess: if config[2] == "none" {
                            None
                        } else {
                            Some(config[2].parse().unwrap())
                        },
                        minimum_temperature: 1.,
                        maximum_temperature: 2000.,
                        minimum_pressure: 1e-12,
                        maximum_pressure: 1e10,
                        ..Default::default()
                    };
                    let flash = if config[0] == "PETS_UV" {
                        model.uv_flash_with_options(target, volume, options)?
                    } else {
                        model.sv_flash_with_options(target, volume, options)?
                    };
                    let state = model.equilibrium_tv(flash.temperature, volume)?;
                    let two = flash.phase == thermopack_rs::FlashPhase::TwoPhase
                        && flash.vapor_fraction > 1e-8
                        && flash.vapor_fraction < 1. - 1e-8;
                    let props = state.properties;
                    return Ok(vec![
                        flash.temperature,
                        flash.pressure,
                        f64::from(two),
                        if two { flash.vapor_fraction } else { 0. },
                        props.volume,
                        props.internal_energy,
                        props.entropy,
                        props.enthalpy,
                    ]);
                }
                if matches!(config[0], "PETS_TP" | "PETS_PH" | "PETS_PS") {
                    let t = fields[2].parse().unwrap();
                    let p = fields[3].parse().unwrap();
                    let options = thermopack_rs::CaloricFlashOptions {
                        temperature_guess: if config[2] == "none" {
                            None
                        } else {
                            Some(config[2].parse().unwrap())
                        },
                        minimum_temperature: 1.,
                        maximum_temperature: 2000.,
                        ..Default::default()
                    };
                    let flash = match config[0] {
                        "PETS_TP" => model.tp_flash(t, p)?,
                        "PETS_PH" => {
                            model.ph_flash_with_options(p, config[1].parse().unwrap(), options)?
                        }
                        _ => model.ps_flash_with_options(p, config[1].parse().unwrap(), options)?,
                    };
                    let state = model.equilibrium_properties(&flash)?;
                    let two = flash.phase == thermopack_rs::FlashPhase::TwoPhase
                        && flash.vapor_fraction > 1e-8
                        && flash.vapor_fraction < 1. - 1e-8;
                    return Ok(vec![
                        flash.temperature,
                        p,
                        f64::from(two),
                        if two { flash.vapor_fraction } else { 0. },
                        if two {
                            flash.liquid_volume.unwrap()
                        } else {
                            state.volume
                        },
                        if two {
                            flash.vapor_volume.unwrap()
                        } else {
                            state.volume
                        },
                        state.enthalpy,
                        state.entropy,
                    ]);
                }
                if config[0] == "PETS_CRITICAL" {
                    let c = model.critical_point()?;
                    return Ok(vec![c.temperature, c.pressure, c.volume, c.compressibility]);
                }
                if config[0] == "PETS_VOLUME_BASELINE" || config[0] == "PETS_VOLUME" {
                    let t = fields[2].parse().unwrap();
                    let p = fields[3].parse().unwrap();
                    let phase = if config[1] == "L" {
                        Phase::Liquid
                    } else {
                        Phase::Vapor
                    };
                    let v = if config[0] == "PETS_VOLUME_BASELINE" {
                        model.baseline_volume_iteration(t, p, phase)?.volume
                    } else {
                        model.volume(t, p, phase)?
                    };
                    let state = model.pressure_tv(t, v, 1.)?;
                    return Ok(if config[0] == "PETS_VOLUME_BASELINE" {
                        vec![v, state.value, state.dt, state.dv]
                    } else {
                        vec![v, -state.dt / state.dv, 1. / state.dv]
                    });
                }
                if config[0] == "PETS_SAT" || config[0] == "PETS_SAT_T" {
                    let sat = if config[0] == "PETS_SAT" {
                        model.saturation_pressure(fields[2].parse().unwrap())?
                    } else {
                        model.saturation_temperature(fields[3].parse().unwrap())?
                    };
                    return Ok(vec![
                        sat.temperature,
                        sat.pressure,
                        sat.liquid_volume,
                        sat.vapor_volume,
                        sat.dpdt,
                    ]);
                }
                if fields[5] == "PETS_POTENTIAL" {
                    let r: Vec<f64> = fields[2].split(',').map(|v| v.parse().unwrap()).collect();
                    let mut values = vec![0.; r.len()];
                    let mut forces = vec![0.; r.len()];
                    model.potential_into(&r, &mut values, Some(&mut forces))?;
                    values.extend(forces);
                    return Ok(values);
                }
                let t = fields[2].parse().unwrap();
                let v = fields[3].parse().unwrap();
                let n = fields[4].parse().unwrap();
                let prepared = model.prepare(t)?;
                let diameter = prepared.diameter();
                let a = prepared.contributions(n / v)?;
                let f = prepared.residual_helmholtz(v, n)?;
                let c = model.lj_truncation_correction(t, v, n)?;
                let p = prepared.caloric_tv(v, n)?;
                let mut values = vec![diameter.value, diameter.dt, diameter.dtt];
                for a in [a.total, a.hard_sphere, a.dispersion] {
                    values.extend([a.value, a.dt, a.drho, a.dtt, a.dtrho, a.drhorho]);
                }
                for f in [f, c] {
                    values.extend([
                        f.value, f.dt, f.dv, f.dn, f.dtt, f.dtv, f.dtn, f.dvv, f.dvn, f.dnn,
                    ]);
                }
                values.extend([
                    p.pressure.value,
                    p.pressure.dt,
                    p.pressure.dv,
                    p.pressure.dn,
                    p.internal_energy,
                    p.entropy,
                    p.helmholtz_energy,
                    p.enthalpy,
                    p.heat_capacity_v,
                    p.chemical_potential,
                ]);
                let (b, c) = model.baseline_virial_coefficients(t)?;
                values.extend([
                    prepared.baseline_speed_of_sound_tv(v, n)?,
                    b,
                    c,
                    p.ln_fugacity_coefficient.unwrap_or(0.),
                    p.heat_capacity_p.unwrap_or(f64::NAN),
                ]);
                Ok(values)
            })();
            match result {
                Ok(values) => println!(
                    "{}",
                    values
                        .iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if fields[5] == "SYMMETRIC" {
            let mut a: Vec<f64> = fields[1]
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?;
            let mut b: Vec<f64> = fields[2]
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?;
            let mut pivots = vec![0; b.len()];
            let mut work = vec![0.; b.len()];
            match thermopack_rs::numerics::solve_symmetric_in_place(
                &mut a,
                &mut b,
                &mut pivots,
                &mut work,
            ) {
                Ok(success) => {
                    print!("{}", usize::from(success));
                    for value in pivots {
                        print!(",{value}");
                    }
                    for value in a.iter().chain(&b) {
                        print!(",{value:.17e}");
                    }
                    println!();
                }
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if fields[5] == "RR" {
            let z: Vec<f64> = fields[1]
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?;
            let k: Vec<f64> = fields[2]
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?;
            let mut state = thermopack_rs::numerics::RachfordRiceState::new(z.len());
            state.vapor_fraction = fields[3].parse()?;
            match thermopack_rs::numerics::rachford_rice_extended(
                &z,
                &k,
                &mut state,
                fields[4] == "1",
            ) {
                Ok(false) => println!("0"),
                Ok(true) => {
                    print!(
                        "1,{:.17e},{:.17e}",
                        state.vapor_fraction, state.liquid_fraction
                    );
                    for value in state.liquid.iter().chain(&state.vapor) {
                        print!(",{value:.17e}");
                    }
                    println!();
                }
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if fields[5] == "IDEAL" || fields[5] == "IDEAL_CONFIG" {
            let t: f64 = fields[2].parse()?;
            let p: f64 = fields[3].parse()?;
            let mut gas = thermopack_rs::IdealGas::new(thermopack_rs::component(fields[1])?)?;
            if fields[5] == "IDEAL_CONFIG" {
                for action in fields[4].split('/') {
                    let tokens: Vec<_> = action.split(':').collect();
                    match tokens[0] {
                        "cp" => gas.set_heat_capacity(
                            tokens[1].parse()?,
                            &tokens[2]
                                .split(',')
                                .filter(|v| !v.is_empty())
                                .map(str::parse)
                                .collect::<Result<Vec<_>, _>>()?,
                        )?,
                        "h" => gas.set_formation_enthalpy(tokens[1].parse()?)?,
                        "s" => gas.set_standard_entropy(
                            tokens[1].parse()?,
                            if tokens[2] == "101325.0" {
                                thermopack_rs::EntropyReference::Atmosphere
                            } else {
                                thermopack_rs::EntropyReference::Bar
                            },
                        )?,
                        _ => return Err("unknown ideal action".into()),
                    }
                }
            }
            match gas.properties(t, p) {
                Ok(properties) => {
                    print!(
                        "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                        properties.enthalpy,
                        properties.entropy,
                        properties.cp,
                        properties.cp / t,
                        -thermopack_rs::R / p
                    );
                    if fields[5] == "IDEAL_CONFIG" {
                        let (corr, params) = gas.heat_capacity_correlation();
                        print!(
                            ",{:.17e},{:.17e},{corr}",
                            gas.formation_enthalpy(),
                            gas.standard_entropy()
                        );
                        for value in params {
                            print!(",{value:.17e}");
                        }
                    }
                    println!();
                }
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        let kind = match fields[0] {
            "PR" => CubicKind::PengRobinson,
            "SRK" => CubicKind::SoaveRedlichKwong,
            "RK" => CubicKind::RedlichKwong,
            "VDW" => CubicKind::VanDerWaals,
            _ => return Err("unknown model".into()),
        };
        if fields[5] == "ALPHA" {
            let result = (|| -> Result<Vec<f64>, thermopack_rs::Error> {
                let selection: Vec<_> = fields[4].split('@').collect();
                let mut model = Cubic::with_alpha(
                    kind,
                    &fields[1].split(',').collect::<Vec<_>>(),
                    thermopack_rs::AlphaCorrelation::from_name(fields[3])?,
                    selection[0],
                )?;
                for setter in &selection[1..] {
                    let data: Vec<_> = setter.split(':').collect();
                    let parameters = data[2]
                        .split(',')
                        .map(str::parse)
                        .collect::<Result<Vec<f64>, _>>()
                        .map_err(|_| {
                            thermopack_rs::Error::InvalidInput("invalid alpha parameter")
                        })?;
                    model.set_alpha(
                        data[0].parse().map_err(|_| {
                            thermopack_rs::Error::InvalidInput("invalid component index")
                        })?,
                        thermopack_rs::AlphaCorrelation::from_name(data[1])?,
                        &parameters,
                    )?;
                }
                let mut values = Vec::new();
                for (i, component) in model.components().iter().enumerate() {
                    let alpha = model.alpha(i)?;
                    let value = alpha
                        .evaluate(fields[2].parse().unwrap(), component.critical_temperature)?;
                    values.extend([alpha.correlation as i32 as f64]);
                    values.extend(alpha.parameters);
                    values.extend([value.value, value.dt, value.dtt]);
                }
                Ok(values)
            })();
            match result {
                Ok(values) => println!(
                    "{}",
                    values
                        .iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        let property_config: Vec<_> = fields[5].split('@').collect();
        let alpha_properties = property_config[0].starts_with("A_");
        let translated = property_config[0].starts_with("V_");
        // The UNIFAC mixing rule is selected with a "unifac" entry in the
        // property configuration.
        let unifac = property_config
            .iter()
            .any(|entry| entry.eq_ignore_ascii_case("unifac"));
        let mut model = if translated {
            match translated_model(kind, fields[1], &property_config) {
                Ok(model) => model,
                Err(error) => {
                    println!("ERROR: {error}");
                    continue;
                }
            }
        } else if alpha_properties {
            Cubic::with_alpha(
                kind,
                &fields[1].split(',').collect::<Vec<_>>(),
                thermopack_rs::AlphaCorrelation::from_name(property_config[1])?,
                property_config[2],
            )?
        } else if unifac {
            match Cubic::with_options(
                kind,
                &fields[1].split(',').collect::<Vec<_>>(),
                thermopack_rs::CubicOptions {
                    mixing_rule: thermopack_rs::MixingRule::Unifac,
                    ..Default::default()
                },
            ) {
                Ok(model) => model,
                Err(error) => {
                    println!("ERROR: {error}");
                    continue;
                }
            }
        } else {
            Cubic::new(kind, &fields[1].split(',').collect::<Vec<_>>())?
        };
        if alpha_properties {
            for setter in &property_config[3..] {
                let data: Vec<_> = setter.split(':').collect();
                let parameters = data[2]
                    .split(',')
                    .map(str::parse)
                    .collect::<Result<Vec<f64>, _>>()?;
                model.set_alpha(
                    data[0].parse()?,
                    thermopack_rs::AlphaCorrelation::from_name(data[1])?,
                    &parameters,
                )?;
            }
        }
        let t: f64 = fields[2].parse()?;
        if property_config[0] == "V_BETA" {
            let mut values = Vec::new();
            for (i, component) in model.components().iter().enumerate() {
                let beta = model.beta(i)?;
                let value = beta.evaluate(t, component.critical_temperature)?;
                values.push(beta.correlation as i32 as f64);
                values.extend(
                    if beta.correlation == thermopack_rs::BetaCorrelation::Classic {
                        [0.; 2]
                    } else {
                        beta.parameters
                    },
                );
                values.extend([
                    value.value,
                    value.dt,
                    value.dtt,
                    component.critical_temperature,
                    component.critical_pressure,
                    component.molar_mass,
                ]);
            }
            println!(
                "{}",
                values
                    .iter()
                    .map(|v| format!("{v:.17e}"))
                    .collect::<Vec<_>>()
                    .join(",")
            );
            continue;
        }
        if property_config[0] == "V_CI" {
            let mut values = Vec::new();
            for i in 0..model.components().len() {
                let ci = model.volume_translation(i)?;
                let c = ci.evaluate(t)?;
                values.extend([
                    if model.volume_translation_enabled() {
                        1.
                    } else {
                        0.
                    },
                    ci.kind as u8 as f64,
                ]);
                values.extend(ci.coefficients);
                values.extend([
                    c.value,
                    c.dt,
                    c.dtt,
                    if ci.kind == thermopack_rs::TranslationKind::Constant {
                        0.
                    } else {
                        1.
                    },
                ]);
            }
            println!(
                "{}",
                values
                    .iter()
                    .map(|v| format!("{v:.17e}"))
                    .collect::<Vec<_>>()
                    .join(",")
            );
            continue;
        }
        if fields[3] == "PSEUDO" {
            let x: Vec<f64> = fields[4]
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?;
            let critical = model.pseudo_critical(&x)?;
            println!(
                "{:.17e},{:.17e},{:.17e},{:.17e}",
                critical.temperature, critical.pressure, critical.volume, critical.compressibility
            );
            continue;
        }
        if property_config[0] == "V_SAT" {
            match model.saturation_pressure(t) {
                Ok(sat) => println!(
                    "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                    sat.temperature, sat.pressure, sat.liquid_volume, sat.vapor_volume, sat.dpdt
                ),
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if fields[3] == "SAT_T" {
            let sat = model.saturation_temperature(t)?;
            println!(
                "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                sat.temperature, sat.pressure, sat.liquid_volume, sat.vapor_volume, sat.dpdt
            );
            continue;
        }
        if fields[3] == "SAT" {
            let sat = model.saturation_pressure(t)?;
            println!(
                "{:.17e},{:.17e},{:.17e},{:.17e}",
                sat.pressure, sat.liquid_volume, sat.vapor_volume, sat.dpdt
            );
            continue;
        }
        let p: f64 = fields[3].parse()?;
        let x: Vec<f64> = fields[4]
            .split(',')
            .map(str::parse)
            .collect::<Result<_, _>>()?;
        if property_config[0] == "V_CALORIC_TV" {
            match model.prepare(t, &x).and_then(|state| state.caloric_tv(p)) {
                Ok(value) => println!(
                    "{:.17e},{:.17e},{:.17e}",
                    value.pressure, value.internal_energy, value.entropy
                ),
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if matches!(property_config[0], "V_UV" | "V_SV") {
            let mut options = model.volume_flash_options();
            if let Some(guess) = property_config.iter().find(|op| op.starts_with("guess:")) {
                let values: Vec<_> = guess.split(':').collect();
                options.temperature_guess = Some(values[1].parse()?);
                options.pressure_guess = Some(values[2].parse()?);
            }
            let result = if property_config[0] == "V_UV" {
                model.uv_flash_with_options(&x, t, p, options)
            } else {
                model.sv_flash_with_options(&x, t, p, options)
            };
            match result {
                Ok(flash) => {
                    let two = flash.phase == thermopack_rs::FlashPhase::TwoPhase;
                    let mut values = vec![
                        flash.temperature,
                        flash.pressure,
                        f64::from(two),
                        if two { flash.vapor_fraction } else { 0. },
                        if two { flash.liquid_fraction } else { 1. },
                    ];
                    values.extend(if two { flash.liquid } else { x.clone() });
                    values.extend(if two { flash.vapor } else { vec![0.; x.len()] });
                    println!(
                        "{}",
                        values
                            .iter()
                            .map(|v| format!("{v:.17e}"))
                            .collect::<Vec<_>>()
                            .join(",")
                    );
                }
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if property_config[0] == "V_SOUND_TV" {
            let amount: f64 = property_config
                .iter()
                .find(|op| op.starts_with("amount:"))
                .map(|op| op.split(':').nth(1).unwrap().parse().unwrap())
                .unwrap_or(1.);
            let moles: Vec<_> = x.iter().map(|z| z * amount).collect();
            match model.baseline_speed_of_sound_tv(t, p * amount, &moles) {
                Ok(value) => println!("{value:.17e}"),
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if property_config[0] == "V_SOUND" {
            let flash = if let Some(quality) =
                property_config.iter().find(|op| op.starts_with("quality:"))
            {
                let bv: f64 = quality.split(':').nth(1).unwrap().parse()?;
                Ok(thermopack_rs::FlashResult {
                    temperature: t,
                    pressure: p,
                    feed: x.clone(),
                    phase: thermopack_rs::FlashPhase::TwoPhase,
                    liquid: x.clone(),
                    vapor: x.clone(),
                    vapor_fraction: bv,
                    liquid_fraction: 1. - bv,
                    liquid_volume: None,
                    vapor_volume: None,
                    iterations: 0,
                    fugacity_residual: 0.,
                })
            } else {
                model.baseline_tp_flash(t, p, &x)
            };
            let mut roots = [Phase::Liquid, Phase::Vapor];
            if let Some(data) = property_config.iter().find(|op| op.starts_with("roots:")) {
                for (root, label) in roots.iter_mut().zip(data.split(':').skip(1)) {
                    *root = if label == "L" {
                        Phase::Liquid
                    } else {
                        Phase::Vapor
                    };
                }
            }
            match flash.and_then(|f| {
                Ok((
                    f.phase,
                    model.baseline_equilibrium_speed_of_sound_with_phases(&f, roots)?,
                ))
            }) {
                Ok((phase, value)) => println!("{},{value:.17e}", phase as i32),
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if property_config[0] == "V_RESP" {
            match model.baseline_tp_flash(t, p, &x).and_then(|f| {
                let r = model.baseline_equilibrium_response(&f)?.ok_or(
                    thermopack_rs::Error::NotConverged("undefined equilibrium TP response"),
                )?;
                Ok((f.phase, r))
            }) {
                Ok((phase, r)) => println!(
                    "{},{:.17e},{:.17e},{:.17e}",
                    phase as i32, r.heat_capacity, r.volume_dt, r.volume_dp
                ),
                Err(e) => println!("ERROR: {e}"),
            }
            continue;
        }
        if property_config[0] == "V_RR" {
            let inputs: Vec<_> = property_config
                .iter()
                .find(|v| v.starts_with("rr:"))
                .ok_or("missing RR inputs")?
                .split(':')
                .collect();
            let k = inputs[1]
                .split(',')
                .map(str::parse)
                .collect::<Result<Vec<f64>, _>>()?;
            let mut state = thermopack_rs::RachfordRiceIteration::new(x.len());
            state.split.vapor_fraction = inputs[2].parse()?;
            match model.rr_successive_substitution_iteration(
                t,
                p,
                &x,
                &k,
                &mut state,
                inputs[3] == "1",
            ) {
                Ok(()) => {
                    let mut values =
                        vec![f64::from(state.converged), f64::from(state.has_solution)];
                    if state.has_solution {
                        values.extend([
                            state.split.vapor_fraction,
                            state.split.liquid_fraction,
                            state.gibbs,
                        ]);
                        values.extend(
                            state
                                .split
                                .liquid
                                .iter()
                                .chain(&state.split.vapor)
                                .chain(&state.equilibrium_ratios)
                                .chain(&state.liquid_ln_fugacity)
                                .chain(&state.vapor_ln_fugacity),
                        );
                    }
                    println!(
                        "{}",
                        values
                            .iter()
                            .map(|v| format!("{v:.17e}"))
                            .collect::<Vec<_>>()
                            .join(",")
                    );
                }
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if matches!(property_config[0], "PH" | "PS" | "V_PH" | "V_PS") {
            let mut options = model.caloric_options();
            for operation in &property_config[1..] {
                if let Some(bounds) = operation.strip_prefix("bounds:") {
                    let data: Vec<_> = bounds.split(':').collect();
                    options.minimum_temperature = data[0].parse()?;
                    options.maximum_temperature = data[1].parse()?;
                }
            }
            let result = if property_config[0].ends_with("PH") {
                model.ph_flash_with_options(p, &x, t, options)
            } else {
                model.ps_flash_with_options(p, &x, t, options)
            };
            match result {
                Ok(flash) => {
                    print!(
                        "{:.17e},{},{:.17e},{:.17e}",
                        flash.temperature,
                        flash.phase as i32,
                        flash.vapor_fraction,
                        flash.liquid_fraction
                    );
                    for value in flash.liquid.iter().chain(&flash.vapor) {
                        print!(",{value:.17e}");
                    }
                    println!();
                }
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if matches!(
            property_config[0],
            "TP" | "V_TP" | "TP_BASELINE" | "V_TP_BASELINE"
        ) {
            let result = if property_config[0].ends_with("BASELINE") {
                model.baseline_tp_flash(t, p, &x)
            } else {
                model.tp_flash(t, p, &x)
            };
            match result {
                Ok(flash) => {
                    print!(
                        "{},{:.17e},{:.17e}",
                        flash.phase as i32, flash.vapor_fraction, flash.liquid_fraction
                    );
                    for value in flash.liquid.iter().chain(&flash.vapor) {
                        print!(",{value:.17e}");
                    }
                    println!();
                }
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if matches!(property_config[0], "EQ" | "V_EQ") {
            match model
                .tp_flash(t, p, &x)
                .and_then(|f| Ok((f.phase, model.equilibrium_properties(&f)?)))
            {
                Ok((phase, prop)) => println!(
                    "{},{:.17e},{:.17e},{:.17e},{:.17e}",
                    phase as i32,
                    prop.enthalpy,
                    prop.entropy,
                    prop.heat_capacity.unwrap_or(f64::NAN),
                    prop.volume
                ),
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        let state = model.prepare(t, &x)?;
        let v = state.volume(
            p,
            if property_config[0].ends_with('L') {
                Phase::Liquid
            } else {
                Phase::Vapor
            },
        )?;
        if fields[5].starts_with('D') {
            let (mut dt, mut dp, mut dn) = (
                vec![0.0; x.len()],
                vec![0.0; x.len()],
                vec![0.0; x.len() * x.len()],
            );
            state.fugacity_derivatives_into(v, &mut dt, &mut dp, &mut dn)?;
            let values: Vec<_> = dt.into_iter().chain(dp).chain(dn).collect();
            println!(
                "{}",
                values
                    .iter()
                    .map(|v| format!("{v:.17e}"))
                    .collect::<Vec<_>>()
                    .join(",")
            );
            continue;
        }
        let pressure = state.pressure(v)?;
        let residual = state.residual(v)?;
        let mut phi = vec![0.0; x.len()];
        state.ln_fugacity_into(v, &mut phi)?;
        print!(
            "{v:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            pressure.value,
            pressure.dt,
            pressure.dv,
            residual.enthalpy,
            residual.entropy,
            state.enthalpy(v)?,
            state.entropy(v)?,
            state.baseline_enthalpy_derivative_tp(v)?
        );
        for value in phi {
            print!(",{value:.17e}");
        }
        if alpha_properties || translated {
            let (mut dt, mut dp, mut dn) = (
                vec![0.; x.len()],
                vec![0.; x.len()],
                vec![0.; x.len() * x.len()],
            );
            state.fugacity_derivatives_into(v, &mut dt, &mut dp, &mut dn)?;
            for value in dt.iter().chain(&dp).chain(&dn) {
                print!(",{value:.17e}");
            }
        }
        println!();
    }
    Ok(())
}
