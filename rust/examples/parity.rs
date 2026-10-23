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
        if fields[0] == "PURE_CURVE" {
            use thermopack_rs::pure_saturation_curve::{
                PureSaturationOptions, PureSaturationStart,
            };
            let kind = match fields[1] {
                "SRK" => CubicKind::SoaveRedlichKwong,
                "RK" => CubicKind::RedlichKwong,
                "VDW" => CubicKind::VanDerWaals,
                _ => CubicKind::PengRobinson,
            };
            let names: Vec<_> = fields[2].split(',').collect();
            let model = if fields[3] == "quantum" {
                Cubic::quantum(&names)?
            } else {
                Cubic::with_options(
                    kind,
                    &names,
                    thermopack_rs::CubicOptions {
                        volume_translation: fields[3] == "translated",
                        parameter_reference: fields[5],
                        mixing_rule: match fields[3] {
                            "unifac" => thermopack_rs::MixingRule::Unifac,
                            "hv1" => thermopack_rs::MixingRule::HuronVidal1,
                            "hv2" => thermopack_rs::MixingRule::HuronVidal2,
                            "nrtl" => thermopack_rs::MixingRule::Nrtl,
                            _ => thermopack_rs::MixingRule::VdW,
                        },
                        ..Default::default()
                    },
                )?
            };
            let opts: Vec<_> = fields[4].split(':').collect();
            let options = PureSaturationOptions {
                start: if opts[1] == "T" {
                    PureSaturationStart::Temperature(opts[2].parse()?)
                } else {
                    PureSaturationStart::Pressure(opts[2].parse()?)
                },
                maximum_pressure_step: opts[3].parse()?,
                maximum_points: opts[4].parse()?,
                logarithmic: opts[5] == "1",
            };
            match model.baseline_pure_saturation_curve(opts[0].parse()?, options) {
                Ok(curve) => println!(
                    "{:?}|{}",
                    curve.end,
                    curve
                        .points
                        .iter()
                        .map(|p| [
                            p.temperature,
                            p.pressure,
                            p.liquid_volume,
                            p.vapor_volume,
                            p.dpdt
                        ]
                        .iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(","))
                        .collect::<Vec<_>>()
                        .join(";")
                ),
                Err(error) => println!("ERROR:{error}"),
            }
            continue;
        }
        if fields[0] != "SOLIDCO2"
            && fields[0] != "HYPERDUAL"
            && fields[0] != "PSEUDOC"
            && fields[0] != "ISOLINE"
            && fields[0] != "SPINODAL"
            && fields[0] != "MIXSPIN"
            && fields[0] != "JTI"
            && fields[0] != "VIRIAL"
            && fields[0] != "VIRIALMATRIX"
            && fields[0] != "THIRDVIRIALMATRIX"
            && fields[0] != "LJS"
            && fields[0] != "LJSSPLINE"
            && fields[0] != "LJSHELMHOLTZ"
            && fields[0] != "LJSSTATE"
            && fields[0] != "LJSFUGACITY"
            && fields[0] != "LJSCALORIC"
            && fields[0] != "LJSVOLUME"
            && fields[0] != "LJSPOT"
            && fields[0] != "LJSWCATERM"
            && fields[0] != "WHITEBEAR"
            && fields[0] != "UVPHI"
            && fields[0] != "UVTERMS"
            && fields[0] != "UVVIRIALS"
            && fields[0] != "WCATERM"
            && fields[0] != "UVVIRIALDIFF"
            && fields[0] != "CHAINRULE"
            && fields[0] != "AISLJS"
            && fields[0] != "CPARDF"
            && fields[0] != "PCSAFTRDF"
            && fields[0] != "PCSAFTHS"
            && fields[0] != "PCSAFTCHAIN"
            && fields[0] != "PCSAFTSTATE"
            && fields[0] != "WSMIX"
            && fields[0] != "WSSTATE"
            && fields[0] != "QUAD"
            && fields[0] != "CRITICAL"
            && fields.len() != 6
        {
            return Err("expected six fields".into());
        }
        if fields[0] == "GEDER" {
            let kind = if fields[1] == "SRK" {
                CubicKind::SoaveRedlichKwong
            } else {
                CubicKind::PengRobinson
            };
            let names: Vec<_> = fields[2].split(',').collect();
            let z: Vec<f64> = fields[4].split(',').map(|s| s.parse().unwrap()).collect();
            let spec: Vec<_> = fields[5].split(':').collect();
            let t = spec[0].parse()?;
            let p = spec[1].parse()?;
            let rule = match spec[3] {
                "HV1" => thermopack_rs::MixingRule::HuronVidal1,
                "HV2" => thermopack_rs::MixingRule::HuronVidal2,
                "NRTL" => thermopack_rs::MixingRule::Nrtl,
                _ => thermopack_rs::MixingRule::Unifac,
            };
            let model = Cubic::with_options(
                kind,
                &names,
                thermopack_rs::CubicOptions {
                    mixing_rule: rule,
                    parameter_reference: spec[4],
                    volume_translation: fields[3] == "translated",
                    ..Default::default()
                },
            )?;
            let result = (|| -> thermopack_rs::Result<Vec<f64>> {
                let prepared = model.prepare(t, &z)?;
                let v = prepared.volume(
                    p,
                    if spec[2] == "L" {
                        Phase::Liquid
                    } else {
                        Phase::Vapor
                    },
                )?;
                let n = z.len();
                let mut dt = vec![0.; n];
                let mut dp = vec![0.; n];
                let mut dn = vec![0.; n * n];
                prepared.fugacity_derivatives_into(v, &mut dt, &mut dp, &mut dn)?;
                Ok(dt.into_iter().chain(dp).chain(dn).collect())
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
                Err(e) => println!("ERROR: {e}"),
            }
            continue;
        }
        if fields[0] == "ENVELOPE" {
            use thermopack_rs::envelope::{EnvelopeConstraint, EnvelopeOptions};
            use thermopack_rs::{MixtureSaturation, SaturationKind};
            let kind = match fields[1] {
                "SRK" => CubicKind::SoaveRedlichKwong,
                "RK" => CubicKind::RedlichKwong,
                "VDW" => CubicKind::VanDerWaals,
                _ => CubicKind::PengRobinson,
            };
            let names: Vec<_> = fields[2].split(',').collect();
            let model = if fields[3] == "quantum" {
                Cubic::quantum(&names)?
            } else {
                Cubic::with_options(
                    kind,
                    &names,
                    thermopack_rs::CubicOptions {
                        volume_translation: fields[3] == "translated",
                        mixing_rule: if fields[3] == "unifac" {
                            thermopack_rs::MixingRule::Unifac
                        } else {
                            thermopack_rs::MixingRule::VdW
                        },
                        ..Default::default()
                    },
                )?
            };
            let z: Vec<f64> = fields[4].split(',').map(|s| s.parse().unwrap()).collect();
            let spec: Vec<_> = fields[5].split('@').collect();
            let opts: Vec<f64> = spec[0].split(':').map(|s| s.parse().unwrap()).collect();
            let result = model.baseline_phase_envelope(
                &z,
                EnvelopeOptions {
                    start: if opts.get(5).is_some_and(|v| *v > 0.) {
                        EnvelopeConstraint::Temperature(opts[5])
                    } else {
                        EnvelopeConstraint::Pressure(opts[0])
                    },
                    initial_kind: if opts.get(6) == Some(&0.) {
                        SaturationKind::Bubble
                    } else {
                        SaturationKind::Dew
                    },
                    maximum_pressure: opts[1],
                    minimum_temperature: Some(opts[2]),
                    maximum_step: opts[3],
                    maximum_points: opts.get(4).copied().unwrap_or(1000.) as usize,
                },
            );
            let format_point = |p: &MixtureSaturation| {
                [
                    p.temperature,
                    p.pressure,
                    f64::from(p.kind == SaturationKind::Dew),
                    p.liquid_volume,
                    p.vapor_volume,
                ]
                .into_iter()
                .chain(p.liquid.iter().copied())
                .chain(p.vapor.iter().copied())
                .chain([p.dpdt, p.fugacity_residual])
                .map(|v| format!("{v:.17e}"))
                .collect::<Vec<_>>()
                .join(",")
            };
            match result {
                Err(e) => println!("ERROR: {e}"),
                Ok(curve) => {
                    let mut projected = Vec::new();
                    if spec.len() > 1 && !spec[1].is_empty() {
                        for encoded in spec[1].split(';') {
                            let target: Vec<f64> =
                                encoded.split(',').map(|s| s.parse().unwrap()).collect();
                            let target_kind = if target[0] > 0.5 {
                                SaturationKind::Dew
                            } else {
                                SaturationKind::Bubble
                            };
                            let target_sign = if target_kind == SaturationKind::Bubble {
                                1.
                            } else {
                                -1.
                            };
                            let distance = |p: &MixtureSaturation| {
                                let sign = if p.kind == SaturationKind::Bubble {
                                    1.
                                } else {
                                    -1.
                                };
                                ((p.temperature / target[1]).ln() / 0.1).powi(2)
                                    + ((p.pressure / target[2]).ln() / 0.2).powi(2)
                                    + (0..z.len())
                                        .filter(|&i| z[i] > 0.)
                                        .map(|i| {
                                            (sign * (p.vapor[i] / p.liquid[i]).ln()
                                                - target_sign * target[3 + i])
                                                .powi(2)
                                        })
                                        .sum::<f64>()
                            };
                            let mut seed = curve
                                .points
                                .iter()
                                .filter(|p| {
                                    let sign = if p.kind == SaturationKind::Bubble {
                                        1.
                                    } else {
                                        -1.
                                    };
                                    (0..z.len())
                                        .filter(|&i| z[i] > 0.)
                                        .map(|i| {
                                            sign * (p.vapor[i] / p.liquid[i]).ln()
                                                * target_sign
                                                * target[3 + i]
                                        })
                                        .sum::<f64>()
                                        > 0.
                                })
                                .min_by(|a, b| distance(a).total_cmp(&distance(b)))
                                .or_else(|| {
                                    curve
                                        .points
                                        .iter()
                                        .min_by(|a, b| distance(a).total_cmp(&distance(b)))
                                })
                                .unwrap()
                                .clone();
                            if seed.kind != target_kind {
                                std::mem::swap(&mut seed.liquid, &mut seed.vapor);
                                std::mem::swap(&mut seed.liquid_volume, &mut seed.vapor_volume);
                                seed.kind = target_kind;
                            }
                            let component = (0..z.len())
                                .filter(|&i| z[i] > 0.)
                                .max_by(|&i, &j| {
                                    target[3 + i].abs().total_cmp(&target[3 + j].abs())
                                })
                                .unwrap();
                            let constraint = if let Some(&index) = target.get(3 + z.len()) {
                                let index = index as usize;
                                if index < z.len() {
                                    EnvelopeConstraint::LogK {
                                        component: index,
                                        value: target[3 + index],
                                    }
                                } else if index == z.len() {
                                    EnvelopeConstraint::Temperature(target[1])
                                } else {
                                    EnvelopeConstraint::Pressure(target[2])
                                }
                            } else if target[3 + component].abs() < 0.5 {
                                EnvelopeConstraint::LogK {
                                    component,
                                    value: target[3 + component],
                                }
                            } else if (seed.dpdt * seed.temperature / seed.pressure).abs() > 1. {
                                EnvelopeConstraint::Pressure(target[2])
                            } else {
                                EnvelopeConstraint::Temperature(target[1])
                            };
                            projected.push(match model.refine_envelope_state(&seed, constraint) {
                                Ok(p) => format_point(&p),
                                Err(e) => format!("ERROR: {e}"),
                            });
                        }
                    }
                    println!(
                        "{:?}|{}|{}",
                        curve.end,
                        curve
                            .points
                            .iter()
                            .map(&format_point)
                            .collect::<Vec<_>>()
                            .join(";"),
                        projected.join(";")
                    );
                }
            }
            continue;
        }
        if fields[0] == "WSMIX" {
            let kind = match fields[1] {
                "SRK" => thermopack_rs::CubicKind::SoaveRedlichKwong,
                "RK" => thermopack_rs::CubicKind::RedlichKwong,
                "VDW" => thermopack_rs::CubicKind::VanDerWaals,
                _ => thermopack_rs::CubicKind::PengRobinson,
            };
            let components: Vec<&str> = fields[2].split(',').collect();
            let temperature: f64 = fields[3].parse().unwrap();
            let composition: Vec<f64> = fields[4]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let kij: f64 = fields[5].parse().unwrap();
            let tauij: f64 = fields[6].parse().unwrap();
            let alpha: f64 = fields[7].parse().unwrap();
            use thermopack_rs::wong_sandler::{Fraction, WongSandlerData};
            let mut data = WongSandlerData::new(components.len());
            data.set_kij(0, 1, Fraction::constant(kij)).unwrap();
            data.set_kij(1, 0, Fraction::constant(kij)).unwrap();
            data.set_tauij(0, 1, Fraction::constant(tauij)).unwrap();
            data.set_tauij(1, 0, Fraction::constant(tauij)).unwrap();
            data.set_alphaij(0, 1, alpha).unwrap();
            data.set_alphaij(1, 0, alpha).unwrap();
            let options = thermopack_rs::CubicOptions {
                mixing_rule: thermopack_rs::excess_gibbs::MixingRule::WongSandler,
                wong_sandler: Some(&data),
                ..Default::default()
            };
            let model = thermopack_rs::Cubic::with_options(kind, &components, options)?;
            let coefficients = model.wong_sandler_coefficients(temperature, &composition)?;
            let mut text = vec![
                format!("{:.17e}", coefficients.a),
                format!("{:.17e}", coefficients.a_t),
                format!("{:.17e}", coefficients.a_tt),
            ];
            for value in coefficients.a_i() {
                text.push(format!("{value:.17e}"));
            }
            for value in coefficients.a_it() {
                text.push(format!("{value:.17e}"));
            }
            text.push(format!("{:.17e}", coefficients.b));
            text.push(format!("{:.17e}", coefficients.b_t));
            text.push(format!("{:.17e}", coefficients.b_tt));
            for value in coefficients.b_i() {
                text.push(format!("{value:.17e}"));
            }
            for value in coefficients.b_it() {
                text.push(format!("{value:.17e}"));
            }
            println!("{}", text.join(","));
            continue;
        }
        if fields[0] == "WSSTATE" {
            let kind = match fields[1] {
                "SRK" => thermopack_rs::CubicKind::SoaveRedlichKwong,
                "VDW" => thermopack_rs::CubicKind::VanDerWaals,
                _ => thermopack_rs::CubicKind::PengRobinson,
            };
            let components: Vec<&str> = fields[2].split(',').collect();
            let temperature: f64 = fields[3].parse().unwrap();
            let volume: f64 = fields[4].parse().unwrap();
            let composition: Vec<f64> = fields[5]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let kij: f64 = fields[6].parse().unwrap();
            let tauij: f64 = fields[7].parse().unwrap();
            let alpha: f64 = fields[8].parse().unwrap();
            use thermopack_rs::wong_sandler::{Fraction, WongSandlerData};
            let mut data = WongSandlerData::new(components.len());
            data.set_kij(0, 1, Fraction::constant(kij)).unwrap();
            data.set_kij(1, 0, Fraction::constant(kij)).unwrap();
            data.set_tauij(0, 1, Fraction::constant(tauij)).unwrap();
            data.set_tauij(1, 0, Fraction::constant(tauij)).unwrap();
            data.set_alphaij(0, 1, alpha).unwrap();
            data.set_alphaij(1, 0, alpha).unwrap();
            let options = thermopack_rs::CubicOptions {
                mixing_rule: thermopack_rs::excess_gibbs::MixingRule::WongSandler,
                wong_sandler: Some(&data),
                ..Default::default()
            };
            let model = thermopack_rs::Cubic::with_options(kind, &components, options)?;
            let state = model.prepare(temperature, &composition)?;
            let pressure = state.pressure(volume)?;
            let mut fugacity = vec![0.0; components.len()];
            state.ln_fugacity_into(volume, &mut fugacity)?;
            let mut dt = vec![0.0; components.len()];
            let mut dp = vec![0.0; components.len()];
            let mut dn = vec![0.0; components.len() * components.len()];
            state.fugacity_derivatives_into(volume, &mut dt, &mut dp, &mut dn)?;
            let mut text = vec![
                format!("{:.17e}", pressure.value),
                format!("{:.17e}", pressure.dv),
                format!("{:.17e}", pressure.dt),
                format!("{:.17e}", state.pressure_second_derivative(volume)?),
                format!("{:.17e}", state.enthalpy(volume)?),
                format!("{:.17e}", state.entropy(volume)?),
            ];
            for value in &fugacity {
                text.push(format!("{value:.17e}"));
            }
            for value in &dt {
                text.push(format!("{value:.17e}"));
            }
            for value in &dp {
                text.push(format!("{value:.17e}"));
            }
            for value in &dn {
                text.push(format!("{value:.17e}"));
            }
            println!("{}", text.join(","));
            continue;
        }
        if fields[0] == "CPARDF" {
            let component = fields[1];
            let temperature: f64 = fields[2].parse().unwrap();
            let volume: f64 = fields[3].parse().unwrap();
            let moles: Vec<f64> = fields[4]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let model = thermopack_rs::Cpa::new(&[component])?;
            let rdf = thermopack_rs::rdf::cpa_rdf(&model, temperature, volume, &moles)?;
            let mut text = vec![
                format!("{:.17e}", rdf.g),
                format!("{:.17e}", rdf.temperature),
                format!("{:.17e}", rdf.volume),
            ];
            for value in &rdf.moles {
                text.push(format!("{value:.17e}"));
            }
            text.push(format!("{:.17e}", rdf.temperature_second));
            text.push(format!("{:.17e}", rdf.temperature_volume));
            for value in &rdf.temperature_moles {
                text.push(format!("{value:.17e}"));
            }
            text.push(format!("{:.17e}", rdf.volume_second));
            for value in &rdf.volume_moles {
                text.push(format!("{value:.17e}"));
            }
            for value in &rdf.moles_second {
                text.push(format!("{value:.17e}"));
            }
            println!("{}", text.join(","));
            continue;
        }
        if fields[0] == "PCSAFTRDF" {
            let components: Vec<&str> = fields[1].split(',').collect();
            let temperature: f64 = fields[2].parse().unwrap();
            let volume: f64 = fields[3].parse().unwrap();
            let moles: Vec<f64> = fields[4]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let first: usize = fields[5].parse().unwrap();
            let second: usize = fields[6].parse().unwrap();
            let model = thermopack_rs::PcSaft::new(&components)?;
            let rdf =
                thermopack_rs::rdf::pcsaft_rdf(&model, temperature, volume, &moles, first, second)?;
            let mut text = vec![
                format!("{:.17e}", rdf.g),
                format!("{:.17e}", rdf.temperature),
                format!("{:.17e}", rdf.volume),
            ];
            for value in &rdf.moles {
                text.push(format!("{value:.17e}"));
            }
            text.push(format!("{:.17e}", rdf.temperature_second));
            text.push(format!("{:.17e}", rdf.temperature_volume));
            for value in &rdf.temperature_moles {
                text.push(format!("{value:.17e}"));
            }
            text.push(format!("{:.17e}", rdf.volume_second));
            for value in &rdf.volume_moles {
                text.push(format!("{value:.17e}"));
            }
            for value in &rdf.moles_second {
                text.push(format!("{value:.17e}"));
            }
            println!("{}", text.join(","));
            continue;
        }
        if fields[0] == "PCSAFTHS" {
            let components: Vec<&str> = fields[1].split(',').collect();
            let temperature: f64 = fields[2].parse().unwrap();
            let volume: f64 = fields[3].parse().unwrap();
            let moles: Vec<f64> = fields[4]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let model = thermopack_rs::PcSaft::new(&components)?;
            let energy =
                thermopack_rs::pc_saft::pcsaft_hs_energy(&model, temperature, volume, &moles)?;
            let mut text = vec![
                format!("{:.17e}", energy.value),
                format!("{:.17e}", energy.d_t),
                format!("{:.17e}", energy.d_v),
            ];
            for value in &energy.d_n {
                text.push(format!("{value:.17e}"));
            }
            text.push(format!("{:.17e}", energy.d_tt));
            text.push(format!("{:.17e}", energy.d_tv));
            for value in &energy.d_tn {
                text.push(format!("{value:.17e}"));
            }
            text.push(format!("{:.17e}", energy.d_vv));
            for value in &energy.d_vn {
                text.push(format!("{value:.17e}"));
            }
            for row in &energy.d_nn {
                for value in row {
                    text.push(format!("{value:.17e}"));
                }
            }
            println!("{}", text.join(","));
            continue;
        }
        if fields[0] == "PCSAFTCHAIN" {
            let components: Vec<&str> = fields[1].split(',').collect();
            let temperature: f64 = fields[2].parse().unwrap();
            let volume: f64 = fields[3].parse().unwrap();
            let moles: Vec<f64> = fields[4]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let model = thermopack_rs::PcSaft::new(&components)?;
            let energy =
                thermopack_rs::pc_saft::pcsaft_chain_energy(&model, temperature, volume, &moles)?;
            let mut text = vec![
                format!("{:.17e}", energy.value),
                format!("{:.17e}", energy.d_t),
                format!("{:.17e}", energy.d_v),
            ];
            for value in &energy.d_n {
                text.push(format!("{value:.17e}"));
            }
            text.push(format!("{:.17e}", energy.d_tt));
            text.push(format!("{:.17e}", energy.d_tv));
            for value in &energy.d_tn {
                text.push(format!("{value:.17e}"));
            }
            text.push(format!("{:.17e}", energy.d_vv));
            for value in &energy.d_vn {
                text.push(format!("{value:.17e}"));
            }
            for row in &energy.d_nn {
                for value in row {
                    text.push(format!("{value:.17e}"));
                }
            }
            println!("{}", text.join(","));
            continue;
        }
        if fields[0] == "PCSAFTSTATE" {
            let components: Vec<&str> = fields[1].split(',').collect();
            let temperature: f64 = fields[2].parse().unwrap();
            let volume: f64 = fields[3].parse().unwrap();
            let moles: Vec<f64> = fields[4]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let total: f64 = moles.iter().sum();
            let composition: Vec<f64> = moles.iter().map(|amount| amount / total).collect();
            let rho = total / volume;
            let model = thermopack_rs::PcSaft::new(&components)?;
            let prepared = model.prepare(temperature, &composition)?;
            let alpha = prepared.residual_helmholtz_energy(rho);
            let compressibility = prepared.compressibility_factor(rho);
            let alpha_rho = (compressibility - 1.0) / rho;
            let state = prepared.state(rho)?;
            let alpha_t =
                -state.u_residual / (thermopack_rs::pc_saft::RGAS * temperature * temperature);
            let (alpha_rhorho, alpha_rhot) = prepared.residual_helmholtz_density_derivatives(rho);
            println!(
                "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                alpha, alpha_rho, alpha_t, alpha_rhorho, alpha_rhot
            );
            continue;
        }
        if fields[0] == "CRITICAL" {
            let kind = match fields[1] {
                "SRK" => thermopack_rs::CubicKind::SoaveRedlichKwong,
                "RK" => thermopack_rs::CubicKind::RedlichKwong,
                "VDW" => thermopack_rs::CubicKind::VanDerWaals,
                _ => thermopack_rs::CubicKind::PengRobinson,
            };
            let components: Vec<&str> = fields[2].split(',').collect();
            let composition: Vec<f64> = fields[3]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let config = fields.get(4).copied().unwrap_or("classic");
            let ws = thermopack_rs::wong_sandler::WongSandlerData::new(components.len());
            let model = thermopack_rs::Cubic::with_options(
                kind,
                &components,
                thermopack_rs::CubicOptions {
                    volume_translation: config == "translated",
                    wong_sandler: (config == "ws").then_some(&ws),
                    mixing_rule: if config == "ws" {
                        thermopack_rs::MixingRule::WongSandler
                    } else {
                        thermopack_rs::MixingRule::VdW
                    },
                    ..Default::default()
                },
            )?;
            let result = if fields.len() >= 7 {
                thermopack_rs::mixture_critical::critical_point_with_guess(
                    &model,
                    &composition,
                    fields[5].parse()?,
                    fields[6].parse()?,
                )
            } else {
                thermopack_rs::critical::critical_point(&model, &composition)
            };
            match result {
                Ok(point) => println!(
                    "{:.17e},{:.17e},{:.17e}",
                    point.temperature, point.volume, point.pressure
                ),
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if fields[0] == "LJSPOT" {
            let separations: Vec<f64> = fields[1]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let model = thermopack_rs::ljs::LjsBh::new();
            let (potential, force) = model.reduced_potential(&separations);
            let text: Vec<String> = potential
                .iter()
                .zip(force.iter())
                .map(|(value, slope)| format!("{value:.17e},{slope:.17e}"))
                .collect();
            println!("{}", text.join(";"));
            continue;
        }
        if fields[0] == "LJSVOLUME" {
            let temperature: f64 = fields[1].parse().unwrap();
            let pressure: f64 = fields[2].parse().unwrap();
            let phase = match fields[3].parse::<i32>().unwrap() {
                2 => Phase::Vapor,
                _ => Phase::Liquid,
            };
            let model = thermopack_rs::ljs::LjsBh::new();
            let volume = model.volume(temperature, pressure, phase)?;
            println!("{:.17e}", volume);
            continue;
        }
        if fields[0] == "LJSCALORIC" {
            let temperature: f64 = fields[1].parse().unwrap();
            let volume: f64 = fields[2].parse().unwrap();
            let moles: f64 = fields[3].parse().unwrap();
            let model = thermopack_rs::ljs::LjsBh::new();
            let state = model.state(temperature, volume, moles)?;
            println!("{:.17e},{:.17e}", state.h_total, state.s_total);
            continue;
        }
        if fields[0] == "QUAD" {
            let index: usize = fields[1].parse().unwrap();
            let rules = [
                thermopack_rs::quadrature::Quadrature::Gauss10,
                thermopack_rs::quadrature::Quadrature::Kronrod21,
                thermopack_rs::quadrature::Quadrature::Gauss15,
                thermopack_rs::quadrature::Quadrature::Kronrod31,
                thermopack_rs::quadrature::Quadrature::Gauss20,
                thermopack_rs::quadrature::Quadrature::Kronrod41,
                thermopack_rs::quadrature::Quadrature::Gauss25,
                thermopack_rs::quadrature::Quadrature::Kronrod51,
                thermopack_rs::quadrature::Quadrature::Gauss30,
                thermopack_rs::quadrature::Quadrature::Kronrod61,
                thermopack_rs::quadrature::Quadrature::Gauss40,
                thermopack_rs::quadrature::Quadrature::Kronrod81,
                thermopack_rs::quadrature::Quadrature::Gauss50,
                thermopack_rs::quadrature::Quadrature::Kronrod101,
            ];
            let rule = rules[index];
            let positions = rule.positions();
            let weights = rule.weights();
            let (integral, differential) = thermopack_rs::quadrature::integrate_with_differential(
                rule,
                |point| Ok((point * point, 2.0 * point, 2.0)),
                0.0,
                1.0,
                0.0,
                1.0,
            )?;
            let mut text: Vec<String> = positions
                .iter()
                .chain(weights.iter())
                .map(|value| format!("{value:.17e}"))
                .collect();
            text.push(format!("{integral:.17e}"));
            text.push(format!("{differential:.17e}"));
            println!("{}", text.join(","));
            continue;
        }
        if fields[0] == "AISLJS" {
            let index: usize = fields[1].parse().unwrap();
            let reduced_distance: f64 = fields[2].parse().unwrap();
            let reduced_density: f64 = fields[3].parse().unwrap();
            let value = thermopack_rs::ljs::ais_ljs_ex(index, reduced_distance, reduced_density)?;
            let correction = thermopack_rs::ljs::chi_ljs_ex(reduced_distance, reduced_density)?;
            println!("{value:.17e},{correction:.17e}");
            continue;
        }
        if fields[0] == "CHAINRULE" {
            use thermopack_rs::chain_rule::{
                PackingFractionDerivatives, TermDerivatives, convert_zeta_x_to_tvn,
            };
            let zeta = PackingFractionDerivatives {
                value: 0.1,
                temperature: -0.2,
                temperature_second: 0.3,
                volume: -0.4,
                volume_second: 0.5,
                volume_third: 0.0,
                temperature_volume: 0.6,
                volume_temperature_second: 0.0,
                volume_temperature_third: 0.0,
                moles: vec![-0.7, 0.8],
                moles_second: vec![0.9, -1.1, 1.2, -1.3],
                temperature_moles: vec![0.0, 0.0],
                volume_moles: vec![0.0, 0.0],
                volume_moles_second: vec![0.0; 4],
                volume_temperature_moles: vec![0.0, 0.0],
            };
            let term = TermDerivatives {
                value: 2.0,
                packing_fraction: 0.25,
                reduced_distance: -0.35,
                packing_fraction_second: 0.45,
                reduced_distance_second: -0.55,
                packing_fraction_reduced_distance: 0.65,
                packing_fraction_third: -0.75,
                packing_fraction_second_reduced_distance: 0.85,
                packing_fraction_reduced_distance_second: -0.95,
            };
            let result = convert_zeta_x_to_tvn(&term, 0.3, 0.15, -0.25, &zeta, 2);
            let values = [
                result.temperature,
                result.volume,
                result.moles[0],
                result.moles[1],
                result.temperature_second,
                result.volume_second,
                result.temperature_volume,
                result.temperature_moles[0],
                result.temperature_moles[1],
                result.volume_moles[0],
                result.volume_moles[1],
                result.moles_second[0],
                result.moles_second[1],
                result.moles_second[2],
                result.moles_second[3],
            ];
            let text: Vec<String> = values.iter().map(|value| format!("{value:.17e}")).collect();
            println!("{}", text.join(","));
            continue;
        }
        if fields[0] == "UVVIRIALDIFF" {
            let sigma: f64 = fields[1].parse().unwrap();
            let eps_divk: f64 = fields[2].parse().unwrap();
            let temperature: f64 = fields[3].parse().unwrap();
            let rho: f64 = fields[4].parse().unwrap();
            let (diameter, diameter_t, diameter_tt) =
                thermopack_rs::hardsphere_wca::hard_sphere_diameter(sigma, eps_divk, temperature);
            let (b2, b2_r, b2_t, b2_rr, b2_tt, b2_rt) = thermopack_rs::ljs::uv_delta_b2_overall(
                sigma,
                eps_divk,
                temperature,
                diameter,
                diameter_t,
                diameter_tt,
                rho,
            );
            let (a, a_r, a_t, a_rr, a_tt, a_rt) = thermopack_rs::ljs::calc_uv_wca(
                sigma,
                eps_divk,
                rho,
                temperature,
                diameter,
                diameter_t,
                diameter_tt,
                true,
                false,
            );
            println!(
                "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                b2, b2_r, b2_t, b2_rr, b2_tt, b2_rt, a, a_r, a_t, a_rr, a_tt, a_rt
            );
            continue;
        }
        if fields[0] == "WCATERM" {
            let index: usize = fields[1].parse().unwrap();
            let eps_divk: f64 = fields[2].parse().unwrap();
            let reduced_distance: f64 = fields[3].parse().unwrap();
            let reduced_density: f64 = fields[4].parse().unwrap();
            let (value, value_r, value_y, value_rr, value_yy, value_ry) =
                thermopack_rs::ljs_wca::wca_term(
                    index,
                    eps_divk,
                    reduced_distance,
                    reduced_density,
                )?;
            println!(
                "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                value, value_r, value_y, value_rr, value_yy, value_ry
            );
            continue;
        }
        if fields[0] == "UVVIRIALS" {
            let eps_divk: f64 = fields[1].parse().unwrap();
            let temperature: f64 = fields[2].parse().unwrap();
            let (q, q_t, q_tt) = thermopack_rs::ljs::uv_q_of_temperature(eps_divk, temperature);
            let (hs2, hs2_t, hs2_tt) = thermopack_rs::ljs::uv_hard_sphere_b2(eps_divk, temperature);
            let (ljs2, ljs2_t, ljs2_tt) = thermopack_rs::ljs::ljs_b2(eps_divk, temperature);
            println!(
                "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                q, q_t, q_tt, hs2, hs2_t, hs2_tt, ljs2, ljs2_t, ljs2_tt
            );
            continue;
        }
        if fields[0] == "UVTERMS" {
            let sigma: f64 = fields[1].parse().unwrap();
            let eps_divk: f64 = fields[2].parse().unwrap();
            let temperature: f64 = fields[3].parse().unwrap();
            let rho: f64 = fields[4].parse().unwrap();
            let (diameter, diameter_t, diameter_tt) =
                thermopack_rs::hardsphere_wca::hard_sphere_diameter(sigma, eps_divk, temperature);
            let (a1, a1_r, a1_t, a1_rr, a1_tt, a1_rt) = thermopack_rs::ljs::uv_a1_u(
                sigma,
                eps_divk,
                temperature,
                diameter,
                diameter_t,
                diameter_tt,
                rho,
            );
            let (b2, b2_t, b2_tt) = thermopack_rs::ljs::uv_a1_b2(
                sigma,
                eps_divk,
                temperature,
                diameter,
                diameter_t,
                diameter_tt,
            );
            println!(
                "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                a1, a1_r, a1_t, a1_rr, a1_tt, a1_rt, b2, b2_t, b2_tt
            );
            continue;
        }
        if fields[0] == "UVPHI" {
            let eps_divk: f64 = fields[1].parse().unwrap();
            let temperature: f64 = fields[2].parse().unwrap();
            let rho: f64 = fields[3].parse().unwrap();
            let (value, value_r, value_t, value_rr, value_tt, value_rt) =
                thermopack_rs::ljs::uv_phi(eps_divk, temperature, rho);
            println!(
                "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                value, value_r, value_t, value_rr, value_tt, value_rt
            );
            continue;
        }
        if fields[0] == "WHITEBEAR" {
            let temperature: f64 = fields[1].parse().unwrap();
            let eps_divk: f64 = fields[2].parse().unwrap();
            let eta: f64 = fields[3].parse().unwrap();
            let r_star: f64 = fields[4].parse().unwrap();
            let r_star_t: f64 = fields[5].parse().unwrap();
            let r_star_tt: f64 = fields[6].parse().unwrap();
            let (r_min, r_min_t, r_min_tt) =
                thermopack_rs::hardsphere_wca::minimum_cavity_radius(temperature, eps_divk);
            let distribution = thermopack_rs::hardsphere_wca::cavity_distribution_function(
                eta, r_star, r_star_t, r_star_tt,
            );
            let (e, e_t, e_tt) = thermopack_rs::hardsphere_wca::lj_boltzmann_factor(
                r_star,
                r_star_t,
                r_star_tt,
                temperature,
                eps_divk,
            );
            let (diameter, diameter_t, diameter_tt) =
                thermopack_rs::hardsphere_wca::hard_sphere_diameter(1.0, eps_divk, temperature);
            let values = [
                r_min,
                r_min_t,
                r_min_tt,
                distribution.value,
                distribution.packing_fraction,
                distribution.packing_fraction_second,
                distribution.temperature,
                distribution.temperature_second,
                distribution.packing_fraction_temperature,
                e,
                e_t,
                e_tt,
                diameter,
                diameter_t,
                diameter_tt,
            ];
            let text: Vec<String> = values.iter().map(|value| format!("{value:.17e}")).collect();
            println!("{}", text.join(","));
            continue;
        }
        if fields[0] == "LJSWCATERM" {
            let index: usize = fields[1].parse().unwrap();
            let reduced_temperature: f64 = fields[2].parse().unwrap();
            let reduced_density: f64 = fields[3].parse().unwrap();
            let (value, _, _, _, _, _) = thermopack_rs::ljs_wca::perturbation_term(
                index,
                reduced_temperature,
                reduced_density,
            )?;
            println!("{value:.17e}");
            continue;
        }
        if fields[0] == "LJSFUGACITY" {
            let temperature: f64 = fields[1].parse().unwrap();
            let volume: f64 = fields[2].parse().unwrap();
            let moles: f64 = fields[3].parse().unwrap();
            let model = thermopack_rs::ljs::LjsBh::new();
            let coefficient = model.fugacity_coefficient(temperature, volume, moles)?;
            let (pressure, _) = model.pressure(temperature, volume, moles)?;
            println!("{:.17e},{:.17e}", coefficient, coefficient + pressure.ln());
            continue;
        }
        if fields[0] == "LJSSTATE" {
            let temperature: f64 = fields[1].parse().unwrap();
            let volume: f64 = fields[2].parse().unwrap();
            let moles: f64 = fields[3].parse().unwrap();
            let model = thermopack_rs::ljs::LjsBh::new();
            let (pressure, compressibility) = model.pressure(temperature, volume, moles)?;
            println!("{:.17e},{:.17e}", pressure, compressibility);
            continue;
        }
        if fields[0] == "LJSHELMHOLTZ" {
            let temperature: f64 = fields[1].parse().unwrap();
            let volume: f64 = fields[2].parse().unwrap();
            let moles: f64 = fields[3].parse().unwrap();
            let mut model = thermopack_rs::ljs::LjsBh::new();
            if fields.len() > 4 && fields[4] == "1" {
                model.set_control(thermopack_rs::ljs::LjsControl {
                    use_lafitte_a3: true,
                    ..thermopack_rs::ljs::LjsControl::default()
                });
            }
            let (energy, derivative) =
                model.residual_helmholtz_energy_derivative(temperature, volume, moles)?;
            println!("{:.17e},{:.17e}", energy * moles, derivative * moles);
            continue;
        }
        if fields[0] == "LJSSPLINE" {
            let reduced_temperature: f64 = fields[1].parse().unwrap();
            let reduced_density: f64 = fields[2].parse().unwrap();
            let index: usize = fields[3].parse().unwrap();
            let model = thermopack_rs::ljs::LjsBh::new();
            let (value, value_r, value_x, value_rr, value_xx, value_rx) =
                model.spline_polynomial(index, reduced_temperature, reduced_density)?;
            println!(
                "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                value, value_r, value_x, value_rr, value_xx, value_rx
            );
            continue;
        }
        if fields[0] == "LJS" {
            let reduced_temperature: f64 = fields[1].parse().unwrap();
            let reduced_density: f64 = fields[2].parse().unwrap();
            let model = thermopack_rs::ljs::LjsBh::new();
            let (first, second, third) =
                model.reduced_perturbation_terms(reduced_temperature, reduced_density)?;
            let diameter = model.diameter_div_sigma(reduced_temperature)?;
            println!(
                "{:.17e},{:.17e},{:.17e},{:.17e}",
                first, second, third, diameter
            );
            continue;
        }
        if fields[0] == "THIRDVIRIALMATRIX" {
            let kind = match fields[1] {
                "SRK" => thermopack_rs::CubicKind::SoaveRedlichKwong,
                "RK" => thermopack_rs::CubicKind::RedlichKwong,
                "VDW" => thermopack_rs::CubicKind::VanDerWaals,
                _ => thermopack_rs::CubicKind::PengRobinson,
            };
            let components: Vec<&str> = fields[2].split(',').collect();
            let temperature: f64 = fields[3].parse().unwrap();
            let model = thermopack_rs::Cubic::new(kind, &components)?;
            let matrix = model.binary_third_virial_coefficient_matrix(temperature)?;
            let text: Vec<String> = matrix.iter().map(|value| format!("{value:.17e}")).collect();
            println!("{}", text.join(","));
            continue;
        }
        if fields[0] == "VIRIALMATRIX" {
            let kind = match fields[1] {
                "SRK" => thermopack_rs::CubicKind::SoaveRedlichKwong,
                "RK" => thermopack_rs::CubicKind::RedlichKwong,
                "VDW" => thermopack_rs::CubicKind::VanDerWaals,
                _ => thermopack_rs::CubicKind::PengRobinson,
            };
            let components: Vec<&str> = fields[2].split(',').collect();
            let temperature: f64 = fields[3].parse().unwrap();
            let model = thermopack_rs::Cubic::new(kind, &components)?;
            let matrix = model.second_virial_coefficient_matrix(temperature)?;
            let text: Vec<String> = matrix.iter().map(|value| format!("{value:.17e}")).collect();
            println!("{}", text.join(","));
            continue;
        }
        if fields[0] == "VIRIAL" {
            let kind = match fields[1] {
                "SRK" => thermopack_rs::CubicKind::SoaveRedlichKwong,
                "RK" => thermopack_rs::CubicKind::RedlichKwong,
                "VDW" => thermopack_rs::CubicKind::VanDerWaals,
                _ => thermopack_rs::CubicKind::PengRobinson,
            };
            let components: Vec<&str> = fields[2].split(',').collect();
            let temperature: f64 = fields[3].parse().unwrap();
            let composition: Vec<f64> = fields[4]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let model = thermopack_rs::Cubic::new(kind, &components)?;
            let state = model.prepare(temperature, &composition)?;
            let (second, third) = state.virial_coefficients()?;
            println!("{second:.17e},{third:.17e}");
            continue;
        }
        if fields[0] == "JTI" {
            let kind = match fields[1] {
                "SRK" => thermopack_rs::CubicKind::SoaveRedlichKwong,
                "RK" => thermopack_rs::CubicKind::RedlichKwong,
                "VDW" => thermopack_rs::CubicKind::VanDerWaals,
                _ => thermopack_rs::CubicKind::PengRobinson,
            };
            let components: Vec<&str> = fields[2].split(',').collect();
            let volume: f64 = fields[3].parse().unwrap();
            let guess: f64 = fields[4].parse().unwrap();
            let model = thermopack_rs::Cubic::new(kind, &components)?;
            let composition = vec![1.0; components.len()];
            let point = thermopack_rs::joule_thompson::inversion_point(
                &model,
                volume,
                guess,
                &composition,
            )?;
            println!(
                "{:.17e},{:.17e},{:.17e}",
                point.temperature, point.volume, point.pressure
            );
            continue;
        }
        if fields[0] == "SPINODAL" {
            let kind = match fields[1] {
                "SRK" => thermopack_rs::CubicKind::SoaveRedlichKwong,
                "RK" => thermopack_rs::CubicKind::RedlichKwong,
                "VDW" => thermopack_rs::CubicKind::VanDerWaals,
                _ => thermopack_rs::CubicKind::PengRobinson,
            };
            let components: Vec<&str> = fields[2].split(',').collect();
            let pressure: f64 = fields[3].parse().unwrap();
            let phase = match fields[4].parse::<i32>().unwrap() {
                2 => Phase::Vapor,
                _ => Phase::Liquid,
            };
            let model = thermopack_rs::Cubic::new(kind, &components)?;
            let composition = vec![1.0; components.len()];
            let point =
                thermopack_rs::spinodal::spinodal_point(&model, pressure, phase, &composition)?;
            println!(
                "{:.17e},{:.17e},{:.17e}",
                point.temperature, point.volume, point.pressure
            );
            continue;
        }
        if fields[0] == "MIXSPIN" {
            use thermopack_rs::mixture_spinodal::*;
            let kind = match fields[1] {
                "SRK" => CubicKind::SoaveRedlichKwong,
                "RK" => CubicKind::RedlichKwong,
                "VDW" => CubicKind::VanDerWaals,
                _ => CubicKind::PengRobinson,
            };
            let names: Vec<_> = fields[2].split(',').collect();
            let z: Vec<f64> = fields[3]
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?;
            let model = Cubic::with_options(
                kind,
                &names,
                thermopack_rs::CubicOptions {
                    volume_translation: fields[4] == "translated",
                    ..Default::default()
                },
            )?;
            let spec: f64 = fields[6].parse()?;
            let result = match fields[5] {
                "BASE_L" | "BASE_V" => baseline_mixture_spinodal_point(
                    &model,
                    spec,
                    if fields[5] == "BASE_L" {
                        Phase::Liquid
                    } else {
                        Phase::Vapor
                    },
                    &z,
                )
                .map(|p| vec![p]),
                "L" | "V" => thermopack_rs::spinodal::spinodal_point(
                    &model,
                    spec,
                    if fields[5] == "L" {
                        Phase::Liquid
                    } else {
                        Phase::Vapor
                    },
                    &z,
                )
                .map(|p| vec![p]),
                "T" => spinodal_volume_at_temperature(&model, spec, &z, fields[7].parse()?)
                    .map(|p| vec![p]),
                "VOLUME" => spinodal_temperature_at_volume(&model, spec, &z, fields[7].parse()?)
                    .map(|p| vec![p]),
                "CURVE" | "BASE_CURVE" => (if fields[5] == "BASE_CURVE" {
                    baseline_trace_spinodal
                } else {
                    trace_spinodal
                })(
                    &model,
                    &z,
                    SpinodalTraceOptions {
                        initial_pressure: spec,
                        log_volume_step: fields[7].parse()?,
                        minimum_vapor_temperature: fields.get(8).map(|s| s.parse()).transpose()?,
                        ..Default::default()
                    },
                )
                .map(|c| c.points),
                _ => Err(thermopack_rs::Error::InvalidInput("spinodal mode")),
            };
            match result {
                Ok(points) => println!(
                    "{}",
                    points
                        .iter()
                        .flat_map(|p| [p.temperature, p.volume, p.pressure])
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
                Err(error) => println!("ERROR: {error}"),
            }
            continue;
        }
        if fields[0] == "META" {
            let kind = match fields[1] {
                "SRK" => CubicKind::SoaveRedlichKwong,
                "RK" => CubicKind::RedlichKwong,
                "VDW" => CubicKind::VanDerWaals,
                _ => CubicKind::PengRobinson,
            };
            let names: Vec<_> = fields[2].split(',').collect();
            let model = if fields[3] == "quantum" {
                Cubic::quantum(&names)?
            } else {
                Cubic::with_options(
                    kind,
                    &names,
                    thermopack_rs::CubicOptions {
                        volume_translation: fields[3] == "translated",
                        mixing_rule: match fields[3] {
                            "unifac" => thermopack_rs::MixingRule::Unifac,
                            "hv1" => thermopack_rs::MixingRule::HuronVidal1,
                            _ => thermopack_rs::MixingRule::VdW,
                        },
                        ..Default::default()
                    },
                )?
            };
            let z: Vec<f64> = fields[5]
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?;
            let cfg: Vec<_> = fields[4].split(':').collect();
            let p = cfg[1].parse()?;
            let s = cfg[2].parse()?;
            let result: thermopack_rs::Result<Vec<f64>> = if cfg[0] == "PS" {
                model
                    .homogeneous_ps(p, s, &z, cfg[3].parse()?, cfg[4].parse()?)
                    .map(|p| {
                        vec![
                            p.temperature,
                            p.volume,
                            p.pressure,
                            p.entropy,
                            p.pressure_residual,
                            p.entropy_residual,
                            p.iterations as f64,
                        ]
                    })
            } else {
                model
                    .baseline_metastable_isentrope(&z, p, s, cfg[3].parse()?, cfg[4].parse()?)
                    .map(|c| {
                        c.points
                            .iter()
                            .flat_map(|p| [p.temperature, p.volume, p.pressure])
                            .collect()
                    })
            };
            match result {
                Ok(v) => println!(
                    "{}",
                    v.iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
                Err(e) => println!("ERROR: {e}"),
            };
            continue;
        }
        if fields[0] == "SPINPROP" {
            use thermopack_rs::spinodal_properties::*;
            let kind = match fields[1] {
                "SRK" => CubicKind::SoaveRedlichKwong,
                "RK" => CubicKind::RedlichKwong,
                "VDW" => CubicKind::VanDerWaals,
                _ => CubicKind::PengRobinson,
            };
            let names: Vec<_> = fields[2].split(',').collect();
            let mut z = vec![0.0; names.len()];
            z[0] = 1.0;
            let model = Cubic::with_options(
                kind,
                &names,
                thermopack_rs::CubicOptions {
                    volume_translation: fields[3] == "translated",
                    ..Default::default()
                },
            )?;
            let config: Vec<_> = fields[4].split(':').collect();
            let p = fields[5].parse()?;
            let result: thermopack_rs::Result<Vec<f64>> = if config[0] == "GRID" {
                baseline_pure_spinodal_curve(
                    &model,
                    &z,
                    p,
                    None,
                    config.get(1).map(|s| s.parse()).transpose()?,
                )
                .map(|c| {
                    c.points
                        .iter()
                        .flat_map(|p| [p.temperature, p.volume, p.pressure])
                        .collect()
                })
            } else {
                let property = match config[1] {
                    "1" => SpinodalProperty::Entropy,
                    "2" => SpinodalProperty::Volume,
                    "3" => SpinodalProperty::Enthalpy,
                    "4" => SpinodalProperty::Temperature,
                    "5" => SpinodalProperty::Pressure,
                    _ => SpinodalProperty::InternalEnergy,
                };
                if config[0] == "EXTREMA" {
                    baseline_pure_spinodal_extrema(&model, &z, property, p).map(|e| {
                        [e.minimum, e.maximum]
                            .into_iter()
                            .flat_map(|p| [p.temperature, p.volume, p.pressure])
                            .chain([f64::from(e.ordering_ambiguous)])
                            .collect()
                    })
                } else {
                    let target = config[3].parse()?;
                    let phase = match config[2] {
                        "1" => Some(Phase::Liquid),
                        "2" => Some(Phase::Vapor),
                        _ => None,
                    };
                    baseline_pure_spinodal_intersection(
                        &model, &z, property, target, p, phase, None,
                    )
                    .map(|found| match found {
                        Some(p) => vec![1.0, p.temperature, p.volume, p.pressure],
                        None => vec![0.0, 0.0, 0.0, 0.0],
                    })
                }
            };
            match result {
                Ok(values) => println!(
                    "{}",
                    values
                        .iter()
                        .map(|v| format!("{v:.17e}"))
                        .collect::<Vec<_>>()
                        .join(",")
                ),
                Err(e) => println!("ERROR: {e}"),
            };
            continue;
        }
        if fields[0] == "ISOLINE" {
            let kind = match fields[1] {
                "SRK" => thermopack_rs::CubicKind::SoaveRedlichKwong,
                "RK" => thermopack_rs::CubicKind::RedlichKwong,
                "VDW" => thermopack_rs::CubicKind::VanDerWaals,
                _ => thermopack_rs::CubicKind::PengRobinson,
            };
            let components: Vec<&str> = fields[2].split(',').collect();
            let composition: Vec<f64> = fields[3]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let specification: f64 = fields[4].parse().unwrap();
            let minimum: f64 = fields[5].parse().unwrap();
            let maximum: f64 = fields[6].parse().unwrap();
            let (lower, upper, points) = if fields.len() > 9 {
                (
                    fields[7].parse().unwrap(),
                    fields[8].parse().unwrap(),
                    fields[9].parse().unwrap(),
                )
            } else {
                (0.0, 0.0, fields[7].parse().unwrap())
            };
            let model = thermopack_rs::Cubic::new(kind, &components)?;
            let line = match fields[fields.len() - 1] {
                "isobar" => thermopack_rs::isolines::isobar(
                    &model,
                    specification,
                    minimum,
                    maximum,
                    &composition,
                    points,
                )?,
                "isenthalp" => thermopack_rs::isolines::isenthalp(
                    &model,
                    thermopack_rs::isolines::CaloricLine {
                        specification,
                        minimum_pressure: minimum,
                        maximum_pressure: maximum,
                        minimum_temperature: lower,
                        maximum_temperature: upper,
                    },
                    &composition,
                    points,
                )?,
                "isentrope" => thermopack_rs::isolines::isentrope(
                    &model,
                    thermopack_rs::isolines::CaloricLine {
                        specification,
                        minimum_pressure: minimum,
                        maximum_pressure: maximum,
                        minimum_temperature: lower,
                        maximum_temperature: upper,
                    },
                    &composition,
                    points,
                )?,
                _ => thermopack_rs::isolines::isotherm(
                    &model,
                    specification,
                    minimum,
                    maximum,
                    &composition,
                    points,
                )?,
            };
            let mut text = String::new();
            for index in 0..line.temperature.len() {
                if index > 0 {
                    text.push(';');
                }
                text.push_str(&format!(
                    "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                    line.temperature[index],
                    line.pressure[index],
                    line.volume[index],
                    line.entropy[index],
                    line.enthalpy[index]
                ));
            }
            println!("{text}");
            continue;
        }
        if fields[0] == "PSEUDOC" {
            let kind = match fields[1] {
                "SRK" => thermopack_rs::CubicKind::SoaveRedlichKwong,
                "RK" => thermopack_rs::CubicKind::RedlichKwong,
                "VDW" => thermopack_rs::CubicKind::VanDerWaals,
                _ => thermopack_rs::CubicKind::PengRobinson,
            };
            let components: Vec<&str> = fields[2].split(',').collect();
            let critical_temperatures: Vec<f64> = fields[3]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let critical_pressures: Vec<f64> = fields[4]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let acentric_factors: Vec<f64> = fields[5]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let molar_masses: Option<Vec<f64>> = if fields[6].is_empty() {
                None
            } else {
                Some(
                    fields[6]
                        .split(',')
                        .map(|value| value.parse().unwrap())
                        .collect(),
                )
            };
            let temperature: f64 = fields[7].parse().unwrap();
            let pressure: f64 = fields[8].parse().unwrap();
            let composition: Vec<f64> = fields[9]
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect();
            let phase = match fields[10].parse::<i32>().unwrap() {
                2 => Phase::Vapor,
                _ => Phase::Liquid,
            };
            let model = thermopack_rs::Cubic::new_pseudo(
                kind,
                &components,
                &critical_temperatures,
                &critical_pressures,
                &acentric_factors,
                molar_masses.as_deref(),
            )?;
            let state = model.prepare(temperature, &composition)?;
            let volume = state.volume(pressure, phase)?;
            let enthalpy = state.enthalpy(volume)?;
            let entropy = state.entropy(volume)?;
            let computed = state.pressure(volume)?.value;
            let compressibility = computed * volume / (thermopack_rs::R * temperature);
            println!(
                "{:.17e},{:.17e},{:.17e},{:.17e}",
                volume, enthalpy, entropy, compressibility
            );
            continue;
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
        if matches!(
            property_config[0].trim_start_matches("V_"),
            "MIX_BP" | "MIX_DP" | "MIX_BT" | "MIX_DT"
        ) {
            let z: Vec<f64> = fields[4]
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?;
            let result = match property_config[0].trim_start_matches("V_") {
                "MIX_BP" => model.bubble_pressure(t, &z),
                "MIX_DP" => model.dew_pressure(t, &z),
                "MIX_BT" => model.bubble_temperature(t, &z),
                _ => model.dew_temperature(t, &z),
            };
            match result {
                Ok(s) => {
                    let mut values = vec![
                        s.temperature,
                        s.pressure,
                        s.liquid_volume,
                        s.vapor_volume,
                        s.dpdt,
                        s.fugacity_residual,
                        s.iterations as f64,
                    ];
                    values.extend(s.liquid);
                    values.extend(s.vapor);
                    println!(
                        "{}",
                        values
                            .iter()
                            .map(|v| format!("{v:.17e}"))
                            .collect::<Vec<_>>()
                            .join(",")
                    );
                }
                Err(e) => println!("ERROR {e}"),
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
