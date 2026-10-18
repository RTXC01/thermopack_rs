use thermopack_rs::{Cubic, CubicKind, CubicOptions, TranslationKind, VolumeTranslation};

pub fn fixture_model(fields: &[&str]) -> Cubic {
    let kind = match fields[0] {
        "PR" => CubicKind::PengRobinson,
        "SRK" => CubicKind::SoaveRedlichKwong,
        "RK" => CubicKind::RedlichKwong,
        "VDW" => CubicKind::VanDerWaals,
        _ => panic!("model"),
    };
    let config: Vec<_> = fields[5].split('@').collect();
    let mut model = if config.contains(&"qhelper") {
        Cubic::quantum(&fields[1].split(',').collect::<Vec<_>>()).unwrap()
    } else {
        Cubic::with_options(
            kind,
            &fields[1].split(',').collect::<Vec<_>>(),
            CubicOptions {
                alpha: Some(thermopack_rs::AlphaCorrelation::from_name(config[1]).unwrap()),
                parameter_reference: config[2],
                volume_translation: config[3] == "1",
            },
        )
        .unwrap()
    };
    for operation in &config[4..] {
        let data: Vec<_> = operation.split(':').collect();
        match data[0] {
            "init" => model.init_peneloux_volume_translation(data[1]).unwrap(),
            "disable" => model.disable_volume_translation(),
            "bounds" => {} // Applied by the caloric fixture driver.
            "qhelper" => {}
            "rr" => {}      // Inputs to the standalone RR iteration fixture.
            "quality" => {} // Pure coexistence acoustic fixture input.
            "roots" => {}   // Optional acoustic phase-root overrides.
            "amount" => {}  // Extensive TV acoustic input.
            "beta" => {
                let params: Vec<f64> = data[3].split(',').map(|v| v.parse().unwrap()).collect();
                model
                    .set_beta(
                        data[1].parse().unwrap(),
                        thermopack_rs::Beta::new(
                            thermopack_rs::BetaCorrelation::from_name(data[2]).unwrap(),
                            params.try_into().unwrap(),
                        )
                        .unwrap(),
                    )
                    .unwrap();
            }
            "kij" => model
                .set_kij(
                    data[1].parse().unwrap(),
                    data[2].parse().unwrap(),
                    data[3].parse().unwrap(),
                )
                .unwrap(),
            "lij" => model
                .set_lij(
                    data[1].parse().unwrap(),
                    data[2].parse().unwrap(),
                    data[3].parse().unwrap(),
                )
                .unwrap(),
            "set" => {
                let coefficients = data[3]
                    .split(',')
                    .map(|v| v.parse().unwrap())
                    .collect::<Vec<f64>>();
                model
                    .set_volume_translation(
                        data[1].parse().unwrap(),
                        VolumeTranslation::new(
                            TranslationKind::from_code(data[2].parse().unwrap()).unwrap(),
                            coefficients.try_into().unwrap(),
                        )
                        .unwrap(),
                    )
                    .unwrap();
            }
            _ => panic!("operation"),
        }
    }
    model
}
