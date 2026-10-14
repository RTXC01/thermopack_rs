use thermopack_rs::{Cubic, CubicKind, Phase};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let model = Cubic::new(CubicKind::PengRobinson, &["CO2", "C1"])?;
    let state = model.prepare(300.0, &[0.8, 0.2])?;
    let volume = state.volume(5e6, Phase::Vapor)?;
    let mut phi = [0.0; 2];
    state.ln_fugacity_into(volume, &mut phi)?;
    println!("v = {volume:.12e} m3/mol");
    println!("{:?}", state.residual(volume)?);
    println!("ln(phi) = {phi:?}");
    Ok(())
}
