use std::fmt;

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Error {
    UnknownComponent(String),
    InvalidInput(&'static str),
    Unsupported(&'static str),
    /// An unsupported label with a descriptive message.
    UnsupportedDetail(String),
    NoPhysicalRoot,
    NotConverged(&'static str),
}

pub type Result<T> = std::result::Result<T, Error>;

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownComponent(s) => write!(f, "unknown component: {s}"),
            Self::InvalidInput(s) => write!(f, "invalid input: {s}"),
            Self::Unsupported(s) => write!(f, "not implemented: {s}"),
            Self::UnsupportedDetail(s) => write!(f, "not implemented: {s}"),
            Self::NoPhysicalRoot => f.write_str("no physical volume root"),
            Self::NotConverged(s) => write!(f, "solver did not converge: {s}"),
        }
    }
}
impl std::error::Error for Error {}

pub(crate) fn positive(value: f64, name: &'static str) -> Result<()> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(Error::InvalidInput(name))
    }
}
