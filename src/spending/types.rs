//----------------------------------------
// spending mod types
//----------------------------------------
use crate::error::CtsimErr;
use thiserror::Error;

pub enum SpendingFcn {
    LDOF, // Lan-Demets Obrien Fleming
}

#[derive(Error, Debug)]
pub enum SpendingFcnErr {
    #[error("arguments to spending function should be in [0, 1]; got {0}")]
    OutOfBounds(f64),
    #[error("spending function vector should end with 1.0; got {0}")]
    BadLastSpend(f64),
    #[error("time vector was empty")]
    TimeVectorEmpty,
    #[error("total alpha spent should be in (0, 1); got {0}")]
    BadAlpha(f64),
    #[error("must specify at least one spending function")]
    NoSpendingFunctionSpecified,
}

impl Into<CtsimErr> for SpendingFcnErr {
    fn into(self) -> CtsimErr {
        CtsimErr::SpendingFcn(self)
    }
}

#[derive(Debug)]
pub enum AlphaSpendingValues {
    OneSidedUpper(Vec<f64>),
    OneSidedLower(Vec<f64>),
    TwoSided((Vec<f64>, Vec<f64>)),
}
