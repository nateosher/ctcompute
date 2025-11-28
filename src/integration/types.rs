//----------------------------------------
// integration mod types
//----------------------------------------

use crate::error::CtsimErr;
use thiserror::Error;

#[derive(PartialEq, Clone, Copy)]
pub enum IntegralType {
    Upper,
    Lower,
}

#[derive(Error, Debug)]
pub enum TrialBoundsError {
    #[error("failed to converge (computed alpha: {0}, target alpha: {1}, tolerance: {2}")]
    FailedToConverge(f64, f64, f64),
}

impl Into<CtsimErr> for TrialBoundsError {
    fn into(self) -> CtsimErr {
        CtsimErr::TrialBounds(self)
    }
}

#[derive(Error, Debug)]
pub enum NormalDistErr {
    #[error("arguments to quantile function should be in [0, 1]; got {0}")]
    QuantileOutOfBounds(f64),
}

impl Into<CtsimErr> for NormalDistErr {
    fn into(self) -> CtsimErr {
        CtsimErr::NormalDist(self)
    }
}
