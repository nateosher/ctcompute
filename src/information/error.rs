//----------------------------------------
// information errors
//----------------------------------------
use crate::error::CtcomputeErr;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum TrialBoundsError {
    #[error("failed to converge (computed alpha: {0}, target alpha: {1}, tolerance: {2}")]
    FailedToConverge(f64, f64, f64),
}

impl Into<CtcomputeErr> for TrialBoundsError {
    fn into(self) -> CtcomputeErr {
        CtcomputeErr::TrialBounds(self)
    }
}

#[derive(Error, Debug)]
pub enum NormalDistErr {
    #[error("arguments to quantile function should be in [0, 1]; got {0}")]
    QuantileOutOfBounds(f64),
}

impl Into<CtcomputeErr> for NormalDistErr {
    fn into(self) -> CtcomputeErr {
        CtcomputeErr::NormalDist(self)
    }
}

#[derive(Error, Debug)]
pub enum RootFindErr {
    #[error("f(lower_bound) is larger than target; use smaller lower bound")]
    BadLowerBound,
}

impl Into<CtcomputeErr> for RootFindErr {
    fn into(self) -> CtcomputeErr {
        CtcomputeErr::RootFind(self)
    }
}
