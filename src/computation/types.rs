//----------------------------------------
// computation mod types
//----------------------------------------
use crate::error::CtsimErr;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum ThetaComputeError {
    #[error("failed to converge while computing theta")]
    FailedToConverge,
    #[error("no valid value found while computing theta")]
    NoValueFound,
}

impl Into<CtsimErr> for ThetaComputeError {
    fn into(self) -> CtsimErr {
        CtsimErr::ThetaCompute(self)
    }
}

#[derive(Default, Debug, PartialEq, Clone, Copy)]
pub enum ComputationTarget {
    #[default]
    SampleSize,
    Alpha,
    Beta,
    EffectSize,
}

pub struct EnrollmentSimResult {
    pub enrollment_times: Vec<f64>,
    pub expected_time: f64,
}
