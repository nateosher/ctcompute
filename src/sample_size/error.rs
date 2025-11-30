use crate::error::CtcomputeErr;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum InformationComputeError {
    #[error("failed to converge while computing total information")]
    FailedToConverge,
    #[error("no valid value found while computing total information")]
    NoValueFound,
}

impl Into<CtcomputeErr> for InformationComputeError {
    fn into(self) -> CtcomputeErr {
        CtcomputeErr::InformationCompute(self)
    }
}
