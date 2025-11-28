//----------------------------------------
// computation errors
//----------------------------------------
use crate::error::CtsimErr;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum InformationComputeError {
    #[error("failed to converge while computing total information")]
    FailedToConverge,
    #[error("no valid value found while computing total information")]
    NoValueFound,
}

impl Into<CtsimErr> for InformationComputeError {
    fn into(self) -> CtsimErr {
        CtsimErr::InformationCompute(self)
    }
}
