use crate::computation::types::ThetaComputeError;
use crate::integration::types::{NormalDistErr, TrialBoundsError};
use crate::spending::types::SpendingFcnErr;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum CtsimErr {
    #[error("while evaluating spending function: {0}")]
    SpendingFcn(SpendingFcnErr),
    #[error("while evaluating normal distribution: {0}")]
    NormalDist(NormalDistErr),
    #[error("while computing trial bounds: {0}")]
    TrialBounds(TrialBoundsError),
    #[error("while computing TTE characteristic: {0}")]
    ThetaCompute(ThetaComputeError),
}
