use crate::duration::error::EnrollmentComputationError;
//----------------------------------------
// Crate error type
//----------------------------------------
pub use crate::information::error::*;
pub use crate::sample_size::error::*;
pub use crate::spending::error::*;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum CtcomputeErr {
    #[error("while evaluating spending function: {0}")]
    SpendingFcn(SpendingFcnErr),
    #[error("while evaluating normal distribution: {0}")]
    NormalDist(NormalDistErr),
    #[error("while computing trial bounds: {0}")]
    TrialBounds(TrialBoundsError),
    #[error("while computing TTE characteristic: {0}")]
    InformationCompute(InformationComputeError),
    #[error("while computing enrollment: {0}")]
    EnrollmentCompute(EnrollmentComputationError),
    #[error("while computing sample size: {0}")]
    SampleSizeCompute(SampleSizeComputeError),
    #[error("while performing root find: {0}")]
    RootFind(RootFindErr),
}
