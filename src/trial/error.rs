use crate::error::CtcomputeErr;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum TrialComputeError {
    #[error("sample size ({n_patients}) is smaller than required # of events ({n_events_needed})")]
    InsufficientSampleSize {
        n_patients: usize,
        n_events_needed: usize,
    },
}

impl Into<CtcomputeErr> for TrialComputeError {
    fn into(self) -> CtcomputeErr {
        CtcomputeErr::TrialCompute(self)
    }
}
