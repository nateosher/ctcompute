//----------------------------------------
// Enrollment errors
//----------------------------------------

use crate::error::CtcomputeErr;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum EnrollmentComputationError {
    #[error(
        "lengths of enrollment times and rates don't match (times length {time_length}, \
        rates length {rate_length}"
    )]
    TimeRateLengths {
        time_length: usize,
        rate_length: usize,
    },
}

impl Into<CtcomputeErr> for EnrollmentComputationError {
    fn into(self) -> CtcomputeErr {
        CtcomputeErr::EnrollmentCompute(self)
    }
}
