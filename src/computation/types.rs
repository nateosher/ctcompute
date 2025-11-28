//----------------------------------------
// computation mod types
//----------------------------------------
#[derive(Default, Debug, PartialEq, Clone, Copy)]
pub enum ComputationTarget {
    #[default]
    SampleSize,
    Alpha,
    Power,
    EffectSize,
}

pub struct EnrollmentSimResult {
    pub enrollment_times: Vec<f64>,
    pub expected_time: f64,
}
