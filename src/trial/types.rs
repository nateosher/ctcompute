#[derive(Debug)]
pub struct Trial {
    pub accrual_duration: f64,
    pub trial_duration: f64,
    pub n_events: usize,
    pub n_patients: usize,
    pub maybe_h0_expected_accrual_duration: Option<f64>, // TODO: these actually
    pub maybe_h1_expected_accrual_duration: Option<f64>, // are not optional, since
    pub maybe_h0_expected_trial_duration: Option<f64>,   // trial length can differ
    pub maybe_h1_expected_trial_duration: Option<f64>,   // even with one look
    pub maybe_h0_expected_sample_size: Option<f64>,
    pub maybe_h1_expected_sample_size: Option<f64>,
}
