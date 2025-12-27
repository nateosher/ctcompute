#[derive(Debug)]
pub struct Trial {
    pub accrual_duration: f64,
    pub trial_duration: f64,
    pub n_events: usize,
    pub n_patients: usize,
    pub h0_expected_accrual_duration: f64,
    pub h1_expected_accrual_duration: f64,
    pub h0_expected_trial_duration: f64,
    pub h1_expected_trial_duration: f64,
    pub h0_expected_sample_size: f64,
    pub h1_expected_sample_size: f64,
}
