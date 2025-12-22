#[derive(Debug, Copy, Clone)]
pub struct SampleSizeCalculation {
    pub min_accrual_dur: f64,
    pub max_accrual_dur: f64,
    pub min_followup_dur: f64,
    pub max_followup_dur: f64,
    pub min_sample_size: f64,
    pub max_sample_size: f64,
}
