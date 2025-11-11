use crate::computation_target::ComputationTarget;
use crate::hypothesis_type::HypothesisType;

#[derive(Debug, Clone, Copy)]
pub struct SurvSimSettings {
    pub seed: u64,
    pub sample_size_trt: usize,
    pub sample_size_ctrl: usize,
    pub lambda_trt: f64,
    pub lambda_ctrl: f64,
    pub simulate_dropout: bool,
    pub lambda_trt_dropout: f64,
    pub lambda_ctrl_dropout: f64,
    pub left_ss_bound: i64,
    pub right_ss_bound: i64,
    pub left_es_bound: f64,
    pub es_tol: f64,
    pub es_ensure_tol: bool,
    pub right_es_bound: f64,
    pub alpha: f64,
    pub beta: f64,
    pub hypothesis: HypothesisType,
    pub reps_per_sim: usize,
    pub target: ComputationTarget,
    pub diff_arm_sizes: bool,
}

impl Default for SurvSimSettings {
    fn default() -> Self {
        Self {
            seed: 24601,
            sample_size_trt: 50,
            sample_size_ctrl: 50,
            lambda_trt: 0.5,
            lambda_ctrl: 0.5,
            simulate_dropout: false,
            lambda_trt_dropout: 0.000000000001,
            lambda_ctrl_dropout: 0.000000000001,
            left_ss_bound: 2,
            right_ss_bound: 100,
            left_es_bound: 0.0,
            right_es_bound: 0.0,
            es_tol: 0.001,
            es_ensure_tol: false,
            alpha: 0.05,
            beta: 0.1,
            hypothesis: HypothesisType::NotEqual,
            reps_per_sim: 10000,
            target: ComputationTarget::SampleSize,
            diff_arm_sizes: false,
        }
    }
}
