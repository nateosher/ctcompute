use crate::computation_target::ComputationTarget;
use crate::surv_sim_settings::SurvSimSettings;

fn tte_compute(surv_sim_settings: SurvSimSettings) {
    match surv_sim_settings.target {
        ComputationTarget::SampleSize => tte_compute_sample_size(surv_sim_settings),
        ComputationTarget::Alpha => tte_compute_alpha(surv_sim_settings),
        ComputationTarget::Beta => tte_compute_beta(surv_sim_settings),
        ComputationTarget::EffectSize => tte_compute_effect_size(surv_sim_settings),
    }
}

fn tte_compute_sample_size(surv_sim_settings: SurvSimSettings) {
    unimplemented!()
}

fn tte_compute_alpha(surv_sim_settings: SurvSimSettings) {
    unimplemented!()
}

fn tte_compute_beta(surv_sim_settings: SurvSimSettings) {
    unimplemented!()
}

fn tte_compute_effect_size(surv_sim_settings: SurvSimSettings) {
    unimplemented!()
}
