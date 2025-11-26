use std::f64;

use crate::{
    computation_target::ComputationTarget,
    ctsim_err::CtsimErr,
    integrate::{IntegralType, find_bounds, psi_k},
    spending_fcns::{SpendingFcn, compute_spending_vec},
};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum TTEComputeError {
    #[error("failed to converge while computing TTE {0:?}")]
    FailedToConverge(ComputationTarget),
    #[error("no valid value found while computing TTE {0:?}")]
    NoValueFound(ComputationTarget),
}

impl Into<CtsimErr> for TTEComputeError {
    fn into(self) -> CtsimErr {
        CtsimErr::TTECompute(self)
    }
}

// pub fn tte_compute(n_per_arm: Option<usize>) {
//     match surv_sim_settings.target {
//         ComputationTarget::SampleSize => tte_compute_sample_size(surv_sim_settings),
//         ComputationTarget::Alpha => tte_compute_alpha(surv_sim_settings),
//         ComputationTarget::Beta => tte_compute_beta(surv_sim_settings),
//         ComputationTarget::EffectSize => tte_compute_effect_size(surv_sim_settings),
//     }
// }

pub fn tte_compute_sample_size(
    alpha: f64,
    target_power: f64,
    hr: f64,
    maybe_lower_spending_fcn_type: Option<SpendingFcn>,
    maybe_upper_spending_fcn_type: Option<SpendingFcn>,
    look_fractions: &Vec<f64>,
    tol: f64,
) -> Result<usize, CtsimErr> {
    let log_hr = hr.ln();
    // TODO: handle 1-look case separately
    let alpha_spend = compute_spending_vec(
        look_fractions,
        alpha,
        maybe_lower_spending_fcn_type,
        maybe_upper_spending_fcn_type,
    )?;
    let bounds = find_bounds(&alpha_spend, look_fractions, 32, tol)?;

    let mut lower_I: f64 = 1.0;
    let mut upper_I: f64 = 1000.0;
    let mut cur_I: f64 = (lower_I + upper_I) / 2.0;
    let mut theta = log_hr * cur_I.sqrt();
    let mut cur_power_lower: f64 = psi_k(&bounds, look_fractions, theta, IntegralType::Lower, 32)
        .iter()
        .sum();
    let mut cur_power_upper: f64 = psi_k(&bounds, look_fractions, theta, IntegralType::Upper, 32)
        .iter()
        .sum();
    let mut cur_power = cur_power_lower + cur_power_upper;
    let mut sufficient_total_informations: Vec<f64> = vec![];
    if cur_power > target_power {
        sufficient_total_informations.push(cur_I);
    }

    while (cur_power - target_power).abs() > tol && (lower_I - upper_I).abs() > tol {
        if cur_power > target_power {
            upper_I = cur_I;
        } else {
            lower_I = cur_I;
        }
        cur_I = (lower_I + upper_I) / 2.0;
        theta = log_hr * cur_I.sqrt();
        cur_power_lower = psi_k(&bounds, look_fractions, theta, IntegralType::Lower, 32)
            .iter()
            .sum();
        cur_power_upper = psi_k(&bounds, look_fractions, theta, IntegralType::Upper, 32)
            .iter()
            .sum();
        cur_power = cur_power_lower + cur_power_upper;

        if cur_power > target_power {
            sufficient_total_informations.push(cur_I);
        }
    }

    if (cur_power - target_power).abs() > tol {
        return Err(TTEComputeError::FailedToConverge(ComputationTarget::SampleSize).into());
    }

    let minimal_I = match sufficient_total_informations.into_iter().reduce(f64::min) {
        Some(i) => i,
        None => Err(TTEComputeError::NoValueFound(ComputationTarget::SampleSize).into())?,
    };

    Ok((minimal_I * 4.0).ceil() as usize)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tte_ss_compute() {
        let computed_ss = tte_compute_sample_size(
            0.025,
            0.9,
            0.5,
            Some(SpendingFcn::LDOF),
            Some(SpendingFcn::LDOF),
            &vec![0.7, 1.0],
            0.0001,
        );
        assert_eq!(computed_ss.unwrap(), 89);
    }
}
