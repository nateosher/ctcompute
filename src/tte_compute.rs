use std::f64;

use crate::{
    computation_target::ComputationTarget,
    ctsim_err::CtsimErr,
    integrate::{find_bounds, psi_k},
    spending_fcns::{SpendingFcn, lan_demets_obrien_fleming_vec},
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
    spending_fcn: SpendingFcn,
    look_fractions: &Vec<f64>,
    tol: f64,
) -> Result<usize, CtsimErr> {
    let log_hr = hr.ln();
    // TODO: handle 1-look case separately
    let alpha_spend = match spending_fcn {
        SpendingFcn::LDOF => lan_demets_obrien_fleming_vec(look_fractions, alpha),
    }?;
    let bounds = find_bounds(&alpha_spend, look_fractions, 32, tol)?;

    let mut lower_I: f64 = 1.0;
    let mut upper_I: f64 = 1000.0;
    let mut cur_I: f64 = (lower_I + upper_I) / 2.0;
    let mut theta = log_hr * cur_I.sqrt();
    let mut cur_power: f64 = psi_k(&bounds, look_fractions, theta, 32).iter().sum();
    let mut sufficient_total_informations: Vec<f64> = vec![];
    if cur_power > target_power {
        sufficient_total_informations.push(cur_I);
    }
    println!("cur power: {cur_power}");
    println!("cur I: {cur_I}");
    while (cur_power - target_power).abs() > tol && (lower_I - upper_I).abs() > tol {
        if cur_power > target_power {
            upper_I = cur_I;
        } else {
            lower_I = cur_I;
        }
        cur_I = (lower_I + upper_I) / 2.0;
        theta = log_hr * cur_I.sqrt();
        cur_power = psi_k(&bounds, look_fractions, theta, 32).iter().sum();
        if cur_power > target_power {
            sufficient_total_informations.push(cur_I);
        }

        println!("cur power: {cur_power}");
        println!("cur I: {cur_I}");
    }

    if (cur_power - target_power).abs() > tol {
        return Err(TTEComputeError::FailedToConverge(ComputationTarget::SampleSize).into());
    }

    println!("{sufficient_total_informations:?}");
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
        let computed_ss =
            tte_compute_sample_size(0.05, 0.9, 0.5, SpendingFcn::LDOF, &vec![0.7, 1.0], 0.0001);
        println!("{computed_ss:?}");
    }
}
