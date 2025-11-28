use std::f64;

use crate::computation::error::InformationComputeError;
use crate::error::CtsimErr;
use crate::integration::{
    integrate::{find_bounds, psi_k},
    std_normal::std_normal_quantile,
    types::IntegralType,
};
use crate::spending::{spending_fcns::compute_spending_vec, types::SpendingFcn};

pub fn compute_information(
    alpha: f64,
    target_power: f64,
    delta: f64,
    maybe_lower_spending_fcn_type: Option<SpendingFcn>,
    maybe_upper_spending_fcn_type: Option<SpendingFcn>,
    maybe_look_fractions: Option<&Vec<f64>>,
    tol: f64,
) -> Result<f64, CtsimErr> {
    // TODO: handle 1-look case separately
    if maybe_look_fractions.is_none() {
        let z_alpha = std_normal_quantile(1.0 - alpha)?;
        let z_beta = std_normal_quantile(target_power)?;
        let theta = z_alpha + z_beta;
        //    theta = delta * \sqrt(I(1))
        // => (theta / delta)^2 = I(1)
        #[allow(non_snake_case)]
        let target_I = (theta / delta).powf(2.0);
        return Ok(target_I);
    }

    let look_fractions = maybe_look_fractions.unwrap();
    let alpha_spend = compute_spending_vec(
        look_fractions,
        alpha,
        maybe_lower_spending_fcn_type,
        maybe_upper_spending_fcn_type,
    )?;
    let bounds = find_bounds(&alpha_spend, look_fractions, 32, tol)?;

    #[allow(non_snake_case)]
    let mut lower_I: f64 = 1.0;
    #[allow(non_snake_case)]
    let mut upper_I: f64 = 1000.0;
    #[allow(non_snake_case)]
    let mut cur_I: f64 = (lower_I + upper_I) / 2.0;

    let mut theta = delta * cur_I.sqrt();
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
        theta = delta * cur_I.sqrt();
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
        return Err(InformationComputeError::FailedToConverge.into());
    }

    #[allow(non_snake_case)]
    let minimal_I = match sufficient_total_informations.into_iter().reduce(f64::min) {
        Some(i) => i,
        None => Err(InformationComputeError::NoValueFound.into())?,
    };

    Ok(minimal_I)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tte_ss_compute() {
        let computed_information = compute_information(
            0.025,
            0.9,
            0.5_f64.ln(),
            Some(SpendingFcn::LDOF),
            None, // Some(SpendingFcn::LDOF),
            Some(&vec![0.7, 1.0]),
            0.0001,
        )
        .unwrap();
        assert_eq!((computed_information * 4.0).ceil(), 89.0);
    }

    #[test]
    fn single_look_info() {
        let computed_information = compute_information(
            0.025,
            0.9,
            0.5_f64.ln(),
            Some(SpendingFcn::LDOF),
            None, // Some(SpendingFcn::LDOF),
            None,
            0.0001,
        )
        .unwrap();
        assert_eq!((computed_information * 4.0).ceil(), 88.0);
    }
}
