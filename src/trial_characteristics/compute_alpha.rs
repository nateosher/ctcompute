use std::f64;

use crate::duration::types::EnrollmentRate;
use crate::error::CtcomputeErr;
use crate::events::expected_events::expected_events_piecewise_arms;
use crate::information::std_normal::{std_normal_cdf, std_normal_quantile};

pub fn compute_alpha(
    power: f64,
    prop_treated: f64,
    lambda_event_trt: f64,
    lambda_event_ctrl: f64,
    maybe_lambda_dropout: Option<f64>,
    enrollment_rate: &EnrollmentRate,
    accrual_end_time: f64,
    followup_end_time: f64,
) -> Result<f64, CtcomputeErr> {
    let log_hr = (lambda_event_trt / lambda_event_ctrl).ln();
    let expected_events = expected_events_piecewise_arms(
        prop_treated,
        maybe_lambda_dropout,
        lambda_event_trt,
        lambda_event_ctrl,
        enrollment_rate,
        accrual_end_time,
        followup_end_time,
    );
    // Just formula in events_needed, solved for I
    let r = prop_treated / (1. - prop_treated);
    #[allow(non_snake_case)]
    let implied_I = (expected_events * (4. * r) / ((1. + r) * (1. + r))) / 4.;
    let z_beta = std_normal_quantile(1. - power)?;

    //    theta = z_alpha + z_beta
    //    theta = delta * \sqrt(I(1))
    // => z_alpha = delta*\sqrt(I(1)) - z_beta
    let z_one_minus_alpha = log_hr * implied_I.sqrt() - z_beta;

    Ok(std_normal_cdf(z_one_minus_alpha))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compute_alpha_1() {
        // See compute_trial test 1
        let er = EnrollmentRate::new(vec![0., 1., 2.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");

        let alpha = compute_alpha(
            0.9,               // Power
            0.5,               // prop treated
            0.018,             // hazard treat
            0.036,             // hazard control
            None,              // dropout
            &er,               // enrollment
            4. + 1. / 3.,      // accrual end time
            91.99448649088544, // followup end time
        )
        .expect("failed to compute alpha");

        assert!((alpha - 0.025) < 0.001);
        // Since we "round up" events needed, information needed,
        // etc., alpha should be slightly smaller than nominal
        // trial alpha
        assert!(alpha < 0.025);
    }

    #[test]
    fn compute_alpha_2() {
        // See compute_trial test 2
        let er = EnrollmentRate::new(vec![0., 10., 20.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");

        let alpha = compute_alpha(
            0.9,                // Power
            0.5,                // prop treated
            0.018,              // hazard treat
            0.036,              // hazard control
            None,               // dropout
            &er,                // enrollment
            23. + 1. / 3.,      // accrual end time
            24.459918657938637, // followup end time
        )
        .expect("failed to compute alpha");

        assert!((alpha - 0.025) < 0.001);
        // Since we "round up" events needed, information needed,
        // etc., alpha should be slightly smaller than nominal
        // trial alpha
        assert!(alpha < 0.025);
    }
}
