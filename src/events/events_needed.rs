use std::f64;

use crate::error::CtcomputeErr;
use crate::information::compute_information::compute_information;
use crate::spending::types::SpendingFcn;

pub fn events_needed(
    alpha: f64,
    power: f64,
    maybe_lower_spending_fcn_type: Option<&SpendingFcn>,
    maybe_upper_spending_fcn_type: Option<&SpendingFcn>,
    maybe_look_fractions: Option<&Vec<f64>>,
    prop_treated: f64,
    lambda_event_trt: f64,
    lambda_event_ctrl: f64,
    tol: f64,
) -> Result<f64, CtcomputeErr> {
    //----------------------------------------
    // Amount of required information
    //----------------------------------------
    #[allow(non_snake_case)]
    let I = compute_information(
        alpha,
        power,
        (lambda_event_trt / lambda_event_ctrl).ln(),
        maybe_lower_spending_fcn_type,
        maybe_upper_spending_fcn_type,
        maybe_look_fractions,
        tol,
    )?;

    //----------------------------------------
    // # of events needed to achieve desired information fraction
    //----------------------------------------
    let r = prop_treated / (1. - prop_treated);
    // Base # events is I / 0.25 = 4I under equal allocation + null hypothesis
    // TODO: add option to compute under variance under alternative hypothesis
    // Inflate I according to randomization ratio r:
    // if p = probability of being in treatment group, r = p / (1 - p)
    let n_events_needed = ((I * 4.) * (1. + r) * (1. + r) / (4. * r)).ceil();

    Ok(n_events_needed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::duration::{expected_enrollment::expected_enrollment, types::EnrollmentRate};
    use crate::events::expected_events::expected_events_piecewise_arms;
    use crate::util::root_find::root_find_monotonic;

    #[test]
    fn events_needed_1() {
        let er = EnrollmentRate::new(vec![0., 5., 10.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");

        let n_needed = events_needed(
            0.025,
            0.9,
            Some(&SpendingFcn::LDOF),
            None,
            Some(&vec![0.5, 0.7, 1.0]),
            2. / 3.,
            0.04,
            0.06,
            0.000001,
        )
        .expect("could not compute necessary number of events");

        let f_expected_by_accrual = |x| {
            expected_events_piecewise_arms(
                2. / 3.,
                Some(0.00878),
                0.04,
                0.06,
                &er,
                x, // x is accrual duration
                f64::INFINITY,
            )
        };

        let min_accrual_dur = root_find_monotonic(f_expected_by_accrual, 0., n_needed, 0.0001)
            .expect("failed to compute minimum accrual duration");
        let min_sample_size = expected_enrollment(min_accrual_dur, &er).unwrap().round();

        assert_eq!(min_sample_size, 349.);
        // For some reason, these values differ slightly, regardless of tol I set
        assert!((min_accrual_dur - 16.663).abs() < 0.1);
    }
}
