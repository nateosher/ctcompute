use std::f64;

use crate::duration::min_followup::min_followup;
use crate::duration::{expected_enrollment::expected_enrollment, types::EnrollmentRate};
use crate::error::CtcomputeErr;
use crate::events::events_needed::events_needed;
use crate::events::expected_events::expected_events_piecewise_arms;
use crate::spending::types::SpendingFcn;
use crate::util::root_find::root_find_monotonic;

/// Computes sample size range given alpha, power, spending function,
/// look fractions, and probability of randomization to treatment
/// For now assumes TTE endpoint
/// Assumes fixed enrollment rates, where last rate is extended as long as necessary
/// Given this, computes range of possible accrual durations (and thus sample sizes and
/// total trial durations)
/// delta indicates how fine the grid of accrual times to check should be
/// perc_change_stop indicates the percentage reduction of followup time
/// which should be considered "diminishing returns." In other words,
/// perc_change_stop = 0.01 means that when an additional increment of delta
/// results in a percent change in total necessary followup time of 1% or less,
/// we consider the corresponding accrual time to be the largest one worth
/// considering.
pub fn compute_ss_range(
    alpha: f64,
    power: f64,
    maybe_lower_spending_fcn_type: Option<SpendingFcn>,
    maybe_upper_spending_fcn_type: Option<SpendingFcn>,
    maybe_look_fractions: Option<&Vec<f64>>,
    prop_treated: f64,
    lambda_event_trt: f64,
    lambda_event_ctrl: f64,
    maybe_lambda_dropout: Option<f64>,
    enrollment_rate: &EnrollmentRate,
    tol: f64,
    delta: f64,
    perc_change_stop: f64,
) -> Result<(usize, usize), CtcomputeErr> {
    let n_events_needed = events_needed(
        alpha,
        power,
        maybe_lower_spending_fcn_type,
        maybe_upper_spending_fcn_type,
        maybe_look_fractions,
        prop_treated,
        lambda_event_trt,
        lambda_event_ctrl,
        tol,
    )?;

    //----------------------------------------
    // Finding smallest amount of subjects
    //----------------------------------------
    // Assuming infinite followup, i.e. observe all events barring dropout
    // This would be the absolute minimum # patients needed
    // Right now this is redundant, since assumign constant dropout hazard and piecewise
    // enrollment it can be computed directly but this is the most general form, and
    // will allow the incorporation of time varying dropout rates

    // Treat as root find, since expected events is monotonic in accrual duration (x)
    let f_expected_by_accrual = |x| {
        expected_events_piecewise_arms(
            prop_treated,
            maybe_lambda_dropout,
            lambda_event_trt,
            lambda_event_ctrl,
            enrollment_rate,
            x, // x is accrual duration
            f64::INFINITY,
        )
    };

    let min_accrual_dur = root_find_monotonic(f_expected_by_accrual, 0., n_events_needed, tol)?;

    //----------------------------------------
    // Finding largest amount of subjects
    //----------------------------------------
    // Obviously you could enroll an arbitrarily large number of subjects; the goal
    // here is to find the point of diminishing returns
    // The added difficulty is that we also need to find the trial duration for each
    // corresponding enrollment duration

    // Convenience closure to compute minimum necessary followup time according to
    // accrual time with everything else fixed, since this is the slope we are
    // interested in
    let min_followup_by_accrual = |accrual_time: f64| {
        min_followup(
            n_events_needed,
            prop_treated,
            lambda_event_trt,
            lambda_event_ctrl,
            maybe_lambda_dropout,
            enrollment_rate,
            accrual_time,
            tol,
        )
    };

    let mut cur_accrual_time = min_accrual_dur * 1.01;
    let mut accrual_times = vec![cur_accrual_time];
    accrual_times.reserve(1_000);
    let mut followup_times = vec![min_followup_by_accrual(cur_accrual_time)?];
    followup_times.reserve(1_000);
    let mut perc_change = f64::INFINITY;

    while perc_change > perc_change_stop {
        cur_accrual_time += delta;
        accrual_times.push(cur_accrual_time);
        followup_times.push(min_followup_by_accrual(cur_accrual_time)?);

        let cur_followup_time = followup_times[followup_times.len() - 1];
        let prev_followup_time = followup_times[followup_times.len() - 2];
        perc_change = (prev_followup_time - cur_followup_time) / prev_followup_time;
    }

    //----------------------------------------
    // Compute final enrollment + return
    //----------------------------------------

    let max_accrual_dur = accrual_times[accrual_times.len() - 1];
    let min_sample_size = expected_enrollment(min_accrual_dur, &enrollment_rate)
        .unwrap()
        .round();
    let max_sample_size = expected_enrollment(max_accrual_dur, &enrollment_rate)
        .unwrap()
        .round();

    Ok((min_sample_size as usize, max_sample_size as usize))
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn ss_range_comparison() {
        let er = EnrollmentRate::new(vec![0., 1., 2.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");
        let range_1 = compute_ss_range(
            0.025,  // alpha
            0.9,    // power
            None,   // lower spending function
            None,   // upper spending function
            None,   // look fractions
            0.5,    // prop treated
            0.018,  // hazard for treated group
            0.036,  // hazard for control group
            None,   // hazard of dropout
            &er,    // enrollment rate
            0.001,  // tol
            0.1,    // delta
            0.0001, // min_perc_change
        )
        .expect("failed to compute first sample size range");

        // Based on the range EAST gives
        assert_eq!(range_1.0, 88);
        assert_eq!(range_1.1, 473);

        // Larger minimum percent change means we should consider smaller
        // maximum accrual durations; however, minimum accrual durations
        // should be the same
        let range_2 = compute_ss_range(
            0.025, // alpha
            0.9,   // power
            None,  // lower spending function
            None,  // upper spending function
            None,  // look fractions
            0.5,   // prop treated
            0.018, // hazard for treated group
            0.036, // hazard for control group
            None,  // hazard of dropout
            &er,   // enrollment rate
            0.001, // tol
            0.1,   // delta
            0.05,  // min_perc_change
        )
        .expect("failed to compute second sample size range");

        // Lower accrual/sample sizes should be the same
        assert_eq!(range_1.0, range_2.0);

        // However, maximum sample size of second range should be smaller, since
        // we stop at steeper dropoff (5% vs. 0.01%)
        assert!(range_1.1 > range_2.1);
    }

    #[test]
    fn ss_range_comparison_2() {
        let er = EnrollmentRate::new(vec![0., 5., 10.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");
        let range = compute_ss_range(
            0.025,                     // alpha
            0.9,                       // power
            Some(SpendingFcn::LDOF),   // lower spending function
            None,                      // upper spending function
            Some(&vec![0.5, 0.7, 1.]), // look fractions
            2. / 3.,                   // prop treated
            0.04,                      // hazard for treated group
            0.06,                      // hazard for control group
            Some(0.00878),             // hazard of dropout
            &er,                       // enrollment rate
            0.00000001,                // tol
            0.1,                       // delta
            0.0001,                    // min_perc_change
        )
        .expect("failed to compute first sample size range");

        assert_eq!(range.0, 349);
        // Again, no set upper bound, but ideally want it to be within
        // ~ 10 patients or so of EAST
        assert!((range.1 as i64 - 752_i64).abs() < 10);
    }

    #[test]
    fn ss_range_comparison_3() {
        let er = EnrollmentRate::new(vec![0., 10., 15.], vec![5., 7., 11.])
            .expect("failed to construct enrollment rate object");
        let range = compute_ss_range(
            0.025,                     // alpha
            0.9,                       // power
            Some(SpendingFcn::LDOF),   // lower spending function
            None,                      // upper spending function
            Some(&vec![0.3, 0.8, 1.]), // look fractions
            0.5,                       // prop treated
            0.3,                       // hazard for treated group
            0.4,                       // hazard for control group
            Some(0.00878),             // hazard of dropout
            &er,                       // enrollment rate
            0.00001,                   // tol
            0.1,                       // delta
            0.0001,                    // min_perc_change
        )
        .expect("failed to compute first sample size range");

        assert!((range.0 as i64 - 532).abs() <= 1);
        // Again, no set upper bound, but ideally want it to be within
        // ~ 10 patients or so of EAST
        assert!((range.1 as i64 - 561_i64).abs() < 10);
    }
}
