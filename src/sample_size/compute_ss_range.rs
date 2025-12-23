use std::f64;

use crate::duration::min_followup::min_followup;
use crate::duration::{expected_enrollment::expected_enrollment, types::EnrollmentRate};
use crate::error::CtcomputeErr;
use crate::events::expected_events::expected_events_piecewise_arms;
use crate::util::root_find::root_find_monotonic;

/// Computes sample size given total necessary information
/// and probability of randomization to treatment
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
#[allow(non_snake_case)]
pub fn compute_ss_range(
    I: f64,
    prop_treated: f64,
    lambda_event_trt: f64,
    lambda_event_ctrl: f64,
    maybe_lambda_dropout: Option<f64>,
    enrollment_rate: &EnrollmentRate,
    tol: f64,
    delta: f64,
    perc_change_stop: f64,
) -> Result<(usize, usize), CtcomputeErr> {
    //----------------------------------------
    // # of events needed to achieve desired information fraction
    //----------------------------------------
    let r = prop_treated / (1. - prop_treated);
    // Base # events is I / 0.25 = 4I under equal allocation + null hypothesis
    // TODO: add option to compute under variance under alternative hypothesis
    // Inflate I according to randomization ratio r:
    // if p = probability of being in treatment group, r = p / (1 - p)
    let n_events_needed = ((I * 4.) * (1. + r) * (1. + r) / (4. * r)).ceil();

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

    let min_accrual_dur = accrual_times[0];
    let max_accrual_dur = accrual_times[accrual_times.len() - 1];
    let min_sample_size = expected_enrollment(min_accrual_dur, &enrollment_rate)
        .unwrap()
        .ceil();
    let max_sample_size = expected_enrollment(max_accrual_dur, &enrollment_rate)
        .unwrap()
        .ceil();

    Ok((min_sample_size as usize, max_sample_size as usize))
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn ss_range_comparison() {
        let enrollment_rate = EnrollmentRate::new(vec![0., 1., 2.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");
        let range_1 = compute_ss_range(
            22.,
            0.5,
            0.018,
            0.036,
            None,
            &enrollment_rate,
            0.001, // tol
            0.1,   // delta
            0.01,  // min_perc_change
        )
        .expect("failed to compute first sample size range");

        // Larger minimum percent change means we should consider smaller
        // maximum accrual durations; however, minimum accrual durations
        // should be the same
        let range_2 = compute_ss_range(
            22.,
            0.5,
            0.018,
            0.036,
            None,
            &enrollment_rate,
            0.001, // tol
            0.1,   // delta
            0.05,  // min_perc_change
        )
        .expect("failed to compute second sample size range");

        // Lower accrual/sample sizes should be the same
        assert_eq!(range_1.0, range_2.0);

        // However, maximum sample size of second range should be smaller, since
        // we stop at steeper dropoff (5% vs. 1%)
        assert!(range_1.1 > range_2.1);
    }

    #[test]
    fn ss_range_no_dropout_2() {
        let enrollment_rate = EnrollmentRate::new(vec![0., 5., 10.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");
        let range = compute_ss_range(
            64.,
            0.5,
            0.04,
            0.06,
            Some(0.00878),
            &enrollment_rate,
            0.001,  // tol
            0.1,    // delta
            0.0001, // min_perc_change
        );
        println!("{range:?}");
    }
}
