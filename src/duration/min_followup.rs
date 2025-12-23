use std::f64;

use crate::duration::types::EnrollmentRate;
use crate::error::{CtcomputeErr, RootFindErr};
use crate::events::expected_events::expected_events_piecewise_arms;
use crate::util::root_find::root_find_monotonic;

/// Given everything except followup time, computes minimum followup time required
/// to observe a given number of events
pub fn min_followup(
    n_events_needed: f64,
    prop_treated: f64,
    lambda_event_trt: f64,
    lambda_event_ctrl: f64,
    maybe_lambda_dropout: Option<f64>,
    enrollment_rate: &EnrollmentRate,
    accrual_end_time: f64,
    tol: f64,
) -> Result<f64, CtcomputeErr> {
    let f_expected_by_followup = |followup_dur| {
        expected_events_piecewise_arms(
            prop_treated,
            maybe_lambda_dropout,
            lambda_event_trt,
            lambda_event_ctrl,
            enrollment_rate,
            accrual_end_time,
            followup_dur,
        )
    };

    let min_followup_res = root_find_monotonic(
        f_expected_by_followup,
        accrual_end_time,
        n_events_needed,
        tol,
    );

    if min_followup_res.is_ok() {
        return min_followup_res;
    }
    // If the root find went wrong, then the lower bound was wrong, which means the enrollment
    // time is longer than it takes to observe the events
    match min_followup_res {
        Err(CtcomputeErr::RootFind(RootFindErr::BadLowerBound)) => true,
        _ => unreachable!(),
    };
    let f_expected_by_followup_and_accrual = |followup_dur| {
        expected_events_piecewise_arms(
            prop_treated,
            maybe_lambda_dropout,
            lambda_event_trt,
            lambda_event_ctrl,
            enrollment_rate,
            followup_dur, // accrual time is set to followup time as well
            followup_dur,
        )
    };

    root_find_monotonic(
        f_expected_by_followup_and_accrual,
        0.,
        n_events_needed,
        0.001,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn min_dur_1() {
        let enrollment_rate = EnrollmentRate::new(vec![0., 1., 2.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");
        let min_followup = min_followup(
            88.,
            0.5,
            0.018,
            0.036,
            None,
            &enrollment_rate,
            10.367,
            0.0001,
        )
        .expect("failed to compute minimum followup time");
        assert!((min_followup - 19.974).abs() < 0.001);
    }

    #[test]
    fn min_dur_2() {
        let enrollment_rate = EnrollmentRate::new(vec![0., 1., 2.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");
        let min_followup = min_followup(
            99.,
            2. / 3.,
            0.018,
            0.036,
            None,
            &enrollment_rate,
            11.533,
            0.0001,
        )
        .expect("failed to compute minimum followup time");
        assert!((min_followup - 22.39).abs() < 0.001);
    }

    #[test]
    fn min_dur_3() {
        let enrollment_rate = EnrollmentRate::new(vec![0., 5., 10.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");
        let min_followup = min_followup(
            88.,
            0.5,
            0.018,
            0.036,
            Some(0.00878),
            &enrollment_rate,
            14.9,
            0.0001,
        )
        .expect("failed to compute minimum followup time");

        assert!((min_followup - 23.652).abs() < 0.001);
    }

    #[test]
    fn min_dur_4() {
        let enrollment_rate = EnrollmentRate::new(vec![0., 1., 2.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");
        let min_followup = min_followup(
            88.,
            0.5,
            0.018,
            0.036,
            None,
            &enrollment_rate,
            4. + 1. / 3.,
            0.0001,
        )
        .expect("failed to compute minimum followup time");
        assert!((min_followup - 91.994).abs() < 0.001);
    }

    #[test]
    fn min_dur_excessive_enrollment() {
        let enrollment_rate = EnrollmentRate::new(vec![0.], vec![8.])
            .expect("failed to construct enrollment rate object");
        let min_followup =
            min_followup(88., 0.5, 0.018, 0.036, None, &enrollment_rate, 500., 0.0001)
                .expect("failed to compute minimum followup time");

        assert!((min_followup - 33.168).abs() < 0.001);
    }

    #[test]
    fn problem_case() {
        let er = EnrollmentRate::new(vec![0., 10., 20.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");

        let min_followup = min_followup(91., 0.5, 0.018, 0.036, None, &er, 10., 0.0001);
        println!("min_followup: {min_followup:?}");
    }
}
