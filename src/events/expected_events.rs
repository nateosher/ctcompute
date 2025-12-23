use crate::duration::types::EnrollmentRate;

pub fn expected_events_piecewise_arms(
    prop_treated: f64,
    maybe_lambda_dropout: Option<f64>,
    lambda_event_trt: f64,
    lambda_event_ctrl: f64,
    enrollment_rate: &EnrollmentRate,
    accrual_end_time: f64,
    followup_end_time: f64,
) -> f64 {
    prop_treated
        * expected_events_piecewise(
            lambda_event_trt,
            maybe_lambda_dropout,
            enrollment_rate,
            accrual_end_time,
            followup_end_time,
        )
        + (1. - prop_treated)
            * expected_events_piecewise(
                lambda_event_ctrl,
                maybe_lambda_dropout,
                enrollment_rate,
                accrual_end_time,
                followup_end_time,
            )
}

fn expected_events_piecewise(
    lambda_event: f64,
    maybe_lambda_dropout: Option<f64>,
    enrollment_rate: &EnrollmentRate,
    accrual_end_time: f64,
    followup_end_time: f64,
) -> f64 {
    let EnrollmentRate::PiecewiseConstant {
        enrollment_durations,
        enrollment_rates,
    } = enrollment_rate;

    // First element is duration, second is total expected events
    let (_, expected_events) = enrollment_durations
        .into_iter()
        .zip(enrollment_rates.into_iter())
        .fold(
            (0., 0.),
            |(cumulative_dur, cumulative_expected), (cur_dur, cur_rate)| {
                let new_events = expected_events_uniform(
                    lambda_event,
                    maybe_lambda_dropout,
                    *cur_rate,
                    cumulative_dur,
                    (cumulative_dur + cur_dur).min(accrual_end_time),
                    followup_end_time,
                );

                (cumulative_dur + cur_dur, cumulative_expected + new_events)
            },
        );

    expected_events
}

/// Probability of observing an event, assuming uniform accrual
fn expected_events_uniform(
    lambda_event: f64,
    maybe_lambda_dropout: Option<f64>,
    accrual_rate: f64,
    accrual_start_time: f64,
    accrual_end_time: f64,
    followup_end_time: f64,
) -> f64 {
    let lambda_dropout = match maybe_lambda_dropout {
        Some(l) => l,
        None => 0.0,
    };
    let lambda_total = lambda_event + lambda_dropout;

    let num_1 = (-lambda_total * (followup_end_time - accrual_end_time)).exp();
    let num_2 = (-lambda_total * (followup_end_time - accrual_start_time)).exp();
    let part_1 = (accrual_end_time - accrual_start_time) - (num_1 - num_2) / lambda_total;
    let part_2 = accrual_rate * lambda_event / lambda_total;

    part_1 * part_2
}

#[cfg(test)]
mod tests {

    use std::f64;

    use super::*;

    #[test]
    fn infinite_followup_no_drop_1() {
        // If there's infinite followup, should just recruit rate x recruitment time
        let expected_events = expected_events_uniform(0.001, None, 8., 0., 10., f64::INFINITY);
        assert_eq!(expected_events, 80.);
    }

    #[test]
    fn infinite_followup_no_drop_2() {
        // Same deal, but showing that event rate doesn't matter when followup is infinite
        let expected_events_1 = expected_events_uniform(0.001, None, 8., 0., 10., f64::INFINITY);
        let expected_events_2 = expected_events_uniform(10., None, 8., 0., 10., f64::INFINITY);

        assert_eq!(expected_events_1, expected_events_2);
    }

    #[test]
    fn infinite_followup_with_drop_1() {
        // With infinite followup, still lose to dropout in proportion to event rate
        let expected_events = expected_events_uniform(0.1, Some(0.1), 8., 0., 10., f64::INFINITY);

        assert_eq!(expected_events, 40.);
    }

    #[test]
    fn infinite_followup_with_drop_2() {
        // With infinite followup, still lose to dropout in proportion to event rate
        let expected_events = expected_events_uniform(0.75, Some(0.25), 8., 0., 10., f64::INFINITY);

        assert_eq!(expected_events, 60.);
    }

    #[test]
    fn uniform_enrollment_1() {
        // 1:1 randomization, no dropout
        let expected_events_control =
            0.5 * expected_events_uniform(0.036, None, 8., 0., 22., 38.296);
        let expected_events_treat = 0.5 * expected_events_uniform(0.018, None, 8., 0., 22., 38.296);
        assert!((expected_events_control + expected_events_treat - 88.).abs() < 0.001)
    }

    #[test]
    fn uniform_enrollment_2() {
        // 2:1 randomization, no dropout
        let expected_events_control =
            (1. / 3.) * expected_events_uniform(0.036, None, 8., 0., 24.75, 43.155);
        let expected_events_treat =
            (2. / 3.) * expected_events_uniform(0.018, None, 8., 0., 24.75, 43.155);
        assert!((expected_events_control + expected_events_treat - 99.).abs() < 0.001)
    }

    #[test]
    fn uniform_with_dropout_1() {
        // 2:1 randomization, with dropout
        let expected_events_control =
            (1. / 3.) * expected_events_uniform(0.036, Some(0.00878), 8., 0., 28.25, 44.369);
        let expected_events_treat =
            (2. / 3.) * expected_events_uniform(0.018, Some(0.00878), 8., 0., 28.25, 44.369);
        assert!((expected_events_control + expected_events_treat - 99.).abs() < 0.001)
    }

    #[test]
    fn piecewise_no_dropout_1() {
        let enrollment_rate = EnrollmentRate::new(vec![0., 1., 2.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");
        let expected_events_ctrl =
            0.5 * expected_events_piecewise(0.036, None, &enrollment_rate, 10.367, 19.974);
        let expected_events_trt =
            0.5 * expected_events_piecewise(0.018, None, &enrollment_rate, 10.367, 19.974);

        assert!((expected_events_ctrl + expected_events_trt - 88.).abs() < 0.01)
    }

    #[test]
    fn piecewise_no_dropout_2() {
        let enrollment_rate = EnrollmentRate::new(vec![0., 1., 2.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");
        let expected_events_ctrl =
            (1. / 3.) * expected_events_piecewise(0.036, None, &enrollment_rate, 11.533, 22.39);
        let expected_events_trt =
            (2. / 3.) * expected_events_piecewise(0.018, None, &enrollment_rate, 11.533, 22.39);

        assert!((expected_events_ctrl + expected_events_trt - 99.).abs() < 0.01)
    }

    #[test]
    fn piecewise_with_dropout_1() {
        let enrollment_rate = EnrollmentRate::new(vec![0., 1., 2.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");
        let expected_events_ctrl = (1. / 3.)
            * expected_events_piecewise(0.036, Some(0.00878), &enrollment_rate, 12.433, 22.545);
        let expected_events_trt = (2. / 3.)
            * expected_events_piecewise(0.018, Some(0.00878), &enrollment_rate, 12.433, 22.545);

        assert!((expected_events_ctrl + expected_events_trt - 99.).abs() < 0.01)
    }

    #[test]
    fn piecewise_with_dropout_2() {
        let enrollment_rate = EnrollmentRate::new(vec![0., 5., 10.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");
        let expected_events_ctrl =
            0.5 * expected_events_piecewise(0.036, Some(0.00878), &enrollment_rate, 14.9, 23.652);
        let expected_events_trt =
            0.5 * expected_events_piecewise(0.018, Some(0.00878), &enrollment_rate, 14.9, 23.652);

        assert!((expected_events_ctrl + expected_events_trt - 88.).abs() < 0.001);
    }

    #[test]
    fn piecewise_with_dropout_3() {
        let enrollment_rate = EnrollmentRate::new(vec![0., 5., 10.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");
        let expected_events_ctrl =
            0.5 * expected_events_piecewise(0.06, Some(0.00878), &enrollment_rate, 21.267, 30.075);
        let expected_events_trt =
            0.5 * expected_events_piecewise(0.04, Some(0.00878), &enrollment_rate, 21.267, 30.075);

        assert!((expected_events_ctrl + expected_events_trt - 256.).abs() < 0.01);
    }
}
