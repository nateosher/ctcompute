use crate::duration::expected_enrollment::expected_enrollment;
use crate::duration::{
    expected_enrollment, expected_enrollment_dur::expected_enrollment_dur,
    min_followup::min_followup, types::EnrollmentRate,
};
use crate::error::{CtcomputeErr, TrialComputeError};
use crate::events::events_needed::events_needed;
use crate::information::{
    compute_information::compute_information, exit_probability::exit_probability,
    trial_bounds::find_bounds, types::IntegralType,
};
use crate::spending::{spending_fcns::compute_spending_vec, types::SpendingFcn};
use crate::trial::types::Trial;

pub fn compute_trial(
    n_patients: usize,
    alpha: f64,
    power: f64,
    maybe_lower_spending_fcn_type: Option<&SpendingFcn>,
    maybe_upper_spending_fcn_type: Option<&SpendingFcn>,
    maybe_look_fractions: Option<&Vec<f64>>,
    prop_treated: f64,
    lambda_event_trt: f64,
    lambda_event_ctrl: f64,
    maybe_lambda_dropout: Option<f64>,
    enrollment_rate: &EnrollmentRate,
    r: usize,
    tol: f64,
) -> Result<Trial, CtcomputeErr> {
    //----------------------------------------
    // Check information + sample size + event requirements
    //----------------------------------------
    let hr = lambda_event_trt / lambda_event_ctrl;
    let log_hr = hr.ln();

    #[allow(non_snake_case)]
    let I = compute_information(
        alpha,
        power,
        log_hr,
        maybe_lower_spending_fcn_type,
        maybe_upper_spending_fcn_type,
        maybe_look_fractions,
        tol,
    )?;

    let n_events = events_needed(
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

    if n_patients < n_events as usize {
        return Err(TrialComputeError::InsufficientSampleSize {
            n_patients: n_patients,
            n_events_needed: n_events as usize,
        }
        .into());
    }

    //----------------------------------------
    // Accrual duration
    //----------------------------------------
    let accrual_duration = expected_enrollment_dur(n_patients, enrollment_rate)?;

    //----------------------------------------
    // Followup duration
    //----------------------------------------

    let trial_duration = min_followup(
        n_events,
        prop_treated,
        lambda_event_trt,
        lambda_event_ctrl,
        maybe_lambda_dropout,
        enrollment_rate,
        accrual_duration,
        tol,
    )?;

    //----------------------------------------
    // Exit probabilities by look
    //----------------------------------------

    let maybe_events_by_look: Option<Vec<f64>> = match maybe_look_fractions {
        None => None,
        Some(look_fractions) => Some(
            look_fractions
                .iter()
                .map(|f| (f * n_events).ceil())
                .collect(),
        ),
    };
    dbg!(&maybe_events_by_look);

    let maybe_h1_times_to_looks = if let Some(ref events_by_look) = maybe_events_by_look {
        let times_to_looks = events_by_look
            .iter()
            .map(|n_events| {
                min_followup(
                    *n_events,
                    prop_treated,
                    lambda_event_trt,
                    lambda_event_ctrl,
                    maybe_lambda_dropout,
                    enrollment_rate,
                    accrual_duration,
                    tol,
                )
            })
            .collect::<Result<Vec<f64>, CtcomputeErr>>()?;
        Some(times_to_looks)
    } else {
        None
    };
    dbg!(&maybe_h1_times_to_looks);

    let maybe_h0_times_to_looks = if let Some(ref events_by_look) = maybe_events_by_look {
        let times_to_looks = events_by_look
            .iter()
            .map(|n_events| {
                min_followup(
                    *n_events,
                    prop_treated,
                    lambda_event_ctrl, // assume both have event rate of control
                    lambda_event_ctrl,
                    maybe_lambda_dropout,
                    enrollment_rate,
                    accrual_duration,
                    tol,
                )
            })
            .collect::<Result<Vec<f64>, CtcomputeErr>>()?;
        Some(times_to_looks)
    } else {
        None
    };

    let (maybe_h1_exit_probs, maybe_h0_exit_probs) =
        if let Some(look_fractions) = maybe_look_fractions {
            let alpha_spend = compute_spending_vec(
                look_fractions,
                alpha,
                maybe_lower_spending_fcn_type,
                maybe_upper_spending_fcn_type,
            )?;
            let bounds = find_bounds(&alpha_spend, look_fractions, r, tol)?;

            let hr = lambda_event_trt / lambda_event_ctrl;
            let log_hr = hr.ln();
            let theta = log_hr * I.sqrt();
            let h1_lower_exit_probs =
                exit_probability(&bounds, look_fractions, theta, IntegralType::Lower, r);
            let h1_upper_exit_probs =
                exit_probability(&bounds, look_fractions, theta, IntegralType::Upper, r);
            let mut h1_exit_probs = h1_lower_exit_probs
                .iter()
                .zip(h1_upper_exit_probs.iter())
                .map(|(p_l, p_u)| p_l + p_u)
                .collect::<Vec<f64>>();

            let h0_lower_exit_probs =
                exit_probability(&bounds, look_fractions, 0., IntegralType::Lower, r);
            let h0_upper_exit_probs =
                exit_probability(&bounds, look_fractions, 0., IntegralType::Upper, r);
            let mut h0_exit_probs = h0_lower_exit_probs
                .iter()
                .zip(h0_upper_exit_probs.iter())
                .map(|(p_l, p_u)| p_l + p_u)
                .collect::<Vec<f64>>();

            // Final "exit" probabilities are complement of exiting
            // at all other times, since trial must end at last look
            let n_looks = h1_exit_probs.len();
            h1_exit_probs[n_looks - 1] = 1.0
                - h1_exit_probs[0..h1_exit_probs.len() - 1]
                    .iter()
                    .sum::<f64>();

            h0_exit_probs[n_looks - 1] = 1.0
                - h0_exit_probs[0..h0_exit_probs.len() - 1]
                    .iter()
                    .sum::<f64>();

            (Some(h1_exit_probs), Some(h0_exit_probs))
        } else {
            (None, None)
        };

    dbg!(&maybe_h0_exit_probs);
    dbg!(&maybe_h1_exit_probs);

    //----------------------------------------
    // Accrual duration under H0/H1
    //----------------------------------------

    let maybe_h0_expected_accrual_duration: Option<f64> = if let Some(ref h0_exit_probs) =
        maybe_h0_exit_probs
        && let Some(ref times_to_looks) = maybe_h0_times_to_looks
    {
        Some(
            h0_exit_probs
                .iter()
                .zip(times_to_looks.iter())
                .map(|(&p_exit, t)| p_exit * (t.min(accrual_duration)))
                .sum(),
        )
    } else {
        None
    };

    let maybe_h1_expected_accrual_duration: Option<f64> = if let Some(ref h1_exit_probs) =
        maybe_h1_exit_probs
        && let Some(ref times_to_looks) = maybe_h1_times_to_looks
    {
        Some(
            h1_exit_probs
                .iter()
                .zip(times_to_looks.iter())
                .map(|(&p_exit, t)| p_exit * (t.min(accrual_duration)))
                .sum(),
        )
    } else {
        None
    };

    //----------------------------------------
    // Trial duration under H0/H1
    //----------------------------------------

    let maybe_h0_expected_trial_duration: Option<f64> = if let Some(ref h0_exit_probs) =
        maybe_h0_exit_probs
        && let Some(ref times_to_looks) = maybe_h0_times_to_looks
    {
        Some(
            h0_exit_probs
                .iter()
                .zip(times_to_looks.iter())
                .map(|(&p_exit, t)| p_exit * t)
                .sum(),
        )
    } else {
        None
    };

    let maybe_h1_expected_trial_duration: Option<f64> = if let Some(ref h1_exit_probs) =
        maybe_h1_exit_probs
        && let Some(ref times_to_looks) = maybe_h1_times_to_looks
    {
        Some(
            h1_exit_probs
                .iter()
                .zip(times_to_looks.iter())
                .map(|(&p_exit, t)| p_exit * t)
                .sum(),
        )
    } else {
        None
    };

    //----------------------------------------
    // Sample size under H0/H1
    //----------------------------------------

    let maybe_h0_expected_sample_size: Option<f64> = if let Some(ref h0_exit_probs) =
        maybe_h0_exit_probs
        && let Some(ref times_to_looks) = maybe_h0_times_to_looks
    {
        let h0_expected_ss_vec_res = h0_exit_probs
            .iter()
            .zip(times_to_looks.iter())
            .map(|(&p_exit, t)| {
                match expected_enrollment(t.min(accrual_duration), &enrollment_rate) {
                    Ok(ee) => Ok(p_exit * ee),
                    Err(e) => Err(e),
                }
            })
            .collect::<Result<Vec<f64>, CtcomputeErr>>()?;
        let h0_expected_ss = h0_expected_ss_vec_res.iter().sum::<f64>();
        Ok(Some(h0_expected_ss))
    } else {
        Ok(None)
    }?;

    let maybe_h1_expected_sample_size: Option<f64> = if let Some(ref h1_exit_probs) =
        maybe_h1_exit_probs
        && let Some(ref times_to_looks) = maybe_h1_times_to_looks
    {
        let h0_expected_ss_vec_res = h1_exit_probs
            .iter()
            .zip(times_to_looks.iter())
            .map(|(&p_exit, t)| {
                match expected_enrollment(t.min(accrual_duration), &enrollment_rate) {
                    Ok(ee) => Ok(p_exit * ee),
                    Err(e) => Err(e),
                }
            })
            .collect::<Result<Vec<f64>, CtcomputeErr>>()?;
        let h0_expected_ss = h0_expected_ss_vec_res.iter().sum::<f64>();
        Ok(Some(h0_expected_ss))
    } else {
        Ok(None)
    }?;

    //----------------------------------------
    // Expected events
    //----------------------------------------
    // TODO

    //----------------------------------------
    // Return result
    //----------------------------------------
    Ok(Trial {
        accrual_duration,
        trial_duration,
        n_events: n_events as usize,
        n_patients,
        maybe_h0_expected_accrual_duration,
        maybe_h1_expected_accrual_duration,
        maybe_h0_expected_trial_duration,
        maybe_h1_expected_trial_duration,
        maybe_h0_expected_sample_size,
        maybe_h1_expected_sample_size,
    })
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn basic_one_look() {
        let er = EnrollmentRate::new(vec![0., 1., 2.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");

        let trial = compute_trial(
            100,   // # patients
            0.025, // alpha
            0.9,   // power
            None,  // lower spending fcn
            None,  // upper spending fcn
            None,  // info fractions
            0.5,   // proportion treated
            0.018, // hazard rate treatment
            0.036, // hazard rate control
            None,  // dropout
            &er,   // enrollment rate
            64,    // r (sets integral grid size)
            0.001, // tol
        );

        println!("trial: {trial:#?}");
    }

    #[test]
    fn excessive_enrollment_three_look() {
        let er = EnrollmentRate::new(vec![0., 10., 20.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");

        let trial = compute_trial(
            400,                        // # patients
            0.025,                      // alpha
            0.9,                        // power
            Some(&SpendingFcn::LDOF),   // lower spending fcn
            None,                       // upper spending fcn
            Some(&vec![0.7, 0.9, 1.0]), // info fractions
            0.5,                        // proportion treated
            0.018,                      // hazard rate treatment
            0.036,                      // hazard rate control
            None,                       // dropout
            &er,                        // enrollment rate
            64,                         // r (sets integral grid size)
            0.00001,                    // tol
        );

        println!("trial: {:#?}", trial);
    }

    #[test]
    fn simple_three_look() {
        let er = EnrollmentRate::new(vec![0.], vec![3.])
            .expect("failed to construct enrollment rate object");

        let trial = compute_trial(
            80,                         // # patients
            0.025,                      // alpha
            0.9,                        // power
            Some(&SpendingFcn::LDOF),   // lower spending fcn
            None,                       // upper spending fcn
            Some(&vec![0.4, 0.8, 1.0]), // info fractions
            0.5,                        // proportion treated
            0.009,                      // hazard rate treatment
            0.036,                      // hazard rate control
            Some(0.00878),              // dropout
            &er,                        // enrollment rate
            32,                         // r (sets integral grid size)
            0.0001,                     // tol
        );

        println!("trial: {:#?}", trial);
    }

    #[test]
    fn simple_three_look_2() {
        let er = EnrollmentRate::new(vec![0., 5.], vec![3., 6.])
            .expect("failed to construct enrollment rate object");

        let trial = compute_trial(
            80,                          // # patients
            0.025,                       // alpha
            0.9,                         // power
            Some(&SpendingFcn::LDOF),    // lower spending fcn
            None,                        // upper spending fcn
            Some(&vec![0.2, 0.75, 1.0]), // info fractions
            2. / 3.,                     // proportion treated
            0.018,                       // hazard rate treatment
            0.036,                       // hazard rate control
            Some(0.00878),               // dropout
            &er,                         // enrollment rate
            32,                          // r (sets integral grid size)
            0.0001,                      // tol
        );

        println!("trial: {:#?}", trial);
    }
}
