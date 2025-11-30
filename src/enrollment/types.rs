use std::f64;

use crate::enrollment::error::EnrollmentComputationError;
use crate::error::CtcomputeErr;

#[derive(Debug)]
pub enum EnrollmentRate {
    PiecewiseConstant {
        enrollment_durations: Vec<f64>,
        enrollment_rates: Vec<f64>,
    },
}

impl EnrollmentRate {
    pub fn new(times: Vec<f64>, rates: Vec<f64>) -> Result<Self, CtcomputeErr> {
        // Rates should have equal length
        if rates.len() != times.len() {
            return Err(EnrollmentComputationError::TimeRateLengths {
                time_length: times.len(),
                rate_length: rates.len(),
            }
            .into());
        }

        // TODO: Also check:
        // - total specified enrollment > 0
        // - Final enrollment rate > 0
        // - All rates/times should be >= 0
        // - Various things related to infinite values

        //----------------------------------------
        // If time doesn't start at 0, starting rate is implicitly zero
        let enrollment_rates = match times[0] {
            0.0 => rates,
            _ => [vec![0.0], rates.clone()].concat(),
        };

        // Set up enrollment rates/times
        // If first time isn't 0, add 0
        let enrollment_times = match times[0] {
            0.0 => times,
            _ => [vec![0.0], times.clone()].concat(),
        };

        let mut durations: Vec<f64> = enrollment_times.windows(2).map(|w| w[1] - w[0]).collect();
        // Last enrollment rate gets extended forever
        durations.push(f64::INFINITY);

        Ok(EnrollmentRate::PiecewiseConstant {
            enrollment_durations: durations,
            enrollment_rates: enrollment_rates,
        })
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn basic_piecewise_enrollment() {
        let piecewise = EnrollmentRate::new(vec![0., 1., 2.], vec![4., 5., 6.]);
        if let Ok(EnrollmentRate::PiecewiseConstant {
            enrollment_durations: durations,
            enrollment_rates: rates,
        }) = piecewise
        {
            assert_eq!(durations, vec![1., 1., f64::INFINITY]);
            assert_eq!(rates, vec![4., 5., 6.]);
        } else {
            panic!()
        }
    }

    #[test]
    fn uniform_enrollment() {
        let piecewise = EnrollmentRate::new(vec![0.], vec![4.]);
        if let Ok(EnrollmentRate::PiecewiseConstant {
            enrollment_durations: durations,
            enrollment_rates: rates,
        }) = piecewise
        {
            assert_eq!(durations, vec![f64::INFINITY]);
            assert_eq!(rates, vec![4.]);
        } else {
            panic!()
        }
    }

    #[test]
    fn mismatched_lengths_error() {
        if let Err(e) = EnrollmentRate::new(vec![0.0, 1.0], vec![1.0]) {
            assert_eq!(
                String::from(
                    "while computing enrollment: \
                     lengths of enrollment times and rates don't match (times length 2, \
                     rates length 1"
                ),
                format!("{}", e)
            );
        } else {
            panic!()
        }
    }
}
