use crate::duration::types::EnrollmentRate;
use crate::error::CtcomputeErr;

/// Given a desired sample size and EnrollmentRate object, computes expected
/// enrollment duration
pub fn expected_enrollment_dur(
    n: usize,
    enrollment_rate: &EnrollmentRate,
) -> Result<f64, CtcomputeErr> {
    let n_float = n as f64;
    // For now, this is the only type of enrollment
    let EnrollmentRate::PiecewiseConstant {
        enrollment_durations,
        enrollment_rates,
    } = enrollment_rate;
    let mut cur_time = 0.;
    let mut cur_n = 0.;

    // This assumes final duration is infinite
    for (cur_duration, cur_rate) in enrollment_durations.iter().zip(enrollment_rates.iter()) {
        if cur_n + cur_duration * cur_rate < n_float {
            cur_n += cur_duration * cur_rate;
            cur_time += cur_duration;
            continue;
        }

        let n_still_needed = n_float - cur_n;

        let time_still_needed = n_still_needed / cur_rate;
        cur_time += time_still_needed;

        break;
    }
    return Ok(cur_time);
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::duration::types::EnrollmentRate;

    #[test]
    fn expected_enrollment_dur_1() {
        let enrollment_rate =
            EnrollmentRate::new(vec![0.], vec![4.]).expect("failed to construct EnrollmentRate");
        let expected_enrollment_dur = expected_enrollment_dur(4, &enrollment_rate)
            .expect("failed to compute expected enrollment duration");
        assert_eq!(expected_enrollment_dur, 1.0);
    }

    #[test]
    fn expected_enrollment_dur_2() {
        let enrollment_rate =
            EnrollmentRate::new(vec![5.], vec![4.]).expect("failed to construct EnrollmentRate");
        let expected_enrollment_dur = expected_enrollment_dur(4, &enrollment_rate)
            .expect("failed to compute expected enrollment duration");
        assert_eq!(expected_enrollment_dur, 6.);
    }

    #[test]
    fn excessive_enrollment_dur() {
        let enrollment_rate = EnrollmentRate::new(
            vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        )
        .expect("failed to construct EnrollmentRate");
        let expected_enrollment_dur = expected_enrollment_dur(10, &enrollment_rate)
            .expect("failed to compute expected enrollment duration");
        assert_eq!(expected_enrollment_dur, 4.);
    }
}
