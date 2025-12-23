use crate::duration::types::EnrollmentRate;
use crate::error::CtcomputeErr;

/// Given enrollment duration and EnrollmentRate object, computes expected
/// number of patients enrolled
pub fn expected_enrollment(
    dur: f64,
    enrollment_rate: &EnrollmentRate,
) -> Result<f64, CtcomputeErr> {
    let EnrollmentRate::PiecewiseConstant {
        enrollment_durations,
        enrollment_rates,
    } = enrollment_rate;
    let mut cur_time = 0.;
    let mut cur_n = 0.;

    // This assumes final duration is infinite
    for (cur_duration, cur_rate) in enrollment_durations.iter().zip(enrollment_rates.iter()) {
        if cur_time + cur_duration < dur {
            cur_n += cur_duration * cur_rate;
            cur_time += cur_duration;
            continue;
        }
        let remaining_time = dur - cur_time;
        cur_n += remaining_time * cur_rate;
        break;
    }
    return Ok(cur_n);
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn expected_enrollment_1() {
        let enrollment_rate = EnrollmentRate::new(
            vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        )
        .expect("failed to construct EnrollmentRate");
        let expected_n = expected_enrollment(2.0, &enrollment_rate)
            .expect("failed to compute expected enrollment");

        assert_eq!(expected_n, 3.0);
    }

    #[test]
    fn expected_enrollment_2() {
        let enrollment_rate = EnrollmentRate::new(
            vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        )
        .expect("failed to construct EnrollmentRate");
        let expected_n = expected_enrollment(4.0, &enrollment_rate)
            .expect("failed to compute expected enrollment");

        assert_eq!(expected_n, 10.);
    }

    #[test]
    fn expected_enrollment_3() {
        let enrollment_rate = EnrollmentRate::new(vec![0., 5., 10.], vec![10., 20., 30.])
            .expect("failed to construct EnrollmentRate");
        let expected_n = expected_enrollment(14.9, &enrollment_rate)
            .expect("failed to compute expected enrollment");

        assert_eq!(expected_n, 297.);
    }
}
