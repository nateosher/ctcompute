use crate::error::{CtcomputeErr, RootFindErr};

// TODO: add failsafes to while loops

/// Given a monotonically increasing function f(x) and lower bound, finds
/// value x' to the right of the lower bound such that f(x') = target
pub fn root_find_monotonic<F>(
    f: F,
    lower_bound: f64,
    target: f64,
    tol: f64,
) -> Result<f64, CtcomputeErr>
where
    F: Fn(f64) -> f64,
{
    if f(lower_bound) >= target {
        return Err(RootFindErr::BadLowerBound.into());
    }
    // Set window for search
    let mut lower_bound = lower_bound;
    let mut upper_bound = lower_bound;
    let mut safety = 0;
    let mut f_upper_bound = f(upper_bound);
    while f_upper_bound < target && safety < 10 {
        upper_bound *= 2.;
        upper_bound += 1.; // In case lower_bound is zero
        f_upper_bound = f(upper_bound);
        safety += 1;
    }

    // Perform search
    let mut x = (lower_bound + upper_bound) / 2.;
    let mut y = f(x);
    while (lower_bound - upper_bound).abs() > tol / 2. && (y - target).abs() > tol {
        if y <= target {
            lower_bound = x;
        } else {
            upper_bound = x;
        }
        x = (lower_bound + upper_bound) / 2.;
        y = f(x);
    }
    Ok(x)
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::duration::types::EnrollmentRate;
    use crate::events::expected_events::expected_events_piecewise_arms;

    #[test]
    fn basic_linear_root_find() {
        let f = |x| x;
        let res =
            root_find_monotonic(f, 0.0, 3., 0.001).expect("failed to perform linear root find");
        assert!((res - 3.0).abs() < 0.001);
    }

    #[test]
    fn basic_quadratic_root_find() {
        let f = |x| x * x;
        let res =
            root_find_monotonic(f, 0.0, 9., 0.001).expect("failed to perform quadratic root find");
        assert!((res - 3.0).abs() < 0.001);
    }

    #[test]
    fn problem_case() {
        let er = EnrollmentRate::new(vec![0., 5., 10.], vec![10., 20., 30.])
            .expect("failed to construct enrollment rate object");

        let f_expected_by_followup_and_accrual = |followup_dur| {
            expected_events_piecewise_arms(
                2. / 3.,
                Some(0.00878),
                0.04,
                0.06,
                &er,
                followup_dur, // accrual time is set to followup time as well
                followup_dur,
            )
        };
        let res = f_expected_by_followup_and_accrual(17.072359168529523);
        println!("res: {res:?}");
        let res_2 = f_expected_by_followup_and_accrual(0.);
        println!("res_2: {res_2:?}");

        let rf = root_find_monotonic(f_expected_by_followup_and_accrual, 0., 292., 0.001);
        println!("rf: {rf:?}");
    }
}
