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
    while f(upper_bound) < target {
        upper_bound *= 2.;
        upper_bound += 1.; // In case lower_bound is zero
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
}
