use crate::error::CtcomputeErr;
use crate::information::exit_probability::exit_probability;
use crate::information::{
    error::TrialBoundsError, std_normal::std_normal_quantile, types::IntegralType,
};
use crate::spending::types::AlphaSpendingValues;

// Information fractions, alpha spend, and quadrature
// size (r), finds bounds that satisfy constraints
// tol is tolerance, i.e. how close to true alpha should
// we shoot for
#[allow(non_snake_case)]
pub fn find_bounds(
    alpha: &AlphaSpendingValues,
    look_fractions: &Vec<f64>,
    r: usize,
    tol: f64,
) -> Result<Vec<(f64, f64)>, CtcomputeErr> {
    let mut bounds: Vec<(f64, f64)> =
        vec![(f64::NEG_INFINITY, f64::INFINITY); look_fractions.len()];

    let (maybe_lower_alpha, maybe_upper_alpha) = match alpha {
        AlphaSpendingValues::OneSidedLower(v) => (Some(v), None),
        AlphaSpendingValues::OneSidedUpper(v) => (None, Some(v)),
        AlphaSpendingValues::TwoSided((v_l, v_u)) => (Some(v_l), Some(v_u)),
    };

    for integral_type in [IntegralType::Lower, IntegralType::Upper] {
        let maybe_alpha = match integral_type {
            IntegralType::Lower => maybe_lower_alpha,
            IntegralType::Upper => maybe_upper_alpha,
        };
        // If no bound was specified, default is +- infinity,
        // so just skip to "other side"
        if maybe_alpha.is_none() {
            continue;
        }
        let alpha = maybe_alpha.unwrap();

        let mut alpha_increments = vec![];
        alpha_increments.reserve(look_fractions.len());

        // First alpha increment is just the first amount of alpha spent
        alpha_increments.push(alpha[0]);
        for i in 1..alpha.len() {
            alpha_increments.push(alpha[i] - alpha[i - 1]);
        }

        // First bound is just \Phi^-1(alpha_increments[0]) or
        // \Phi^-1(1 - alpha_increments[0]), depending on lower vs
        // upper
        if integral_type == IntegralType::Lower {
            let bound_1 = std_normal_quantile(alpha_increments[0])?;
            bounds[0].0 = bound_1;
        } else {
            let bound_1 = std_normal_quantile(1.0 - alpha_increments[0])?;
            bounds[0].1 = bound_1;
        }

        for i in 1..alpha_increments.len() {
            // Assuming we are searching between 0, 7, so start
            // in the middle
            const ABS_MAX_Z: f64 = 7.;
            let (mut lower_bound, mut upper_bound) = match integral_type {
                IntegralType::Lower => (-ABS_MAX_Z, 0.0),
                IntegralType::Upper => (0.0, ABS_MAX_Z),
            };
            let mut mid = (lower_bound + upper_bound) / 2.0;
            match integral_type {
                IntegralType::Lower => bounds[i].0 = mid,
                IntegralType::Upper => bounds[i].1 = mid,
            };

            let mut cur_alpha: f64 = exit_probability(
                &bounds[0..=i],
                &look_fractions[0..=i],
                0.0,
                integral_type,
                r,
            )[i];
            let target_alpha = alpha_increments[i];
            let mut diff: f64 = target_alpha - cur_alpha;

            // Iterate until we get close enought to target alpha
            while diff.abs() > tol && (lower_bound - upper_bound).abs() > tol {
                if integral_type == IntegralType::Upper && cur_alpha <= target_alpha {
                    upper_bound = mid;
                } else if integral_type == IntegralType::Upper && cur_alpha > target_alpha {
                    lower_bound = mid;
                } else if integral_type == IntegralType::Lower && cur_alpha <= target_alpha {
                    lower_bound = mid;
                } else if integral_type == IntegralType::Lower && cur_alpha > target_alpha {
                    upper_bound = mid;
                } else {
                    unreachable!()
                }
                mid = (lower_bound + upper_bound) / 2.0;

                // Update bounds
                match integral_type {
                    IntegralType::Lower => bounds[i].0 = mid,
                    IntegralType::Upper => bounds[i].1 = mid,
                };

                // compute new alpha
                cur_alpha = exit_probability(
                    &bounds[0..=i],
                    &look_fractions[0..=i],
                    0.0,
                    integral_type,
                    r,
                )[i];

                diff = target_alpha - cur_alpha;
            }

            if diff.abs() > tol {
                return Err(
                    TrialBoundsError::FailedToConverge(cur_alpha, target_alpha, tol).into(),
                );
            }
        }
    }

    Ok(bounds)
}

#[cfg(test)]
mod tests {
    use crate::spending::{spending_fcns::compute_spending_vec, types::SpendingFcn};

    use super::*;

    #[test]
    fn pocock_2() {
        let alpha: f64 = exit_probability(
            &vec![(-2.17878788, 2.17878788), (-2.17878788, 2.17878788)],
            &vec![0.5, 1.0],
            0.0,
            IntegralType::Lower,
            32,
        )
        .iter()
        .sum();

        assert!((alpha - 0.025).abs() < 0.0001)
    }

    #[test]
    fn pocock_3() {
        let alpha: f64 = exit_probability(
            &vec![(-2.289, 2.289), (-2.289, 2.289), (-2.289, 2.289)],
            &vec![1.0 / 3.0, 2.0 / 3.0, 1.0],
            0.0,
            IntegralType::Lower,
            32,
        )
        .iter()
        .sum();
        assert!((alpha - 0.025).abs() < 0.0001)
    }

    #[test]
    fn alpha_ld_of_1() {
        let alpha: f64 = exit_probability(
            &vec![(-2.437995, 2.437995), (-1.999873, 1.999873)],
            &vec![0.7, 1.0],
            0.0,
            IntegralType::Lower,
            32,
        )
        .iter()
        .sum();
        assert!((alpha - 0.025).abs() < 0.0001)
    }
    #[test]
    fn alpha_ld_of_2() {
        let alpha: f64 = exit_probability(
            &vec![(-3.356869, 3.356869), (-1.962261, 1.962261)],
            &vec![0.4, 1.0],
            0.0,
            IntegralType::Lower,
            32,
        )
        .iter()
        .sum();
        assert!((alpha - 0.025).abs() < 0.0001)
    }

    #[test]
    fn ldof_bounds_2_looks() {
        // This is the cumulative alpha spent for LDOF per ldBounds(c(0.7, 1.0))
        let alpha_spend = compute_spending_vec(
            &vec![0.7, 1.0],
            0.025,
            Some(&SpendingFcn::LDOF),
            Some(&SpendingFcn::LDOF),
        )
        .unwrap();
        let bounds = find_bounds(&alpha_spend, &vec![0.7, 1.0], 32, 0.0001).unwrap();

        assert!((bounds[0].1 - 2.437995).abs() < 0.001);
        assert!((bounds[1].1 - 1.999930).abs() < 0.001);
    }

    #[test]
    fn ldof_bounds_3_looks() {
        // This is the cumulative alpha spent for LDOF per ldBounds(c(0.3, 0.6, 1.0))
        let alpha_spend = compute_spending_vec(
            &vec![0.3, 0.6, 1.0],
            0.025,
            Some(&SpendingFcn::LDOF),
            Some(&SpendingFcn::LDOF),
        )
        .unwrap();
        let bounds = find_bounds(&alpha_spend, &vec![0.3, 0.6, 1.0], 32, 0.00001).unwrap();

        assert!((bounds[0].1 - 3.928573).abs() < 0.001);
        assert!((bounds[1].1 - 2.669967).abs() < 0.001);
        assert!((bounds[2].1 - 1.981004).abs() < 0.001);
    }

    #[test]
    fn ldof_bounds_3_looks_2() {
        // This is the cumulative alpha spent for LDOF per ldBounds(c(0.3, 0.6, 1.0))
        let alpha_spend = compute_spending_vec(
            &vec![1. / 3., 2. / 3., 1.0],
            0.025,
            Some(&SpendingFcn::LDOF),
            None,
        )
        .unwrap();

        let bounds = find_bounds(&alpha_spend, &vec![1. / 3., 2. / 3., 1.], 32, 0.00001).unwrap();

        // These bounds are from gsDesign
        assert!((bounds[0].0 - -3.710303).abs() < 0.001);
        assert!((bounds[1].0 - -2.511427).abs() < 0.001);
        assert!((bounds[2].0 - -1.993048).abs() < 0.001);
    }

    #[test]
    fn ldof_bounds_3_looks_3() {
        // This is the cumulative alpha spent for LDOF per ldBounds(c(0.3, 0.6, 1.0))
        let alpha_spend =
            compute_spending_vec(&vec![0.4, 0.8, 1.0], 0.025, Some(&SpendingFcn::LDOF), None)
                .unwrap();

        let bounds = find_bounds(&alpha_spend, &vec![0.4, 0.8, 1.], 32, 0.00001).unwrap();

        // These are from EAST; changed info fracs due to weird numerical issues
        // with 1/3, 2/3
        assert!((bounds[0].0 - -3.356869).abs() < 0.001);
        assert!((bounds[1].0 - -2.254642).abs() < 0.001);
        assert!((bounds[2].0 - -2.025815).abs() < 0.001);
    }
}
