use crate::error::CtsimErr;
use crate::integration::{
    error::TrialBoundsError,
    quadrature::Quadrature,
    std_normal::{std_normal_cdf, std_normal_pdf, std_normal_quantile},
    types::IntegralType,
};
use crate::spending::types::AlphaSpendingValues;

// J&T p. 354
#[allow(non_snake_case)]
pub fn f_1(z_1: f64, I_1: f64, theta: f64) -> f64 {
    std_normal_pdf(z_1 - theta * I_1.sqrt())
}

#[allow(non_snake_case)]
pub fn f_k(z_k_1: f64, I_k_1: f64, z_k: f64, theta: f64, I_k: f64) -> f64 {
    let delta_k = I_k - I_k_1;
    let sqrt_delta_k = delta_k.sqrt();
    (I_k.sqrt() / sqrt_delta_k)
        * std_normal_pdf((z_k * I_k.sqrt() - z_k_1 * I_k_1.sqrt() - theta * delta_k) / sqrt_delta_k)
}

// e_k in J&T's notation
#[allow(non_snake_case)]
pub fn upper_exit_prob(z_k_1: f64, I_k_1: f64, b_k: f64, I_k: f64, theta: f64) -> f64 {
    let delta_k = I_k - I_k_1;
    std_normal_cdf((z_k_1 * I_k_1.sqrt() + theta * delta_k - b_k * I_k.sqrt()) / delta_k.sqrt())
}

// Defined to capture probability of being less than a; no
// corresponding definition in J&T
#[allow(non_snake_case)]
pub fn lower_exit_prob(z_k_1: f64, I_k_1: f64, a_k: f64, I_k: f64, theta: f64) -> f64 {
    let delta_k = I_k - I_k_1;
    std_normal_cdf((a_k * I_k.sqrt() - z_k_1 * I_k_1.sqrt() - theta * delta_k) / delta_k.sqrt())
}

// J&T, p. 355
// TODO: rename, since this is technically psi_k and xi_k
// Computes exit probabilities given bounds, information fractions,
// theta, and integral type (upper or lower). r controls quadrature size.
#[allow(non_snake_case)]
pub fn psi_k(
    bounds: &[(f64, f64)],
    look_fractions: &[f64],
    theta: f64,
    integral_type: IntegralType,
    r: usize,
) -> Vec<f64> {
    //----------------------------------------
    // Quadratures
    let quadratures: Vec<Quadrature> = bounds
        .iter()
        .zip(look_fractions.iter())
        .map(|(&(a_k, b_k), &I_k)| Quadrature::new(theta, I_k, r, a_k, b_k))
        .collect();

    //----------------------------------------
    // h_k values
    let mut h_k: Vec<Vec<f64>> = vec![];
    h_k.reserve(quadratures.len());
    for i in 0..quadratures.len() {
        let v = if i == 0 {
            quadratures[i]
                .w
                .iter()
                .zip(quadratures[i].z.iter())
                .map(|(&w_k, &z_k)| w_k * f_1(z_k, look_fractions[i], theta))
                .collect()
        } else {
            quadratures[i]
                .w
                .iter()
                .zip(quadratures[i].z.iter())
                .map(|(&w_k, &z_k)| {
                    // Vector of h_{k-1} values
                    let h_k_1_vec = &h_k[i - 1];
                    // Quadrature of look k-1
                    let q_k_1 = &quadratures[i - 1];
                    let I_k_1 = look_fractions[i - 1];
                    let I_k = look_fractions[i];
                    h_k_1_vec
                        .iter()
                        .zip(q_k_1.z.iter())
                        .map(|(h_k_1, &z_k_1)| h_k_1 * w_k * f_k(z_k_1, I_k_1, z_k, theta, I_k))
                        .sum()
                })
                .collect()
        };
        h_k.push(v);
    }

    let mut res = vec![];
    res.reserve(look_fractions.len());

    //----------------------------------------
    // Compute integral

    // First integral just uses cdf of N(theta(sqrt(I_1)), 1)
    res.push(match integral_type {
        IntegralType::Lower => std_normal_cdf(bounds[0].0 - theta * look_fractions[0].sqrt()),
        IntegralType::Upper => 1.0 - std_normal_cdf(bounds[0].1 - theta * look_fractions[0].sqrt()),
    });

    // whether or not the final integral should be probability that
    // statistic is greater than upper bound or lower than lower bound
    let final_integral_type = match integral_type {
        IntegralType::Upper => upper_exit_prob,
        IntegralType::Lower => lower_exit_prob,
    };

    for i in 1..look_fractions.len() {
        let final_bound = match integral_type {
            IntegralType::Upper => bounds[i].1,
            IntegralType::Lower => bounds[i].0,
        };
        let I_k_1 = look_fractions[i - 1];
        let I_k = look_fractions[i];
        let prob: f64 = h_k[i - 1]
            .iter()
            .zip(quadratures[i - 1].z.iter())
            .map(|(&h_k_1, &z_k_1)| {
                h_k_1 * final_integral_type(z_k_1, I_k_1, final_bound, I_k, theta)
            })
            .sum();

        res.push(prob);
    }

    res
}

// Information fractions, alpha spend, and quadrature
// size (r), finds bounds that satisfy constraints
// tol is tolerance, i.e. how close to true alpha should
// we shoot for
// TODO: update alpha to be split into sides, i.e. allow
// for unequal allocation
#[allow(non_snake_case)]
pub fn find_bounds(
    alpha: &AlphaSpendingValues,
    look_fractions: &Vec<f64>,
    r: usize,
    tol: f64,
) -> Result<Vec<(f64, f64)>, CtsimErr> {
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
            let (mut lower_bound, mut upper_bound) = match integral_type {
                IntegralType::Lower => (-7.0, 0.0),
                IntegralType::Upper => (0.0, 7.0),
            };
            let mut mid = (lower_bound + upper_bound) / 2.0;
            match integral_type {
                IntegralType::Lower => bounds[i].0 = mid,
                IntegralType::Upper => bounds[i].1 = mid,
            };

            let mut cur_alpha: f64 = psi_k(
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
                cur_alpha = psi_k(
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
        let alpha: f64 = psi_k(
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
        let alpha: f64 = psi_k(
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
        let alpha: f64 = psi_k(
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
        let alpha: f64 = psi_k(
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
            Some(SpendingFcn::LDOF),
            Some(SpendingFcn::LDOF),
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
            Some(SpendingFcn::LDOF),
            Some(SpendingFcn::LDOF),
        )
        .unwrap();
        let bounds = find_bounds(&alpha_spend, &vec![0.3, 0.6, 1.0], 32, 0.00001).unwrap();

        assert!((bounds[0].1 - 3.928573).abs() < 0.001);
        assert!((bounds[1].1 - 2.669967).abs() < 0.001);
        assert!((bounds[2].1 - 1.981004).abs() < 0.001);
    }
}
