use crate::{ctsim_err::CtsimErr, quadrature::Quadrature, spending_fcns::AlphaSpendingValues};
use probability_rs::{Continuous, Distribution, dist::normal::Normal};
use thiserror::Error;

// TODO: incorporate hypothesis sidedness
#[derive(Error, Debug)]
pub enum NormalDistErr {
    #[error("arguments to quantile function should be in [0, 1]; got {0}")]
    QuantileOutOfBounds(f64),
}

#[derive(PartialEq, Clone, Copy)]
pub enum IntegralType {
    Upper,
    Lower,
}

impl Into<CtsimErr> for NormalDistErr {
    fn into(self) -> CtsimErr {
        CtsimErr::NormalDist(self)
    }
}

pub fn std_normal_pdf(z: f64) -> f64 {
    (-z * z / 2.0).exp() / (2.0 * std::f64::consts::PI).sqrt()
}

pub fn std_normal_cdf(z: f64) -> f64 {
    let std_normal = Normal::new(0.0, 1.0).unwrap();
    std_normal.cdf(z)
}

#[allow(clippy::excessive_precision)]
pub fn std_normal_quantile_helper(p: f64) -> f64 {
    assert!(p > 0.0 && p < 1.0, "p must be in (0,1)");

    // Coefficients (Acklam 2003). See public documentation.
    const A: [f64; 6] = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.383577518672690e+02,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    ];
    const B: [f64; 5] = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ];
    const C: [f64; 6] = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    ];
    const D: [f64; 4] = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    ];
    const P_LOW: f64 = 0.02425;
    const P_HIGH: f64 = 1.0 - P_LOW;
    if p < P_LOW {
        // Lower tail region
        let q = (-2.0 * p.ln()).sqrt();
        let x = (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0);
        return x;
    }
    if p > P_HIGH {
        // Upper tail region
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        let x = (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0);
        return -x;
    }
    // Central region
    let q = p - 0.5;
    let r = q * q;
    (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
        / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
}
pub fn std_normal_quantile(p: f64) -> Result<f64, CtsimErr> {
    if p < 0.0 || p > 1.0 {
        return Err(NormalDistErr::QuantileOutOfBounds(p).into());
    }
    // let std_normal = Normal::new(0.0, 1.0).unwrap();
    // For some reason their quantile function is inverted, i.e.
    // 1st percentile is actually 99th, etc.
    let z_q = std_normal_quantile_helper(p);
    Ok(z_q)
}

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

#[derive(Error, Debug)]
pub enum TrialBoundsError {
    #[error("failed to converge (computed alpha: {0}, target alpha: {1}, tolerance: {2}")]
    FailedToConverge(f64, f64, f64),
}

impl Into<CtsimErr> for TrialBoundsError {
    fn into(self) -> CtsimErr {
        CtsimErr::TrialBounds(self)
    }
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
        // println!("alpha: {alpha:?}");
        // println!("alpha_increments: {alpha_increments:?}");

        // First bound is just \Phi^-1(alpha_increments[0]) or
        // \Phi^-1(1 - alpha_increments[0]), depending on lower vs
        // upper
        // println!("alpha_increments: {:?}", alpha_increments);
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
            // println!("---------------------------------------");
            // println!("mid: {mid}");
            // println!("bounds {:?}", bounds);
            // println!("&I[0..=i] {:?}", &look_fractions[0..=i]);
            // println!("cur_alpha: {cur_alpha}");
            // println!("target_alpha: {target_alpha}");

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
                // println!("---------------------------------------");
                // println!("mid: {mid}");
                // println!("bounds {:?}", &bounds[0..=i]);
                // println!("&I[0..=i] {:?}", &look_fractions[0..=i]);
                // println!("cur_alpha: {cur_alpha}");
                // println!("target_alpha: {target_alpha}");

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
    use crate::spending_fcns::{SpendingFcn, compute_spending_vec};

    use super::*;

    #[test]
    fn standard_normal_pdf_1() {
        assert!((0.3989423 - std_normal_pdf(0.0)).abs() < 0.0000001)
    }

    #[test]
    fn standard_normal_pdf_2() {
        assert!((0.004431848 - std_normal_pdf(3.0)).abs() < 0.0000001)
    }

    #[test]
    fn standard_normal_pdf_3() {
        assert_eq!(std_normal_pdf(-1.0), std_normal_pdf(1.0));
    }

    #[test]
    fn std_normal_quantile_err() {
        if let Err(e) = std_normal_quantile(1.1) {
            assert_eq!(
                String::from(
                    "while evaluating normal distribution: arguments to \
                    quantile function should be in [0, 1]; got 1.1"
                ),
                format!("{}", e)
            );
        } else {
            panic!()
        }
    }

    #[test]
    fn std_normal_quantile_value() {
        assert!((std_normal_quantile(0.975).unwrap() - 1.96).abs() < 0.0001)
    }

    #[test]
    fn std_normal_quantile_value_2() {
        assert!((std_normal_quantile(0.007384489).unwrap() - -2.437995).abs() < 0.0001)
    }

    #[test]
    fn std_normal_quantile_symmetric() {
        assert_eq!(
            std_normal_quantile(0.975).unwrap(),
            -std_normal_quantile(0.025).unwrap()
        )
    }

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
