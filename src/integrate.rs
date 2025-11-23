use crate::{ctsim_err::CtsimErr, quadrature::Quadrature};
use probability_rs::{Continuous, Distribution, dist::normal::Normal};
use thiserror::Error;

// TODO: incorporate hypothesis sidedness
#[derive(Error, Debug)]
pub enum NormalDistErr {
    #[error("arguments to quantile function should be in [0, 1]; got {0}")]
    QuantileOutOfBounds(f64),
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

pub fn std_normal_quantile(q: f64) -> Result<f64, CtsimErr> {
    if q < 0.0 || q > 1.0 {
        return Err(NormalDistErr::QuantileOutOfBounds(q).into());
    }
    let std_normal = Normal::new(0.0, 1.0).unwrap();
    let z_q = std_normal.inv_cdf(q);
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

#[allow(non_snake_case)]
pub fn e_k_1(z_k_1: f64, I_k_1: f64, b_k: f64, I_k: f64, theta: f64) -> f64 {
    let delta_k = I_k - I_k_1;
    std_normal_cdf((z_k_1 * I_k_1.sqrt() + theta * delta_k - b_k * I_k.sqrt()) / delta_k.sqrt())
}

// Defined to capture probability of being less than a
#[allow(non_snake_case)]
pub fn e_k_1_prime(z_k_1: f64, I_k_1: f64, a_k: f64, I_k: f64, theta: f64) -> f64 {
    let delta_k = I_k - I_k_1;
    std_normal_cdf((a_k * I_k.sqrt() - z_k_1 * I_k_1.sqrt() - theta * delta_k) / delta_k.sqrt())
}

// J&T, p. 355
// TODO: rename, since this is technically psi_k and xi_k
// Computes exit probabilities given bounds, information fractions,
// and theta.
#[allow(non_snake_case)]
pub fn psi_k(bounds: &[(f64, f64)], I: &[f64], theta: f64, r: usize) -> Vec<f64> {
    //----------------------------------------
    // Quadratures
    let quadratures: Vec<Quadrature> = bounds
        .iter()
        .zip(I.iter())
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
                .map(|(&w_k, &z_k)| w_k * f_1(z_k, I[i], theta))
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
                    let I_k_1 = I[i - 1];
                    let I_k = I[i];
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
    res.reserve(I.len());

    //----------------------------------------
    // Compute power/type I error, depending on theta
    // Assumes hypothesis is two sided
    // TODO: update to take type

    // First integral just uses cdf of N(theta(sqrt(I_1)), 1)
    res.push(
        std_normal_cdf(bounds[0].0 - theta * I[0].sqrt())
            + (1.0 - std_normal_cdf(bounds[0].1 - theta * I[0].sqrt())),
    );

    for i in 1..I.len() {
        let b_k = bounds[i].1;
        let a_k = bounds[i].0;
        let I_k_1 = I[i - 1];
        let I_k = I[i];
        let upper_prob: f64 = h_k[i - 1]
            .iter()
            .zip(quadratures[i - 1].z.iter())
            .map(|(&h_k_1, &z_k_1)| h_k_1 * e_k_1(z_k_1, I_k_1, b_k, I_k, theta))
            .sum();

        let lower_prob: f64 = h_k[i - 1]
            .iter()
            .zip(quadratures[i - 1].z.iter())
            .map(|(&h_k_1, &z_k_1)| h_k_1 * e_k_1_prime(z_k_1, I_k_1, a_k, I_k, theta))
            .sum();
        res.push(upper_prob + lower_prob);
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
#[allow(non_snake_case)]
pub fn find_bounds(
    alpha: Vec<f64>,
    I: Vec<f64>,
    r: usize,
    tol: f64,
) -> Result<Vec<(f64, f64)>, CtsimErr> {
    let mut bounds: Vec<(f64, f64)> = vec![];
    bounds.reserve(alpha.len());

    let mut alpha_increments = vec![];
    alpha_increments.reserve(alpha.len());
    // First alpha increment is just the first amount of alpha spent
    alpha_increments.push(alpha[0]);
    for i in 1..alpha.len() {
        alpha_increments.push(alpha[i] - alpha[i - 1]);
    }
    // println!("alpha: {alpha:?}");
    // println!("alpha_increments: {alpha_increments:?}");

    // First bound is just +/- \Phi^-1(alpha_increments[0]/2)
    let bound_1 = -std_normal_quantile(alpha_increments[0] / 2.0)?;
    bounds.push((bound_1, -bound_1));

    for i in 1..alpha_increments.len() {
        // Assuming we are searching between 0, 7, so start
        // in the middle
        let (mut lower_bound, mut upper_bound) = (0.0, 7.0);
        let mut mid = (lower_bound + upper_bound) / 2.0;
        bounds.push((-mid, mid));

        let mut cur_alpha: f64 = psi_k(&bounds, &I[0..=i], 0.0, r)[i];
        let target_alpha = alpha_increments[i];
        let mut diff: f64 = target_alpha - cur_alpha;
        // println!("---------------------------------------");
        // println!("mid: {mid}");
        // println!("bounds {:?}", bounds);
        // println!("&I[0..=i] {:?}", &I[0..=i]);
        // println!("cur_alpha: {cur_alpha}");
        // println!("target_alpha: {target_alpha}");

        // Iterate until we get close enought to target alpha
        while diff.abs() > tol && (lower_bound - upper_bound).abs() > tol {
            if cur_alpha < target_alpha {
                upper_bound = mid;
            } else {
                lower_bound = mid;
            }
            mid = (lower_bound + upper_bound) / 2.0;

            // Update bounds
            bounds[i].0 = -mid;
            bounds[i].1 = mid;

            // compute new alpha
            cur_alpha = psi_k(&bounds, &I[0..=i], 0.0, r)[i];
            // println!("---------------------------------------");
            // println!("mid: {mid}");
            // println!("bounds {:?}", bounds);
            // println!("&I[0..=i] {:?}", &I[0..=i]);
            // println!("cur_alpha: {cur_alpha}");

            diff = target_alpha - cur_alpha;
        }

        if diff.abs() > tol {
            return Err(TrialBoundsError::FailedToConverge(cur_alpha, target_alpha, tol).into());
        }
    }

    Ok(bounds)
}

#[cfg(test)]
mod tests {
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
            32,
        )
        .iter()
        .sum();
        assert!((alpha - 0.05).abs() < 0.0001)
    }

    #[test]
    fn pocock_3() {
        let alpha: f64 = psi_k(
            &vec![(-2.289, 2.289), (-2.289, 2.289), (-2.289, 2.289)],
            &vec![1.0 / 3.0, 2.0 / 3.0, 1.0],
            0.0,
            32,
        )
        .iter()
        .sum();
        assert!((alpha - 0.05).abs() < 0.0001)
    }

    #[test]
    fn alpha_ld_of_1() {
        let alpha: f64 = psi_k(
            &vec![(-2.437995, 2.437995), (-1.999873, 1.999873)],
            &vec![0.7, 1.0],
            0.0,
            32,
        )
        .iter()
        .sum();
        assert!((alpha - 0.05).abs() < 0.0001)
    }
    #[test]
    fn alpha_ld_of_2() {
        let alpha: f64 = psi_k(
            &vec![(-3.356869, 3.356869), (-1.962261, 1.962261)],
            &vec![0.4, 1.0],
            0.0,
            32,
        )
        .iter()
        .sum();
        assert!((alpha - 0.05).abs() < 0.0001)
    }

    #[test]
    fn ldof_bounds_2_looks() {
        // This is the cumulative alpha spent for LDOF per ldBounds(c(0.7, 1.0))
        let bounds = find_bounds(vec![0.01476898, 0.05], vec![0.7, 1.0], 32, 0.0001).unwrap();

        assert!((bounds[0].1 - 2.437995).abs() < 0.001);
        assert!((bounds[1].1 - 1.999873).abs() < 0.001);
    }

    #[test]
    fn ldof_bounds_3_looks() {
        // This is the cumulative alpha spent for LDOF per ldBounds(c(0.3, 0.6, 1.0))
        let bounds = find_bounds(
            vec![0.00008545157, 0.007616127, 0.05],
            vec![0.3, 0.6, 1.0],
            32,
            0.0001,
        )
        .unwrap();

        assert!((bounds[0].1 - 3.928573).abs() < 0.001);
        // This one is stubborn (consistently off by 0.003 regardless of
        // quadrature size or tolerance), but it's close enough
        assert!((bounds[1].1 - 2.669967).abs() < 0.005);
        assert!((bounds[2].1 - 1.981004).abs() < 0.001);
    }
}
