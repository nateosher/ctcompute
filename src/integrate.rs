use crate::quadrature::Quadrature;
use probability_rs::{Distribution, dist::normal::Normal};

pub fn std_normal_pdf(z: f64) -> f64 {
    (-z * z / 2.0).exp() / (2.0 * std::f64::consts::PI).sqrt()
}

pub fn std_normal_cdf(z: f64) -> f64 {
    let std_normal = Normal::new(0.0, 1.0).unwrap();
    std_normal.cdf(z)
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

pub struct IntegralResult {
    // cutoffs: Vec<f64>,
    // power: Vec<f64>,
    alpha: f64,
}

// J&T, p. 355
// TODO: rename, since this is technically psi_k and xi_k
// Computes exit probabilities given bounds, information fractions,
// and theta.
#[allow(non_snake_case)]
pub fn psi_k(bounds: Vec<(f64, f64)>, I: Vec<f64>, theta: f64, r: usize) -> Vec<f64> {
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

// Information fractions, alpha spend, and quadrature
// size (r), finds bounds that satisfy
#[allow(non_snake_case)]
pub fn find_bounds(alpha: Vec<f64>, I: Vec<f64>, r: usize) -> Vec<f64> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pocock_2() {
        let alpha: f64 = psi_k(
            vec![(-2.17878788, 2.17878788), (-2.17878788, 2.17878788)],
            vec![0.5, 1.0],
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
            vec![(-2.289, 2.289), (-2.289, 2.289), (-2.289, 2.289)],
            vec![1.0 / 3.0, 2.0 / 3.0, 1.0],
            0.0,
            32,
        )
        .iter()
        .sum();
        assert!((alpha - 0.05).abs() < 0.0001)
    }
    #[test]
    fn ld_of_1() {
        let alpha: f64 = psi_k(
            vec![(-2.437995, 2.437995), (-1.999873, 1.999873)],
            vec![0.7, 1.0],
            0.0,
            32,
        )
        .iter()
        .sum();
        println!("{alpha}");
        assert!((alpha - 0.05).abs() < 0.0001)
    }
    #[test]
    fn ld_of_2() {
        let alpha: f64 = psi_k(
            vec![(-3.356869, 3.356869), (-1.962261, 1.962261)],
            vec![0.4, 1.0],
            0.0,
            32,
        )
        .iter()
        .sum();
        println!("{alpha}");
        assert!((alpha - 0.05).abs() < 0.0001)
    }
}
