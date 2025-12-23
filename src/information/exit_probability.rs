use crate::information::{
    quadrature::Quadrature,
    std_normal::{std_normal_cdf, std_normal_pdf},
    types::IntegralType,
};

// J&T p. 354
#[allow(non_snake_case)]
fn f_1(z_1: f64, I_1: f64, theta: f64) -> f64 {
    std_normal_pdf(z_1 - theta * I_1.sqrt())
}

// TODO: clean this up
#[allow(non_snake_case)]
fn f_k(z_k_1: f64, I_k_1: f64, z_k: f64, theta: f64, I_k: f64) -> f64 {
    #[cfg(debug_assertions)]
    assert!(I_k_1 > 0.);

    #[cfg(debug_assertions)]
    assert!(I_k > 0.);

    let delta_k = I_k - I_k_1;
    #[cfg(debug_assertions)]
    assert!(delta_k > 0.);

    let sqrt_delta_k = delta_k.sqrt();
    #[cfg(debug_assertions)]
    assert!(!sqrt_delta_k.is_nan());

    let prod_part_1 = I_k.sqrt() / sqrt_delta_k;
    let pass_1_1 = z_k * I_k.sqrt();
    #[cfg(debug_assertions)]
    if pass_1_1.is_infinite() {
        println!("z_k: {z_k}");
        println!("I_k.sqrt(): {}", I_k.sqrt());
        panic!("infinite pass_1_1");
    }

    let pass_1_2 = z_k_1 * I_k_1.sqrt();
    #[cfg(debug_assertions)]
    assert!(!pass_1_2.is_infinite());

    let pass_1_3 = theta * delta_k;
    #[cfg(debug_assertions)]
    assert!(!pass_1_3.is_infinite());

    let pass_1 = pass_1_1 - pass_1_2 - pass_1_3;
    #[cfg(debug_assertions)]
    if pass_1.is_nan() {
        println!("pass_1_1: {pass_1_1}");
        println!("pass_1_2: {pass_1_2}");
        println!("pass_1_3: {pass_1_3}");
        panic!("you already know")
    }
    if pass_1.is_nan() {
        return 0.;
    }

    let to_pass_to_std_normal_pdf = pass_1 / sqrt_delta_k;

    #[cfg(debug_assertions)]
    assert!(!to_pass_to_std_normal_pdf.is_nan());

    prod_part_1 * std_normal_pdf(to_pass_to_std_normal_pdf)
}

// e_k in J&T's notation
#[allow(non_snake_case)]
fn upper_exit_prob(z_k_1: f64, I_k_1: f64, b_k: f64, I_k: f64, theta: f64) -> f64 {
    let delta_k = I_k - I_k_1;
    std_normal_cdf((z_k_1 * I_k_1.sqrt() + theta * delta_k - b_k * I_k.sqrt()) / delta_k.sqrt())
}

// Defined to capture probability of being less than a; no
// corresponding definition in J&T
#[allow(non_snake_case)]
fn lower_exit_prob(z_k_1: f64, I_k_1: f64, a_k: f64, I_k: f64, theta: f64) -> f64 {
    let delta_k = I_k - I_k_1;
    std_normal_cdf((a_k * I_k.sqrt() - z_k_1 * I_k_1.sqrt() - theta * delta_k) / delta_k.sqrt())
}

// J&T, p. 355
// TODO: rename, since this is technically psi_k and xi_k
// Computes exit probabilities given bounds, information fractions,
// theta, and integral type (upper or lower). r controls quadrature size.
#[allow(non_snake_case)]
pub fn exit_probability(
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exit_prob_1() {
        let lower_prob = exit_probability(
            &vec![
                (-2.437990437182772, f64::INFINITY),
                (-2.1328125, f64::INFINITY),
                (-2.078125, f64::INFINITY),
            ],
            &vec![0.7, 0.9, 1.0],
            -34.66,
            IntegralType::Lower,
            32,
        );

        println!("{lower_prob:?}");
    }
}
