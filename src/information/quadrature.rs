// z are locations, w are weights
#[derive(Debug)]
pub struct Quadrature {
    pub z: Vec<f64>,
    pub w: Vec<f64>,
}

// Jennison and Turnbull, 356
// TODO: error handling
impl Quadrature {
    #[allow(non_snake_case)]
    pub fn new(theta: f64, I_j: f64, r: usize, a: f64, b: f64) -> Quadrature {
        let theta_sqrt_I_j = theta * I_j.sqrt();
        let r_f = r as f64;

        // J&T, bottom of page 356
        let x_full: Vec<f64> = (1..6 * r)
            .map(|i| {
                let i_f = i as f64;
                if 1 <= i && i <= r - 1 {
                    theta_sqrt_I_j + (-3.0 - 4.0 * ((r_f / i_f).ln()))
                } else if r <= i && i <= 5 * r {
                    theta_sqrt_I_j + (-3.0 + 3.0 * (i_f - r_f) / (2.0 * r_f))
                } else if 5 * r + 1 <= i && i <= 6 * r - 1 {
                    theta_sqrt_I_j + (3.0 + 4.0 * (r_f / (6.0 * r_f - i_f)).ln())
                } else {
                    unreachable!()
                }
            })
            .collect();

        // Trim to be within (a, b)
        let mut x: Vec<f64> = x_full
            .iter()
            .filter(|&&i| (i as f64) > a && (i as f64) < b)
            .map(|&i| i)
            .collect();

        // If there were no points placed between a and b, there is virtually
        // no mass between (a, b), so just put a single point with weight 1
        // between a and b
        if x.is_empty() {
            return Quadrature {
                z: vec![(a + b) / 2.],
                w: vec![1.0],
            };
        }

        // If lower end was less than a, make first element a
        if x[0] != x_full[0] {
            x.insert(0, a);
        }
        // If upper end was greater than b, make last element b
        if x[x.len() - 1] != x_full[x_full.len() - 1] {
            x.push(b);
        }

        //----------------------------------------
        // Z values, i.e. locations
        // J&T, top of page 357

        // This is equivalent to 12r - 3 when no points are trimmed
        let m = 2 * x.len() - 1;
        let mut z: Vec<f64> = vec![f64::INFINITY; m];

        // Odds first
        for i in (0..m).step_by(2) {
            // Textbook indexes from 1, so have to adjust
            let x_index = i / 2;
            z[i] = x[x_index];
        }

        // Then evens
        for i in (1..m).step_by(2) {
            z[i] = (z[i - 1] + z[i + 1]) / 2.0;
        }

        //----------------------------------------
        // Weights
        let mut w: Vec<f64> = vec![f64::INFINITY; m];
        let one_sixth = 1.0 / 6.0;
        let four_sixths = 2.0 / 3.0;

        // TODO: re-write, fix indices?
        for i in 0..m {
            if i == 0 {
                w[i] = one_sixth * (z[2] - z[0]);
            } else if i == m - 1 {
                w[i] = one_sixth * (z[m - 1] - z[m - 3]);
            } else if i % 2 == 0 {
                // Odd/even indices are switched, since book indexes
                // from one instead of 0
                w[i] = one_sixth * (z[i + 2] - z[i - 2]);
            } else if i % 2 == 1 {
                // See comment above, same deal
                w[i] = four_sixths * (z[i + 1] - z[i - 1]);
            } else {
                unreachable!()
            }
        }

        Quadrature { z, w }
    }
}

#[cfg(test)]
mod tests {
    use std::f64;

    use super::*;
    use crate::information::std_normal::std_normal_pdf;

    #[test]
    fn basic_quadrature_indefinite() {
        let test_quad = Quadrature::new(0.0, 1.0, 4, f64::NEG_INFINITY, f64::INFINITY);
        assert!(test_quad.w.iter().all(|w_j| w_j.is_finite()));
        assert!(test_quad.z.iter().all(|z_j| z_j.is_finite()));
    }

    #[test]
    fn basic_quadrature_trimmed() {
        let test_quad = Quadrature::new(0.0, 1.0, 4, -1.0, 1.0);
        assert!(test_quad.z.last().unwrap() <= &1.0);
        assert!(test_quad.z[0] >= -1.0);
    }

    #[test]
    fn basic_quadrature_indefinite_2() {
        let test_quad = Quadrature::new(0.0, 1.0, 5, f64::NEG_INFINITY, f64::INFINITY);
        assert!(test_quad.w.iter().all(|w_j| w_j.is_finite()));
        assert!(test_quad.z.iter().all(|z_j| z_j.is_finite()));
    }

    #[test]
    fn basic_quadrature_trimmed_2() {
        let test_quad = Quadrature::new(0.0, 1.0, 5, -1.0, 1.0);
        assert!(test_quad.z.last().unwrap() <= &1.0);
        assert!(test_quad.z[0] >= -1.0);
    }

    #[test]
    fn standard_normal_integral() {
        let test_quad = Quadrature::new(0.0, 1.0, 16, -1.959964, 1.959964);
        let std_normal_integral: f64 = test_quad
            .z
            .iter()
            .zip(test_quad.w.iter())
            .map(|(z, w)| w * std_normal_pdf(*z))
            .sum();
        assert!((std_normal_integral - 0.95).abs() < 0.00000001);
    }

    #[test]
    fn standard_normal_integral_2() {
        let test_quad = Quadrature::new(0.0, 1.0, 16, -1.5, 2.7);
        let std_normal_integral: f64 = test_quad
            .z
            .iter()
            .zip(test_quad.w.iter())
            .map(|(z, w)| w * std_normal_pdf(*z))
            .sum();
        assert!((std_normal_integral - 0.9297258).abs() < 0.0000001);
    }
}
