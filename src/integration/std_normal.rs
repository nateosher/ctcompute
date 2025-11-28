use crate::error::CtsimErr;
use crate::integration::error::NormalDistErr;
use probability_rs::{Distribution, dist::normal::Normal};

pub fn std_normal_pdf(z: f64) -> f64 {
    (-z * z / 2.0).exp() / (2.0 * std::f64::consts::PI).sqrt()
}

pub fn std_normal_cdf(z: f64) -> f64 {
    let std_normal = Normal::new(0.0, 1.0).unwrap();
    std_normal.cdf(z)
}

// from probability_rs; my temporary workaround to quantile bug
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
    Ok(std_normal_quantile_helper(p))
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
}
