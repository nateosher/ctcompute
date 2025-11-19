use crate::ctsim_err::CtsimErr;
use probability_rs::{Continuous, Distribution, dist::normal::Normal};
use thiserror::Error;

pub enum SpendingFcn {
    LDOF,
}

#[derive(Error, Debug)]
pub enum SpendingFcnErr {
    #[error("arguments to spending function should be in [0, 1]; got {0}")]
    OutOfBounds(f64),
}

fn lan_demets_obrien_fleming_vec(v: &Vec<f64>, alpha: f64) -> Vec<Result<f64, CtsimErr>> {
    v.iter()
        .map(|&t| lan_demets_obrien_fleming(t, alpha))
        .collect()
}

fn lan_demets_obrien_fleming(t: f64, alpha: f64) -> Result<f64, CtsimErr> {
    if t < 0.0 || t > 1.0 {
        Err(CtsimErr::SpendingFcn(SpendingFcnErr::OutOfBounds(t)))
    } else {
        let std_normal = Normal::new(0.0, 1.0).unwrap();
        let z_alpha_d_2 = std_normal.inv_cdf(1.0 - alpha / 2.0);
        let spend = (2.0 - 2.0 * std_normal.cdf(z_alpha_d_2 / t.sqrt())).min(alpha);
        Ok(spend)
    }
}

#[cfg(test)]
mod tests {
    use std::panic;

    use super::*;

    #[test]
    fn spending_fcn_error() {
        if let Err(e) = lan_demets_obrien_fleming(1.1, 0.05) {
            assert_eq!(
                String::from(
                    "while evaluating spending function: arguments to spending function should be in [0, 1]; got 1.1"
                ),
                format!("{}", e)
            );
        } else {
            panic!()
        }
    }

    #[test]
    fn ldof_0_75_0_05() {
        assert!(lan_demets_obrien_fleming(0.75, 0.05).is_ok_and(|x| (x - 0.0236).abs() < 0.0001))
    }
}
