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
    #[error("spending function vector should end with 1.0; got {0}")]
    BadLastSpend(f64),
    #[error("time vector was empty")]
    TimeVectorEmpty,
    #[error("total alpha spent should be in (0, 1); got {0}")]
    BadAlpha(f64),
}

impl Into<CtsimErr> for SpendingFcnErr {
    fn into(self) -> CtsimErr {
        CtsimErr::SpendingFcn(self)
    }
}

fn check_spending_fcn_arguments(t_v: &Vec<f64>, alpha: f64) -> Result<(), CtsimErr> {
    if t_v.is_empty() {
        return Err(SpendingFcnErr::TimeVectorEmpty.into());
    }
    if let Some(&last_spend) = t_v.last()
        && last_spend != 1.0
    {
        return Err(SpendingFcnErr::BadLastSpend(last_spend).into());
    }
    if alpha <= 0.0 || alpha >= 1.0 {
        return Err(SpendingFcnErr::BadAlpha(alpha).into());
    }
    Ok(())
}

// Returns cumulative alpha spent at each look
// Note to self: gsDesign splits alpha before applying spending
// function; keep this in mind when comparing
// TODO: refactor into a general spending vector computation function
// that takes spending function as argument (probably through enum)
pub fn lan_demets_obrien_fleming_vec(t_v: &Vec<f64>, alpha: f64) -> Result<Vec<f64>, CtsimErr> {
    check_spending_fcn_arguments(t_v, alpha)?;
    let timepoint_spend: Vec<f64> = t_v
        .iter()
        .map(|&t| lan_demets_obrien_fleming(t, alpha))
        .collect::<Result<Vec<f64>, CtsimErr>>()?;

    let mut cumulative_spend = timepoint_spend.clone();
    // First and last timepoint spends don't need to be adjusted
    for i in 1..cumulative_spend.len() - 1 {
        println!("cumulative spend updated");
        cumulative_spend[i] += cumulative_spend[i - 1];
    }

    Ok(cumulative_spend)
}

fn lan_demets_obrien_fleming(t: f64, alpha: f64) -> Result<f64, CtsimErr> {
    if t < 0.0 || t > 1.0 {
        Err(SpendingFcnErr::OutOfBounds(t).into())
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
                    "while evaluating spending function: arguments to \
                    spending function should be in [0, 1]; got 1.1"
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
