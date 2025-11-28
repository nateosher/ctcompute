use crate::error::CtsimErr;
use crate::spending::{
    error::SpendingFcnErr,
    types::{AlphaSpendingValues, SpendingFcn},
};
use probability_rs::{Continuous, Distribution, dist::normal::Normal};

// Returns cumulative alpha spent at each look
// TODO: ensure spending is monotonically increasing, jumps
// are > 0, and final value = alpha
pub fn compute_spending_vec(
    look_fractions: &Vec<f64>,
    alpha: f64,
    maybe_lower_spending_fcn_type: Option<SpendingFcn>,
    maybe_upper_spending_fcn_type: Option<SpendingFcn>,
) -> Result<AlphaSpendingValues, CtsimErr> {
    //----------------------------------------
    // Check arguments
    if look_fractions.is_empty() {
        return Err(SpendingFcnErr::TimeVectorEmpty.into());
    }
    if let Some(&last_spend) = look_fractions.last()
        && last_spend != 1.0
    {
        return Err(SpendingFcnErr::BadLastSpend(last_spend).into());
    }
    if alpha <= 0.0 || alpha >= 1.0 {
        return Err(SpendingFcnErr::BadAlpha(alpha).into());
    }
    if maybe_lower_spending_fcn_type.is_none() && maybe_upper_spending_fcn_type.is_none() {
        return Err(SpendingFcnErr::NoSpendingFunctionSpecified.into());
    }

    //----------------------------------------
    // Compute alpha spend
    let maybe_lower_spending_fcn = match maybe_lower_spending_fcn_type {
        Some(SpendingFcn::LDOF) => Some(lan_demets_obrien_fleming),
        None => None,
    };

    let maybe_lower_spending_bound: Option<Vec<f64>> = match maybe_lower_spending_fcn {
        Some(spending_fcn) => look_fractions
            .iter()
            .map(|&t| spending_fcn(t, alpha))
            .collect::<Result<Vec<f64>, CtsimErr>>()
            .map(Some)?, // Turns into Result<Option<Vec<f64>>, CtsimErr>
        None => None,
    };

    let maybe_upper_spending_fcn = match maybe_upper_spending_fcn_type {
        Some(SpendingFcn::LDOF) => Some(lan_demets_obrien_fleming),
        None => None,
    };

    let maybe_upper_spending_bound: Option<Vec<f64>> = match maybe_upper_spending_fcn {
        Some(spending_fcn) => look_fractions
            .iter()
            .map(|&t| spending_fcn(t, alpha))
            .collect::<Result<Vec<f64>, CtsimErr>>()
            .map(Some)?, // Turns into Result<Option<Vec<f64>>, CtsimErr>
        None => None,
    };

    match (maybe_lower_spending_bound, maybe_upper_spending_bound) {
        (Some(lower), Some(upper)) => Ok(AlphaSpendingValues::TwoSided((lower, upper))),
        (None, Some(upper)) => Ok(AlphaSpendingValues::OneSidedUpper(upper)),

        (Some(lower), None) => Ok(AlphaSpendingValues::OneSidedLower(lower)),
        _ => unreachable!(),
    }
}

fn lan_demets_obrien_fleming(t: f64, alpha: f64) -> Result<f64, CtsimErr> {
    if t < 0.0 || t > 1.0 {
        Err(SpendingFcnErr::OutOfBounds(t).into())
    } else {
        let std_normal = Normal::new(0.0, 1.0).unwrap();
        let z_alpha_d_2 = std_normal.inv_cdf(alpha / 2.0);
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
        assert!(
            lan_demets_obrien_fleming(0.75, 0.025).is_ok_and(|x| (x - 0.009649325).abs() < 0.0001)
        )
    }

    #[test]
    fn ldof_2_look_005() {
        let alpha_spend = compute_spending_vec(
            &vec![0.7, 1.0],
            0.025,
            Some(SpendingFcn::LDOF),
            Some(SpendingFcn::LDOF),
        )
        .unwrap();

        if let AlphaSpendingValues::TwoSided((lower, upper)) = alpha_spend {
            assert!((lower[0] - 0.007384489).abs() < 0.0001);
            assert!((upper[0] - 0.007384489).abs() < 0.0001)
        } else {
            panic!()
        }
    }
    #[test]
    fn ldof_3_look_0025() {
        let alpha_spend = compute_spending_vec(
            &vec![0.3, 0.6, 1.0],
            0.025,
            Some(SpendingFcn::LDOF),
            Some(SpendingFcn::LDOF),
        )
        .unwrap();

        if let AlphaSpendingValues::TwoSided((lower, upper)) = alpha_spend {
            assert!((lower[0] - 4.272579e-05).abs() < 0.0001);
            assert!((upper[0] - 4.272579e-05).abs() < 0.0001);
            assert!((lower[1] - lower[0] - 3.765338e-03).abs() < 0.0001);
            assert!((upper[1] - lower[0] - 3.765338e-03).abs() < 0.0001);
        } else {
            panic!()
        }
    }
}
