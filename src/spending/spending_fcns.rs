use crate::error::CtcomputeErr;
use crate::information::std_normal::{std_normal_cdf, std_normal_quantile};
use crate::spending::{
    error::SpendingFcnErr,
    types::{AlphaSpendingValues, SpendingFcn},
};

// Returns cumulative alpha spent at each look
// TODO: ensure spending is monotonically increasing, jumps
// are > 0, and final value = alpha
pub fn compute_spending_vec(
    look_fractions: &Vec<f64>,
    alpha: f64,
    maybe_lower_spending_fcn_type: Option<&SpendingFcn>,
    maybe_upper_spending_fcn_type: Option<&SpendingFcn>,
) -> Result<AlphaSpendingValues, CtcomputeErr> {
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
    let maybe_lower_spending_values = match maybe_lower_spending_fcn_type {
        Some(SpendingFcn::LDOF) => look_fractions
            .iter()
            .map(|&t| lan_demets_obrien_fleming(t, alpha))
            .collect::<Result<Vec<f64>, CtcomputeErr>>()
            .map(Some)?,
        Some(SpendingFcn::Custom { cumulative_spend }) => {
            if cumulative_spend.len() != look_fractions.len() {
                return Err(SpendingFcnErr::MismatchedLengths.into());
            }
            if cumulative_spend[cumulative_spend.len() - 1] != alpha {
                return Err(SpendingFcnErr::BadLastSpend(
                    cumulative_spend[cumulative_spend.len() - 1],
                )
                .into());
            }
            Ok(Some(cumulative_spend.clone()))
        }?,
        None => None,
    };

    let maybe_upper_spending_values = match maybe_upper_spending_fcn_type {
        Some(SpendingFcn::LDOF) => look_fractions
            .iter()
            .map(|&t| lan_demets_obrien_fleming(t, alpha))
            .collect::<Result<Vec<f64>, CtcomputeErr>>()
            .map(Some)?,
        Some(SpendingFcn::Custom { cumulative_spend }) => {
            if cumulative_spend.len() != look_fractions.len() {
                return Err(SpendingFcnErr::MismatchedLengths.into());
            }
            if cumulative_spend[cumulative_spend.len() - 1] != alpha {
                return Err(SpendingFcnErr::BadLastSpend(
                    cumulative_spend[cumulative_spend.len() - 1],
                )
                .into());
            }
            Ok(Some(cumulative_spend.clone()))
        }?,
        None => None,
    };

    match (maybe_lower_spending_values, maybe_upper_spending_values) {
        (Some(lower), Some(upper)) => Ok(AlphaSpendingValues::TwoSided((lower, upper))),
        (None, Some(upper)) => Ok(AlphaSpendingValues::OneSidedUpper(upper)),

        (Some(lower), None) => Ok(AlphaSpendingValues::OneSidedLower(lower)),
        _ => unreachable!(),
    }
}

fn lan_demets_obrien_fleming(t: f64, alpha: f64) -> Result<f64, CtcomputeErr> {
    if t < 0.0 || t > 1.0 {
        Err(SpendingFcnErr::OutOfBounds(t).into())
    } else if t == 1.0 {
        // Hardcode to avoid numerical precision issues
        Ok(alpha)
    } else {
        // let std_normal = Normal::new(0.0, 1.0).unwrap();
        let z_alpha = std_normal_quantile(1. - alpha / 2.).expect("failed to compute z_alpha");
        let spend = (2. - 2. * std_normal_cdf(z_alpha / t.sqrt())).min(alpha);
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
            Some(&SpendingFcn::LDOF),
            Some(&SpendingFcn::LDOF),
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
            Some(&SpendingFcn::LDOF),
            Some(&SpendingFcn::LDOF),
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

    #[test]
    fn ldof_3_look_0025_2() {
        let alpha_spend = compute_spending_vec(
            &vec![1. / 3., 2. / 3., 1.0],
            0.025,
            Some(&SpendingFcn::LDOF),
            None,
        )
        .unwrap();

        if let AlphaSpendingValues::OneSidedLower(lower) = alpha_spend {
            assert!((lower[0] - 0.0001035057).abs() < 0.0001);
            assert!((lower[1] - lower[0] - 0.0059448834).abs() < 0.0001);
            assert!((lower[2] - lower[1] - lower[0] - 0.0189516109).abs() < 0.001);
        } else {
            panic!()
        }
    }

    #[test]
    fn custom_1() {
        let alpha_spend_custom = compute_spending_vec(
            &vec![0.3, 0.6, 1.0],
            0.025,
            Some(&SpendingFcn::Custom {
                cumulative_spend: vec![0.01, 0.02, 0.025],
            }),
            Some(&SpendingFcn::Custom {
                cumulative_spend: vec![0.1, 0.2, 0.025],
            }),
        );

        if let Ok(AlphaSpendingValues::TwoSided((lower_spend, upper_spend))) = alpha_spend_custom {
            assert_eq!(lower_spend, vec![0.01, 0.02, 0.025]);
            assert_eq!(upper_spend, vec![0.1, 0.2, 0.025]);
        } else {
            panic!()
        }
    }
}
