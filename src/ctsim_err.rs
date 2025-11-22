use crate::{
    integrate::{NormalDistErr, TrialBoundsError},
    spending_fcns::SpendingFcnErr,
};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum CtsimErr {
    #[error("while evaluating spending function: {0}")]
    SpendingFcn(SpendingFcnErr),
    #[error("while evaluating normal distribution: {0}")]
    NormalDist(NormalDistErr),
    #[error("while computing trial bounds: {0}")]
    TrialBounds(TrialBoundsError),
}
