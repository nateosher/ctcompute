use crate::spending_fcns::SpendingFcnErr;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum CtsimErr {
    #[error("while evaluating spending function: {0}")]
    SpendingFcn(SpendingFcnErr),
}
