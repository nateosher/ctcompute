//----------------------------------------
// compute mod types
//----------------------------------------

#[derive(Debug, Clone, Copy)]
pub enum ComputationTarget {
    SampleSize,
    Alpha,
    Power,
    EffectSize,
}

#[derive(Debug, Clone, Copy)]
pub enum EndpointType {
    Continuous,
    Binary,
    TimeToEvent,
}

pub use crate::spending::types::AlphaSpendingValues;
