//----------------------------------------
// spending mod types
//----------------------------------------
#[derive(Debug, Clone)]
pub enum SpendingFcn {
    LDOF, // Lan-Demets Obrien Fleming
    Custom { cumulative_spend: Vec<f64> },
}

#[derive(Debug)]
pub enum AlphaSpendingValues {
    OneSidedUpper(Vec<f64>),
    OneSidedLower(Vec<f64>),
    TwoSided((Vec<f64>, Vec<f64>)),
}
