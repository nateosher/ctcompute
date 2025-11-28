//----------------------------------------
// spending mod types
//----------------------------------------
pub enum SpendingFcn {
    LDOF, // Lan-Demets Obrien Fleming
}

#[derive(Debug)]
pub enum AlphaSpendingValues {
    OneSidedUpper(Vec<f64>),
    OneSidedLower(Vec<f64>),
    TwoSided((Vec<f64>, Vec<f64>)),
}
