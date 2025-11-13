#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Delta {
    Observed,
    Censored,
}
