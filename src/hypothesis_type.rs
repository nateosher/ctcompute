#[derive(Default, Debug, PartialEq, Copy, Clone)]
pub enum HypothesisType {
    #[default]
    NotEqual,
    TrtGreater,
    TrtLess,
}
