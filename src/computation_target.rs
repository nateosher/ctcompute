#[derive(Default, Debug, PartialEq, Clone, Copy)]
pub enum ComputationTarget {
    #[default]
    SampleSize,
    Alpha,
    Beta,
    EffectSize,
}
