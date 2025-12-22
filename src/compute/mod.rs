//----------------------------------------
// computation mod
//----------------------------------------
pub mod types;

pub use crate::duration::expected_enrollment::expected_enrollment_dur;
pub use crate::information::compute_information;
pub use crate::information::integrate::find_bounds;
pub use crate::sample_size::compute_ss::compute_ss_range;
pub use crate::spending::spending_fcns::compute_spending_vec;
