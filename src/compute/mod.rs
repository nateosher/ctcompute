//----------------------------------------
// computation mod
//----------------------------------------
pub mod types;

pub use crate::enrollment::enrollment_sim::sim_enrollment_times;
pub use crate::integration::integrate::find_bounds;
pub use crate::sample_size::total_information::compute_information;
pub use crate::spending::spending_fcns::compute_spending_vec;
