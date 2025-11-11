mod computation_target;
mod hypothesis_type;
mod surv_sim_settings;
mod tte_compute;
mod tte_sim;

use std::time::Instant;
use tte_sim::{run_n_tte_sims, run_tte_sim};

fn main() {
    let mut start = Instant::now();
    run_tte_sim(
        50,    // sample_size_trt
        50,    // sample_size_ctrl
        0.2,   // lambda_trt
        0.1,   // lambda_ctrl
        113.9, // lambda_trt_dropout
        113.9, // lambda_ctrl_dropout
        24601, // seed
    );
    let mut duration = start.elapsed();
    println!("Single tte sim (n = 100): {:?}", duration);

    start = Instant::now();
    run_n_tte_sims(
        10000, // n (# simulations)
        50,    // sample_size_trt
        50,    // sample_size_ctrl
        0.2,   // lambda_trt
        0.1,   // lambda_ctrl
        113.9, // lambda_trt_dropout
        113.9, // lambda_ctrl_dropout
        24601, // seed
    );
    duration = start.elapsed();
    println!("10k tte sims (n = 100): {:?}", duration);
}
