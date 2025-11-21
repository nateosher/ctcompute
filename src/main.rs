mod arm;
mod computation_target;
mod ctsim_err;
mod delta;
mod enrollment_sim;
mod hypothesis_type;
mod integrate;
mod quadrature;
mod spending_fcns;
mod surv_sim_settings;
mod tte_compute;
mod tte_sim;

use enrollment_sim::{EnrollmentSimResult, sim_enrollment_times};
use std::time::Instant;
use tte_sim::{run_n_tte_sims, run_tte_sim};

fn main() {
    let test_enrollment_times = vec![0.0, 3.0, 7.0, 12.0];
    let test_enrollment_rates = vec![1.0, 3.0, 7.0, 12.0];
    let mut start = Instant::now();
    let single_tte_result = run_tte_sim(
        50,                     // sample_size_trt
        50,                     // sample_size_ctrl
        0.2,                    // lambda_trt
        0.1,                    // lambda_ctrl
        0.00877,                // lambda_trt_dropout
        0.00877,                // lambda_ctrl_dropout
        &test_enrollment_times, // enrollment_times
        &test_enrollment_rates, // enrollment_rates
        24601,                  // seed
    );
    let mut duration = start.elapsed();
    println!("Single tte sim (n = 100): {:?}", duration);
    println!("Single tte sim result: {:#?}", single_tte_result);

    start = Instant::now();
    run_n_tte_sims(
        10000,                  // n (# simulations)
        50,                     // sample_size_trt
        50,                     // sample_size_ctrl
        0.2,                    // lambda_trt
        0.1,                    // lambda_ctrl
        113.9,                  // lambda_trt_dropout
        113.9,                  // lambda_ctrl_dropout
        &test_enrollment_times, // enrollment_times
        &test_enrollment_rates, // enrollment_rates
        24601,                  // seed
    );
    duration = start.elapsed();
    println!("10k tte sims (n = 100): {:?}", duration);

    let patient_enrollment_sim = sim_enrollment_times(
        20,                         // n
        &vec![0.0, 3.0, 7.0, 12.0], // enrollment_times
        &vec![1.0, 1.0, 1.0, 1.0],  // enrollment_rates
        2460123456,                 // seed
    );

    println!(
        "Patient enrollment times: {:?}",
        patient_enrollment_sim.enrollment_times
    );
}
