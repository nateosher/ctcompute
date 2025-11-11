use itertools::{Itertools, izip};
use rand::distributions::DistIter;
use rand::{SeedableRng, distributions::Distribution, rngs};
use statrs::distribution::{ContinuousCDF, DiscreteUniform, Exp, Normal};
use statrs::statistics::{Data, OrderStatistics, Statistics};

use crate::enrollment_sim::sim_enrollment_times;

pub fn run_n_tte_sims(
    n: usize,
    sample_size_trt: usize,
    sample_size_ctrl: usize,
    lambda_trt: f64,
    lambda_ctrl: f64,
    lambda_trt_dropout: f64,
    lambda_ctrl_dropout: f64,
    enrollment_times: &Vec<f64>,
    enrollment_rates: &Vec<f64>,
    seed: u64,
) -> Vec<f64> {
    (0..n)
        .map(|i| {
            run_tte_sim(
                sample_size_trt,
                sample_size_ctrl,
                lambda_trt,
                lambda_ctrl,
                lambda_trt_dropout,
                lambda_ctrl_dropout,
                enrollment_times,
                enrollment_rates,
                seed + (i as u64),
            )
        })
        .collect()
}

pub fn run_tte_sim(
    sample_size_trt: usize,
    sample_size_ctrl: usize,
    lambda_trt: f64,
    lambda_ctrl: f64,
    lambda_trt_dropout: f64,
    lambda_ctrl_dropout: f64,
    enrollment_times: &Vec<f64>,
    enrollment_rates: &Vec<f64>,
    seed: u64,
) -> f64 {
    //----------------------------------------
    // Choose seeds based on given seed
    // TODO: simplify?
    let master_rng = rngs::StdRng::seed_from_u64(seed);
    let seed_distribution = DiscreteUniform::new(1000000, i64::MAX).unwrap();
    let mut seed_generator: DistIter<_, _, i64> = seed_distribution.sample_iter(master_rng);
    let mut seed_1 = seed_generator.next().unwrap() as u64;
    let mut seed_2 = seed_generator.next().unwrap() as u64;
    let mut seed_3 = seed_generator.next().unwrap() as u64;
    let mut seed_4 = seed_generator.next().unwrap() as u64;

    //----------------------------------------
    // Enrollment times
    let patient_enrollment_times = sim_enrollment_times(
        sample_size_trt + sample_size_ctrl, // n
        enrollment_times,                   // enrollment_times
        enrollment_rates,                   // enrollment_rates
        seed,                               // seed
    );

    //----------------------------------------
    // Treatment survival distribution
    let trt_surv_exp = Exp::new(lambda_trt).unwrap();
    let trt_surv_rng = rngs::StdRng::seed_from_u64(seed_1); // rand::thread_rng();
    let trt_surv_samples = trt_surv_exp
        .sample_iter(trt_surv_rng)
        .take(sample_size_trt)
        .zip(patient_enrollment_times.iter())
        .map(|(t_surv, t_enroll)| t_surv + t_enroll);

    //----------------------------------------
    // Control survival distribution
    let ctrl_surv_exp = Exp::new(lambda_ctrl).unwrap();
    let ctrl_surv_rng = rngs::StdRng::seed_from_u64(seed_2);
    let ctrl_surv_samples = ctrl_surv_exp
        .sample_iter(ctrl_surv_rng)
        .take(sample_size_ctrl)
        .zip(patient_enrollment_times.iter())
        .map(|(t_surv, t_enroll)| t_surv + t_enroll);

    //----------------------------------------
    // Treatment dropout distribution
    let trt_drop_exp = Exp::new(lambda_trt_dropout).unwrap();
    let trt_drop_rng = rngs::StdRng::seed_from_u64(seed_3);
    let trt_drop_samples = trt_drop_exp
        .sample_iter(trt_drop_rng)
        .take(sample_size_trt)
        .zip(patient_enrollment_times.iter())
        .map(|(t_drop, t_enroll)| t_drop + t_enroll);

    //----------------------------------------
    // Control dropout distribution
    let ctrl_drop_exp = Exp::new(lambda_ctrl_dropout).unwrap();
    let ctrl_drop_rng = rngs::StdRng::seed_from_u64(seed_4);
    let ctrl_drop_samples = ctrl_drop_exp
        .sample_iter(ctrl_drop_rng)
        .take(sample_size_ctrl)
        .zip(patient_enrollment_times.iter())
        .map(|(t_drop, t_enroll)| t_drop + t_enroll);

    //----------------------------------------
    // Simulate actual data
    // 1 for treatment, 0 for control
    let arms =
        std::iter::repeat_n(1, sample_size_trt).chain(std::iter::repeat_n(0, sample_size_ctrl));
    let surv_times = trt_surv_samples.chain(ctrl_surv_samples);
    let drop_times = trt_drop_samples.chain(ctrl_drop_samples);
    let raw_data = izip!(arms, surv_times, drop_times);
    let times_and_events: Vec<(usize, f64, usize)> = raw_data
        .map(|(arm, t_surv, t_drp)| {
            let t = t_surv.min(t_drp);
            let delta = if t_surv < t_drp { 1 } else { 0 };
            (arm, t, delta)
        })
        .sorted_by(|a, b| {
            f64::partial_cmp(&a.1, &b.1).expect("attempted to sort survival data with NaNs")
        })
        .collect();

    //----------------------------------------
    // Compute number of events at each time
    let initial_value = (sample_size_trt, sample_size_ctrl);
    let n_at_each_time = std::iter::once(initial_value).chain(times_and_events.iter().scan(
        initial_value,
        |current_count, (arm, _, _)| {
            current_count.0 -= arm;
            current_count.1 -= 1 - arm;
            if *current_count == (0_usize, 0_usize) {
                return None;
            }
            Some(*current_count)
        },
    ));

    #[cfg(debug_assertions)]
    assert_eq!(
        n_at_each_time.clone().collect::<Vec<_>>().len(),
        times_and_events.len()
    );

    //----------------------------------------
    // Computing the actual logrank statistic
    let logrank = times_and_events.iter().zip(n_at_each_time).fold(
        (0.0, 0.0),
        |acc, ((arm, _, delta), (n_treat, n_ctrl))| {
            if *delta == 0 {
                return acc;
            }
            #[allow(non_snake_case)]
            let E = (n_treat as f64) / (n_treat as f64 + n_ctrl as f64);
            #[allow(non_snake_case)]
            let O = if *arm == 1 { 1.0 } else { 0.0 };
            #[allow(non_snake_case)]
            let V = E * (1.0 - E);

            (acc.0 + (O - E), acc.1 + V)
        },
    );
    logrank.0 / logrank.1.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_sim_test() {
        run_tte_sim(
            50,                    // sample_size_trt
            50,                    // sample_size_ctrl
            0.2,                   // lambda_trt
            0.1,                   // lambda_ctrl
            113.9,                 // lambda_trt_dropout
            113.9,                 // lambda_ctrl_dropout
            &vec![1.0, 3.0, 6.0],  // enrollment times
            &vec![1.0, 5.0, 11.0], //enrollment rates
            24601,                 // seed
        );
    }
}
