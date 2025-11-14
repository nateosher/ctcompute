use rand::{
    SeedableRng,
    distributions::{Distribution, WeightedIndex},
    rngs,
};
use statrs::distribution::Uniform;

pub struct EnrollmentSimResult {
    pub enrollment_times: Vec<f64>,
    pub expected_time: f64,
}

pub fn sim_enrollment_times(
    n: usize,
    enrollment_times: &Vec<f64>,
    enrollment_rates: &Vec<f64>,
    seed: u64,
) -> EnrollmentSimResult {
    //----------------------------------------
    // Set up enrollment rates/times
    let mut time_lengths: Vec<f64> = enrollment_times.windows(2).map(|w| w[1] - w[0]).collect();
    let total_expected_enrollment: f64 = enrollment_rates[0..enrollment_rates.len() - 1]
        .iter()
        .zip(time_lengths.iter())
        .map(|(lambda, delta)| lambda * delta)
        .sum();
    let mut additional_time_needed = 0.0;
    if (total_expected_enrollment.floor() as usize) < n {
        let additional_enrollment_needed = (n as f64) - total_expected_enrollment;
        additional_time_needed = additional_enrollment_needed
            / enrollment_rates
                .last()
                .expect("enrollment rates should not be empty");
    }
    time_lengths.push(additional_time_needed);
    let total_time: f64 = time_lengths.iter().sum();
    let time_proportions: Vec<f64> = time_lengths.iter().map(|l| l / total_time).collect();

    //----------------------------------------
    // Set up rngs/distributions
    let mut master_rng = rngs::StdRng::seed_from_u64(seed);
    let interval_distribution = WeightedIndex::new(&time_proportions).unwrap();
    let mut random_interval_uniform_distributions: Vec<_> = time_lengths
        .iter()
        .enumerate()
        .map(|(i, l)| {
            let unif_rng = rngs::StdRng::seed_from_u64(seed + (i as u64));
            Uniform::new(0.0, *l).unwrap().sample_iter(unif_rng)
        })
        .collect();

    //----------------------------------------
    // Simulate enrollment
    let mut patient_enrollment_times: Vec<f64> = (0..n)
        .map(|_| {
            let enrollment_interval_index = interval_distribution.sample(&mut master_rng) as usize;
            let enrollment_interval_entry = random_interval_uniform_distributions
                [enrollment_interval_index]
                .next()
                .unwrap();
            enrollment_times[enrollment_interval_index] + enrollment_interval_entry
        })
        .collect();

    patient_enrollment_times.sort_by(|a, b| a.partial_cmp(b).unwrap());

    EnrollmentSimResult {
        enrollment_times: patient_enrollment_times,
        expected_time: time_lengths.iter().sum(),
    }
}
