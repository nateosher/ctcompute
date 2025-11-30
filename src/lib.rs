//----------------------------------------
// Root lib
//----------------------------------------
//! The purpose of this library is to provide utility functions for computing
//! sample sizes, type I error, power, and stopping boundaries for clinical
//! trials of various kinds. Right now the focus is on group sequential
//! designs.

/// This module houses the public API for computing boundaries, spending
/// function values, and enrollment times
pub mod compute;
mod enrollment;
/// This module contains error types
pub mod error;
mod integration;
mod sample_size;
mod spending;
