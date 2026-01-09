//! Partial charge calculation for molecular systems.
//!
//! This module handles the assignment of partial atomic charges using
//! multiple methods, including global QEq, hybrid (biological + QEq), and
//! no charges.

mod hybrid;
mod qeq;
mod spatial;

use super::config::ChargeMethod;
use super::error::Error;
use super::intermediate::IntermediateSystem;

/// Assigns partial charges to atoms based on the configured method.
///
/// Modifies the `charge` field of each atom in the intermediate system
/// according to the specified charge method.
///
/// # Arguments
///
/// * `system` — Mutable reference to the intermediate system
/// * `method` — Charge calculation method to use
///
/// # Errors
///
/// Returns [`Error`] if:
/// - QEq solver fails to converge ([`Error::ChargeCalculation`])
/// - Hybrid method is used without biological metadata ([`Error::MissingBioMetadata`])
/// - Classical charge lookup fails ([`Error::HybridChargeAssignment`])
pub fn assign_charges(system: &mut IntermediateSystem, method: &ChargeMethod) -> Result<(), Error> {
    match method {
        ChargeMethod::None => Ok(()),
        ChargeMethod::Qeq(config) => qeq::assign_qeq_charges(system, config),
        ChargeMethod::Hybrid(config) => hybrid::assign_hybrid_charges(system, config),
    }
}
