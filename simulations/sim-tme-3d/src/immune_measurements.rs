//! Canonical 3D condition gate for the dimension-neutral passive ledger.

use crate::{Condition, DoseSchedule, Treatment};
pub use ferroptosis_core::immune_measurements::Measurements;

/// The public measured entry point supplies only default Overrides. Restrict
/// its condition as well, so this ledger cannot silently mislabel a different
/// geometry, death route, or runtime seed as the canonical comparison.
pub fn validate_condition(condition: &Condition) {
    let treatment_name = match condition.treatment {
        Treatment::Control => "Control",
        Treatment::RSL3 => "RSL3",
        Treatment::SDT => "SDT",
        _ => panic!("immune measurements support only canonical Control/RSL3/SDT arms"),
    };
    assert!(
        condition.name == format!("immune_{treatment_name}")
            && condition.treatment_name == treatment_name
            && condition.immune_on
            && !condition.stromal_on
            && !condition.ph_on
            && condition.o2_lambda == Some(crate::ZONE_REF_LAMBDA)
            && matches!(condition.dose_schedule, DoseSchedule::Constant),
        "immune measurements require the unchanged canonical immune condition"
    );
}
