#ifndef FERROPTOSIS_CORE_H
#define FERROPTOSIS_CORE_H

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

/**
 * Opaque random number generator. Create with ferro_rng_new, free with ferro_rng_free.
 */
typedef struct FerroRng FerroRng;

/**
 * Simulation parameters (20 f64 + 1 u32 = 21 fields).
 * Use ferro_params_default() or ferro_params_invivo() to create.
 */
typedef struct FerroParams {
  double fenton_rate;
  double gsh_scav_efficiency;
  double gsh_km;
  double nrf2_gsh_rate;
  double lp_rate;
  double lp_propagation;
  double gpx4_rate;
  double fsp1_rate;
  double scd_mufa_rate;
  double scd_mufa_max;
  double initial_mufa_protection;
  double scd_mufa_decay;
  double gpx4_degradation_by_ros;
  double gpx4_nrf2_upregulation;
  double sdt_ros;
  double pdt_ros;
  double rsl3_gpx4_inhib;
  double gsh_max;
  double gpx4_nrf2_target_multiplier;
  double death_threshold;
  uint32_t post_death_steps;
} FerroParams;

/**
 * Cell biochemical parameters (7 fields, all f64).
 */
typedef struct FerroCell {
  double iron;
  double gsh;
  double gpx4;
  double fsp1;
  double basal_ros;
  double lipid_unsat;
  double nrf2;
} FerroCell;

/**
 * Result of a single-cell ferroptosis simulation.
 */
typedef struct FerroResult {
  bool dead;
  double final_lp;
  double final_gsh;
  double final_gpx4;
} FerroResult;

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/**
 * Create a new random number generator seeded with the given value.
 * The caller owns the returned pointer and MUST free it with ferro_rng_free.
 * Each thread should own its own FerroRng instance.
 */
struct FerroRng *ferro_rng_new(uint64_t seed);

/**
 * Free a FerroRng created by ferro_rng_new. Passing NULL is safe (no-op).
 * Do NOT double-free or use after free.
 */
void ferro_rng_free(struct FerroRng *rng);

/**
 * Create default simulation parameters (2D culture).
 * Returns a FerroParams struct by value (all fields populated).
 */
struct FerroParams ferro_params_default(void);

/**
 * Create in-vivo simulation parameters (3D/in-vivo with SCD1/MUFA protection).
 * Returns a FerroParams struct by value.
 */
struct FerroParams ferro_params_invivo(void);

/**
 * Generate a stochastic cell with phenotype-specific biochemical parameters.
 *
 * phenotype: 0=Glycolytic, 1=OXPHOS, 2=Persister, 3=PersisterNrf2, 4=Stromal.
 * Invalid values default to Glycolytic.
 *
 * rng: Must be a valid FerroRng pointer (from ferro_rng_new). Must not be NULL.
 */
struct FerroCell ferro_gen_cell(int32_t phenotype, struct FerroRng *rng);

/**
 * Run a full 180-step ferroptosis simulation for one cell.
 *
 * cell: Pointer to a FerroCell (from ferro_gen_cell). Must not be NULL.
 * treatment: 0=Control, 1=RSL3, 2=SDT, 3=PDT. Invalid values default to Control.
 * params: Pointer to FerroParams (from ferro_params_default/invivo). Must not be NULL.
 * rng: Must be a valid FerroRng pointer. Must not be NULL.
 *
 * Returns a FerroResult with dead status and final LP/GSH/GPX4 values.
 */
struct FerroResult ferro_sim_cell(const struct FerroCell *cell,
                                  int32_t treatment,
                                  const struct FerroParams *params,
                                  struct FerroRng *rng);

#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus

#endif  /* FERROPTOSIS_CORE_H */
