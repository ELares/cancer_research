/**
 * C integration test for ferroptosis-core FFI bindings.
 *
 * Compile and run:
 *   cargo build --release -p ferroptosis-ffi
 *   cc -o test_ffi tests/test_ffi.c -L ../target/release -lferroptosis_ffi -lm
 *   DYLD_LIBRARY_PATH=../target/release ./test_ffi     # macOS
 *   LD_LIBRARY_PATH=../target/release ./test_ffi        # Linux
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/ferroptosis.h"

int main(void) {
    printf("=== ferroptosis-core C FFI test ===\n\n");

    /* Create RNG */
    FerroRng *rng = ferro_rng_new(42);
    if (!rng) {
        fprintf(stderr, "ERROR: ferro_rng_new returned NULL\n");
        return 1;
    }

    /* Get default params */
    FerroParams params = ferro_params_default();
    printf("Params: death_threshold=%.1f, rsl3_gpx4_inhib=%.2f\n",
           params.death_threshold, params.rsl3_gpx4_inhib);

    /* Simulate 1000 Persister + RSL3 cells */
    int n_cells = 1000;
    int n_dead = 0;
    double sum_lp = 0.0;

    for (int i = 0; i < n_cells; i++) {
        FerroCell cell = ferro_gen_cell(2, rng);  /* 2 = Persister */
        FerroResult res = ferro_sim_cell(&cell, 1, &params, rng);  /* 1 = RSL3 */
        if (res.dead) n_dead++;
        sum_lp += res.final_lp;
    }

    double death_rate = (double)n_dead / (double)n_cells;
    double mean_lp = sum_lp / (double)n_cells;
    printf("\nPersister + RSL3 (n=%d):\n", n_cells);
    printf("  Death rate: %.1f%%\n", death_rate * 100.0);
    printf("  Mean LP: %.2f\n", mean_lp);

    /* Verify death rate is in expected range (30-55% for Persister+RSL3) */
    if (death_rate < 0.30 || death_rate > 0.55) {
        fprintf(stderr, "FAIL: Death rate %.1f%% outside expected range [30-55%%]\n",
                death_rate * 100.0);
        ferro_rng_free(rng);
        return 1;
    }

    /* Test invivo params */
    FerroParams invivo = ferro_params_invivo();
    printf("\nIn-vivo params: scd_mufa_rate=%.3f, initial_mufa=%.2f\n",
           invivo.scd_mufa_rate, invivo.initial_mufa_protection);

    /* Test NULL safety */
    FerroResult null_res = ferro_sim_cell(NULL, 0, &params, rng);
    if (null_res.dead) {
        fprintf(stderr, "FAIL: NULL cell should not return dead\n");
        ferro_rng_free(rng);
        return 1;
    }

    /* Cleanup */
    ferro_rng_free(rng);
    ferro_rng_free(NULL);  /* Should be a no-op */

    printf("\nAll tests PASSED.\n");
    return 0;
}
