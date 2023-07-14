/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include <cstdio>
#include "max2sat.hpp"
#include "constants.hpp"
#define NUM_THREADS 1

const Arb hard_lo_b1_range(HARD_LO_B1 - STEP_1_RAD, HARD_LO_B1 + STEP_1_RAD);
const Arb hard_lo_b2_range(HARD_LO_B2 - STEP_1_RAD, HARD_LO_B2 + STEP_1_RAD);
const Arb hard_lo_b12_range(HARD_LO_B12 - STEP_1_RAD, HARD_LO_B12 + STEP_1_RAD);

const Arb hard_hi_b1_range(HARD_HI_B1 - STEP_1_RAD, HARD_HI_B1 + STEP_1_RAD);
const Arb hard_hi_b2_range(HARD_HI_B2 - STEP_1_RAD, HARD_HI_B2 + STEP_1_RAD);
const Arb hard_hi_b12_range(HARD_HI_B12 - STEP_1_RAD, HARD_HI_B12 + STEP_1_RAD);

Arb vol_est(0);
Arb excl_est(0);
Arb tri_est(0);

Arb vol(const Arb &b1_range, const Arb &b2_range, const Arb &b12_range) {
    return b1_range.rad()*b2_range.rad()*b12_range.rad();
}

int check(const Arb &beta_range, const Arb &b1_range, const Arb &b2_range, const Arb &b12_range) {    
    if (hard_lo_b1_range.contains(b1_range) &&
        hard_lo_b2_range.contains(b2_range) &&
        hard_lo_b12_range.contains(b12_range)) {
        excl_est = excl_est + vol(b1_range, b2_range, b12_range);
        return 1;
    }

    if (hard_hi_b1_range.contains(b1_range) &&
        hard_hi_b2_range.contains(b2_range) &&
        hard_hi_b12_range.contains(b12_range)) {
        excl_est = excl_est + vol(b1_range, b2_range, b12_range);
        return 1;
    }

    Max2Sat c(b1_range, b2_range, b12_range);
    
    if (c.triangle_check() == -1) {
        // invalid region, so "good" by default
        tri_est = tri_est + vol(b1_range, b2_range, b12_range);
        return 1;
    }

#ifdef DEBUG
    flint_printf("DEBUG: ");
    beta_range.print();
    flint_printf(": ");
    c.println();
    flint_printf("rho_safe: ");
    c.rho_safe().println();
    flint_printf("VOL : ");
    vol_est.println();
    flint_printf("TRI : ");
    tri_est.println();
    flint_printf("EXCL: ");
    excl_est.println();
#endif

    Arb obj = c.obj(UPPER_CUTOFF_RATIO, beta_range);

#ifdef DEBUG
    flint_printf("obj: ");
    obj.println();
#endif

    if (obj >= 0) {
        vol_est = vol_est + vol(b1_range, b2_range, b12_range);
        return 1;
    }
    if (obj < EPS_LIMIT) {
        flint_printf("FAIL: the check did NOT succeed at ");
        beta_range.println();
        c.ratio(beta_range).println();
        exit(1);
        return 0;
    }

    // otherwise we need to split
    Arb rbeta = beta_range.rad();
    Arb rb1 = b1_range.rad();
    Arb rb2 = b2_range.rad();
    Arb rb12 = b12_range.rad();

    // split beta
    if (rbeta >= rb1 && rbeta >= rb2 && rbeta >= rb12) {
#ifdef DEBUG
        flint_printf("SPLIT: beta\n");
#endif

        return check(beta_range.left_half(), b1_range, b2_range, b12_range) &&
            check(beta_range.right_half(), b1_range, b2_range, b12_range);
    }
    else if (rb1 >= rb2 && rb1 >= rb12) {
#ifdef DEBUG
        flint_printf("SPLIT: b1\n");
#endif

        return check(beta_range, b1_range.left_half(), b2_range, b12_range) &&
            check(beta_range, b1_range.right_half(), b2_range, b12_range);
    }
    else if (rb2 >= rb12) {
#ifdef DEBUG
        flint_printf("SPLIT: b2\n");
#endif
        return check(beta_range, b1_range, b2_range.left_half(), b12_range) &&
            check(beta_range, b1_range, b2_range.right_half(), b12_range);
    }
    else {
#ifdef DEBUG
        flint_printf("SPLIT: b12\n");
#endif

         return check(beta_range, b1_range, b2_range, b12_range.left_half()) &&
            check(beta_range, b1_range, b2_range, b12_range.right_half());
    }
}

int main(int argc, char* argv[]) {
    
    flint_set_num_threads(NUM_THREADS);

    Arb beta_range(BETA_LO, BETA_HI);
    // FIXME to bigger range
    Arb b1_range(-1, 1);
    Arb b2_range(-1, 1);
    Arb b12_range(-1, 1);
    
    flint_printf("Goal ratio: %f\n", UPPER_CUTOFF_RATIO);
    flint_printf("Step 1: rule out configurations outside of\n");
    flint_printf("(");
    hard_lo_b1_range.print();
    flint_printf(" ");
    hard_lo_b2_range.print();
    flint_printf(" ");
    hard_lo_b12_range.print();
    flint_printf(")\n");

    flint_printf("(");
    hard_hi_b1_range.print();
    flint_printf(" ");
    hard_hi_b2_range.print();
    flint_printf(" ");
    hard_hi_b12_range.print();
    flint_printf(")\n");

    int res = check(beta_range, b1_range, b2_range, b12_range);
    
    flint_printf("Result: %d\n", res);
    
    flint_cleanup_master();

    return 0;
}
