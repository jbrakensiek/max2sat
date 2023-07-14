/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include <cstdio>
#include "max2sat.hpp"
#include "constants.hpp"
#define NUM_THREADS 1

int check_lo(Arb beta_range) {
    Max2Sat c(HARD_LO_B1, HARD_LO_B2, HARD_LO_B12);
    Arb obj = c.obj(LOWER_CUTOFF_RATIO, beta_range);

    #ifdef DEBUG
    flint_printf("DEBUG: ");
    beta_range.print();
    flint_printf(", ");
    obj.println();
    #endif
    
    if (obj <= 0) return 1;
    if (obj > (-EPS_LIMIT)) {
        flint_printf("FAIL: the check did NOT succeed at ");
        beta_range.println();
        c.ratio(beta_range).println();
        exit(1);
        return 0;
    }

    return check_lo(beta_range.left_half()) && check_lo(beta_range.right_half());
}

int check_hi(Arb beta_range) {
    Max2Sat c(HARD_HI_B1, HARD_HI_B2, HARD_HI_B12);
    Arb obj = c.obj(LOWER_CUTOFF_RATIO, beta_range);

    #ifdef DEBUG
    flint_printf("DEBUG: ");
    beta_range.print();
    flint_printf(", ");
    obj.println();
    #endif
    
    if (obj <= 0) return 1;
    if (obj > (-EPS_LIMIT)) {
        flint_printf("FAIL: the check did NOT succeed at ");
        beta_range.println();
        c.ratio(beta_range).println();
        exit(1);
        return 0;
    }

    return check_hi(beta_range.left_half()) && check_hi(beta_range.right_half());
}

int main(int argc, char* argv[]) {
    
    flint_set_num_threads(NUM_THREADS);

    // document what we are checking
    Arb beta_range_lo(MINUS_ONE, BETA_LO);
    Arb beta_range_hi(BETA_HI, PLUS_ONE);
    Arb cutoff(LOWER_CUTOFF_RATIO);

    flint_printf("Goal ratio: %f\n", LOWER_CUTOFF_RATIO);

    flint_printf("Step 0A: rule out beta in ");
    beta_range_lo.pretty_println();

    int res1 = check_lo(beta_range_lo);

    flint_printf("Result: %d\n", res1);

    flint_printf("Step 0B: rule out beta in ");
    beta_range_hi.pretty_println();

    int res2 = check_hi(beta_range_hi);

    flint_printf("Result: %d\n", res2);
    flint_printf("Step 0 succeeded.\n");
    
    flint_cleanup_master();

    return 0;
}
