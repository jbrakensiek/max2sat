/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include <cstdio>
#include <cassert>
#include "max2sat.hpp"
#include "constants.hpp"
#define NUM_THREADS 1



Arb eval_norm(const Arb &t1, const Arb &t2) {
    Arb b(TYPE_4_B1_HARD - TYPE_4_HARD_EPS, TYPE_4_B1_HARD + TYPE_4_HARD_EPS);
    Arb b12 = -1 + 2*b;
    Arb rho = Config::rho_safe(b, b, b12);
    Arb ans_1 = TYPE_4_P1 * (1 - biv_norm_cdf_norm_thresh(t1, t1, rho));
    Arb ans_2 = TYPE_4_P2 * (1 - biv_norm_cdf_norm_thresh(t2, t2, rho));
    Arb ans_3 = TYPE_4_P3 * (1 - biv_norm_cdf_norm_thresh(1-t1, t2, rho));
    Arb ans_4 = TYPE_4_P4 * (1 - biv_norm_cdf_norm_thresh(1-t2, t1, rho));
    Arb ans_5 = TYPE_4_P5 * t2;
    Arb ans_6 = TYPE_4_P6 * t1;
    Arb num = ans_1 + ans_2 + ans_3 + ans_4 + ans_5 + ans_6;
    Arb denom = TYPE_4_P1 + TYPE_4_P2 + TYPE_4_P3 + TYPE_4_P4 + TYPE_4_P5 + TYPE_4_P6;
    assert (! (num < -1e-9));
    return num / denom;
}


int check(const Arb &t1, const Arb &t2) {
    if (Arb::abs(t1 - TYPE_4_T1) < TYPE_4_T_EPS &&
        Arb::abs(t2 - TYPE_4_T2) < TYPE_4_T_EPS) {
        // too close, so stop
        return 1;
    }

    if (t1 < TYPE_4_T_EPS && t2 < TYPE_4_T_EPS) {
        return 1;
    }
    
    if (1 - t1 < TYPE_4_T_EPS && 1 - t2 < TYPE_4_T_EPS) {
        return 1;
    }

    if (eval_norm(t1, t2) < TYPE_4_HARD_BOUND) {
        return 1;
    }

    // otherwise we need to split
    Arb rt1 = t1.rad();
    Arb rt2 = t2.rad();

    flint_printf("AT: \n");
    t1.println();
    t2.println();
    eval_norm(t1, t2).println();

    if (rt1 >= rt2) {
#ifdef DEBUG
        flint_printf("SPLIT: t1\n");
#endif
        return check(t1.left_half(), t2) && check(t1.right_half(), t2);
    }
    else {
#ifdef DEBUG
        flint_printf("SPLIT: t2\n");
#endif
        return check(t1, t2.left_half()) &&
                        check(t1, t2.right_half());
    }
}

int main(int argc, char* argv[]) {
    
    flint_set_num_threads(NUM_THREADS);
    
    //(1/(1-eval_low(0.1489, 0.1489))).println();
    //(1/(1-eval_low(-0.1489, -0.1489))).println();
    Arb t1(.2);
    Arb t2(.3);
    eval_norm(0, 0).println();
    eval_norm(TYPE_4_T1, TYPE_4_T2).println();
    eval_norm(1, 1).println();

    /*
    Arb b1(-0.1);
    Arb b2(-0.2);
    Arb rho(0.3);
    Arb eps(0.0000001);

    Arb b12 = Config::b12_from_rel_rho(b1, b2, rho);
    b12.println();
    Max2Sat c(b1, b2, b12);
    flint_printf("%d\n", c.triangle_check());

    Arb ob = obj(b1, b2, rho);
    ob.println();
    (c.prob(1) - c.value()).println();
    Arb ob_d_b1_est = (obj(b1+eps, b2, rho) - obj(b1, b2, rho)) / eps;
    Arb ob_d_b2_est = (obj(b1, b2+eps, rho) - obj(b1, b2, rho)) / eps;
    Arb ob_d_rho_est = (obj(b1, b2, rho+eps) - obj(b1, b2, rho)) / eps;

    ob_d_b1_est.println();
    obj_d_b1(b1, b2, rho).println();

    ob_d_b2_est.println();
    obj_d_b2(b1, b2, rho).println();

    ob_d_rho_est.println();
    obj_d_rho(b1, b2, rho).println();
    */

    Arb t1_range(0, 1);
    Arb t2_range(0, 1);
    
     /*   flint_printf("Goal ratio: %f\n", UPPER_CUTOFF_RATIO);
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
    flint_printf(")\n");*/

    flint_printf("Type 4, hardness\n");

    flint_printf("RESULT: %d\n", check(t2_range, t2_range));
    
    flint_cleanup_master();

    return 0;
}
