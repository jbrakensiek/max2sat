/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include <cstdio>
#include <cassert>
#include "max2sat.hpp"
#include "constants.hpp"
#define NUM_THREADS 1


Arb obj(const Arb &b1, const Arb &b2, const Arb &rho) {
    //Max2Sat::prob_from_rel(b1, b2, rho, 1).pretty_println();
    //Max2Sat::value_from_rel(b1, b2, rho).pretty_println();
    return Max2Sat::prob_from_rel(b1, b2, rho, 1) - Max2Sat::value_from_rel(b1, b2, rho); 
}

Arb obj_d_b1(const Arb &b1, const Arb &b2, const Arb &rho) {
    return Max2Sat::prob_from_rel_d_b1(b1, b2, rho, 1) - Max2Sat::value_from_rel_d_b1(b1, b2, rho); 
}

Arb obj_d_b2(const Arb &b1, const Arb &b2, const Arb &rho) {
    return Max2Sat::prob_from_rel_d_b2(b1, b2, rho, 1) - Max2Sat::value_from_rel_d_b2(b1, b2, rho); 
}

Arb obj_d_rho(const Arb &b1, const Arb &b2, const Arb &rho) {
    return Max2Sat::prob_from_rel_d_rho(b1, b2, rho, 1) - Max2Sat::value_from_rel_d_rho(b1, b2, rho); 
}


Arb eval_low(const Arb &b1, const Arb &b2) {
    Arb b12 = -1 + Arb::abs(b1 + b2);
    Arb rho = Config::rho_safe(b1, b2, b12);
    return obj(b1, b2, rho);
}

Arb eval_low_pos_d_b1(const Arb &b1, const Arb &b2) {
    Arb b12 = -1 + Arb::abs(b1 + b2);
    Arb rho = Config::rho_safe(b1, b2, b12);

    Arb x = Arb::safe_sqrt((1-b2) / (1 + b2));
    Arb y = 1 / ((1 + b1) * Arb::safe_sqrt(1 - b1.sqr()));

    return obj_d_b1(b1, b2, rho) + obj_d_rho(b1, b2, rho) * x * y;
}

Arb eval_low_pos_d_b2(const Arb &b1, const Arb &b2) {
    Arb b12 = -1 + Arb::abs(b1 + b2);
    Arb rho = Config::rho_safe(b1, b2, b12);

    Arb x = Arb::safe_sqrt((1-b1) / (1 + b1));
    Arb y = 1 / ((1 + b2) * Arb::safe_sqrt(1 - b2.sqr()));

    return obj_d_b2(b1, b2, rho) + obj_d_rho(b1, b2, rho) * x * y;
}

Arb eval_low_neg_d_b1(const Arb &b1, const Arb &b2) {
    Arb b12 = -1 + Arb::abs(b1 + b2);
    Arb rho = Config::rho_safe(b1, b2, b12);

    Arb x = Arb::safe_sqrt((1+b2) / (1 - b2));
    Arb y = 1 / ((1 - b1) * Arb::safe_sqrt(1 - b1.sqr()));

    return obj_d_b1(b1, b2, rho) - obj_d_rho(b1, b2, rho) * x * y;
}

Arb eval_low_neg_d_b2(const Arb &b1, const Arb &b2) {
    Arb b12 = -1 + Arb::abs(b1 + b2);
    Arb rho = Config::rho_safe(b1, b2, b12);

    Arb x = Arb::safe_sqrt((1+b1) / (1 - b1));
    Arb y = 1 / ((1 - b2) * Arb::safe_sqrt(1 - b2.sqr()));

    return obj_d_b2(b1, b2, rho) - obj_d_rho(b1, b2, rho) * x * y;
}


int check(const Arb &b1, const Arb &b2) {
    Arb b12 = -1 + Arb::abs(b1 + b2);
    Arb rho = Config::rho_safe(b1, b2, b12);

    if (obj(b1, b2, rho) >= OBJ_HI) {
        //(1/(1-obj(b1, b2, rho))).println();
        return 1;
    }

    if (b1 + b2 < 0) {
        // Can ignore this case
        return 1;
    }

    if (b1 + b2 > 0) {
        Arb d_b1 = eval_low_pos_d_b1(b1, b2);
        Arb d_b2 = eval_low_pos_d_b2(b1, b2);
        if(d_b1.is_nan()) {
            assert(!(d_b1 > 0));
            assert(!(d_b1 < 0));
        }
        if (d_b2.is_nan()) {
            assert(!(d_b2 > 0));
            assert(!(d_b2 < 0));
        }
        if (d_b1 > 0 || d_b1 < 0) {
            return 1;
        }
        else if (d_b2 > 0 || d_b2 < 0) {
            return 1;
        }
    }

    if (Arb::abs(b1 - TYPE_4_B1_HARD) < TYPE_4_EPS &&
        Arb::abs(b2 - TYPE_4_B2_HARD) < TYPE_4_EPS) {
        // SKIP
        return 1;
    }

    /* flint_printf("STUFF\n");
    b1.println();
    b2.println();
    b12.println();
    (1/(1-eval_low(b1, b2))).println(); */

    // otherwise we need to split
    Arb rb1 = b1.rad();
    Arb rb2 = b2.rad();

    if (rb1 >= rb2) {
#ifdef DEBUG
        flint_printf("SPLIT: b1\n");
#endif
        return check(b1.left_half(), b2) && check(b1.right_half(), b2);
    }
    else {
#ifdef DEBUG
        flint_printf("SPLIT: b2\n");
#endif
        return check(b1, b2.left_half()) &&
                        check(b1, b2.right_half());
    }
}

int main(int argc, char* argv[]) {
    
    flint_set_num_threads(NUM_THREADS);
    
    //(1/(1-eval_low(0.1489, 0.1489))).println();
    //(1/(1-eval_low(-0.1489, -0.1489))).println();
    /* Arb b1(.2);
    Arb b2(.3);
    Arb eps(.0000001);
    Arb x = eval_low(b1, b2);
    Arb y = eval_low(b1 + eps, b2);
    Arb z = eval_low(b1, b2 + eps);
    ((y - x) / eps).println();
    eval_low_pos_d_b1(b1, b2).println();

    ((z - x) / eps).println();
    eval_low_pos_d_b2(b1, b2).println();

    Arb bb1(-.2);
    Arb bb2(-.3);
    Arb xx = eval_low(bb1, bb2);
    Arb yy = eval_low(bb1 + eps, bb2);
    Arb zz = eval_low(bb1, bb2 + eps);
    ((yy - xx) / eps).println();
    eval_low_neg_d_b1(bb1, bb2).println();

    ((zz - xx) / eps).println();
    eval_low_neg_d_b2(bb1, bb2).println(); */


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

    Arb b1_range(-1, 1);
    Arb b2_range(-1, 1);
    
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

    flint_printf("Step 2: easy to approx on the boundary away from (%.10f, %.10f):\n", TYPE_4_B1_HARD, TYPE_4_B1_HARD);

    flint_printf("RESULT: %d\n", check(b1_range, b2_range));

    flint_printf("UPPER BOUND on optimal ratio: "); (1/(1-eval_low(TYPE_4_B1_HARD, TYPE_4_B2_HARD))).println();
    //flint_printf("CHECK: %d\n", eval_low(TYPE_4_B1_HARD, TYPE_4_B2_HARD) < 1 - 1/0.9462);
    
    flint_cleanup_master();

    return 0;
}
