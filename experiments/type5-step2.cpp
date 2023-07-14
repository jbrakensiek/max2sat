/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include <cstdio>
#include <cassert>
#include "max2sat.hpp"
#include "constants.hpp"
#define NUM_THREADS 1


Arb obj(const Arb &b1, const Arb &b2, const Arb &rho, const Arb &beta) {
//    Max2Sat::prob_from_rel(b1, b2, rho, 1).pretty_println();
//    Max2Sat::value_from_rel(b1, b2, rho).pretty_println();
    return Max2Sat::prob_from_rel(b1, b2, rho, beta)
        - beta * Max2Sat::value_from_rel(b1, b2, rho); 
}

Arb obj_d_b1(const Arb &b1, const Arb &b2, const Arb &rho, const Arb &beta) {
    return Max2Sat::prob_from_rel_d_b1(b1, b2, rho, beta)
        - beta * Max2Sat::value_from_rel_d_b1(b1, b2, rho); 
}

Arb obj_d_b2(const Arb &b1, const Arb &b2, const Arb &rho, const Arb &beta) {
    return Max2Sat::prob_from_rel_d_b2(b1, b2, rho, beta)
        - beta * Max2Sat::value_from_rel_d_b2(b1, b2, rho); 
}

Arb obj_d_rho(const Arb &b1, const Arb &b2, const Arb &rho, const Arb &beta) {
    return Max2Sat::prob_from_rel_d_rho(b1, b2, rho, beta)
        - beta *  Max2Sat::value_from_rel_d_rho(b1, b2, rho); 
}

Arb eval_low(const Arb &b1, const Arb &b2, const Arb &beta) {
    Arb b12 = -1 + Arb::abs(b1 + b2);
    Arb rho = Config::rho_safe(b1, b2, b12);
    return obj(b1, b2, rho, beta);
}

Arb eval_low_pos_d_b1(const Arb &b1, const Arb &b2, const Arb &beta) {
    Arb b12 = -1 + Arb::abs(b1 + b2);
    Arb rho = Config::rho_safe(b1, b2, b12);

    Arb x = Arb::safe_sqrt((1-b2) / (1 + b2));
    Arb y = 1 / ((1 + b1) * Arb::safe_sqrt(1 - b1.sqr()));

    return obj_d_b1(b1, b2, rho, beta) + obj_d_rho(b1, b2, rho, beta) * x * y;
}

Arb eval_low_pos_d_b2(const Arb &b1, const Arb &b2, const Arb &beta) {
    Arb b12 = -1 + Arb::abs(b1 + b2);
    Arb rho = Config::rho_safe(b1, b2, b12);

    Arb x = Arb::safe_sqrt((1-b1) / (1 + b1));
    Arb y = 1 / ((1 + b2) * Arb::safe_sqrt(1 - b2.sqr()));

    return obj_d_b2(b1, b2, rho, beta) + obj_d_rho(b1, b2, rho, beta) * x * y;
}

int pos_check(const Arb &b1, const Arb &b2, const Arb &beta) {
    Arb b12 = -1 + Arb::abs(b1 + b2);
    Arb rho = Config::rho_safe(b1, b2, b12);
    Arb ob = obj(b1, b2, rho, beta);

    if (ob > 0) {
        return 1;
    }

    /*flint_printf("POS_CHECK\n");
    b1.println();
    b2.println();
    ob.println();*/
    
    Arb rb1 = b1.rad();
    Arb rb2 = b2.rad();

    if (rb1 >= rb2) {
        //flint_printf("SPLIT: b1\n");
        return pos_check(b1.left_half(), b2, beta) &&
            pos_check(b1.right_half(), b2, beta);
    }
    else {
        //flint_printf("SPLIT: b2\n");
        return pos_check(b1, b2.left_half(), beta) &&
            pos_check(b1, b2.right_half(), beta);
    }
}

int check(const Arb &b1, const Arb &b2, const Arb &beta) {
    Arb b12;
    b12 = -1 + Arb::abs(b1 + b2);

    Arb rho = Config::rho_safe(b1, b2, b12);

    if (obj(b1, b2, rho, beta) > TYPE_5_LOWER_BOUND) {
        //(1/(1-obj(b1, b2, rho))).println();
        return 1;
    }

    if (b1 + b2 < 0) {
        // Can ignore this case by symmetry
        return 1;
    }

    if (b1 + b2 > 0) {
        Arb d_b1 = eval_low_pos_d_b1(b1, b2, beta);
        Arb d_b2 = eval_low_pos_d_b2(b1, b2, beta);
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

    if (Arb::abs(b1 - TYPE_5_B1_HARD) < TYPE_5_EPS &&
        Arb::abs(b2 - TYPE_5_B2_HARD) < TYPE_5_EPS) {
        // SKIP as too close
        return 1;
    }

/*    flint_printf("STUFF\n");
    b1.println();
    b2.println();
    b12.println();
    beta.println();
    eval_low(b1, b2, beta).println(); */

    //assert(!(obj(b1, b2, rho, beta) < 0));

    // otherwise we need to split
    Arb rb1 = b1.rad();
    Arb rb2 = b2.rad();
    Arb rbeta = beta.rad();

    if (rb1 >= rb2 && rb1 >= rbeta) {
#ifdef DEBUG
        flint_printf("SPLIT: b1\n");
#endif
        return check(b1.left_half(), b2, beta) &&
            check(b1.right_half(), b2, beta);
    }
    else if (rb2 >= rbeta) {
#ifdef DEBUG
        flint_printf("SPLIT: b2\n");
#endif
        return check(b1, b2.left_half(), beta) &&
            check(b1, b2.right_half(), beta);
    }
    else {
        return check(b1, b2, beta.left_half()) &&
            check(b1, b2, beta.right_half());
    }
}

int main(int argc, char* argv[]) {
    
    flint_set_num_threads(NUM_THREADS);
    
    //(1/(1-eval_low(0.1489, 0.1489))).println();
    //(1/(1-eval_low(-0.1489, -0.1489))).println();
    Arb b1(.2);
    Arb b2(.3);
    Arb eps(.0000001);
    Arb beta(0.901655);
    Arb x = eval_low(b1, b2, beta);
    Arb y = eval_low(b1 + eps, b2, beta);
    Arb z = eval_low(b1, b2 + eps, beta);
    ((y - x) / eps).println();
    eval_low_pos_d_b1(b1, b2, beta).println();

    ((z - x) / eps).println();
    eval_low_pos_d_b2(b1, b2, beta).println();

    /* Arb bb1(-.2);
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
    Arb beta_range(TYPE_5_BETA_LO, TYPE_5_BETA_HI);

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

    flint_printf("Step 2: easy to approx on the boundary away from (%f, %f):\n", TYPE_5_B1_HARD, TYPE_5_B2_HARD);

    flint_printf("RESULT: %d\n", check(b1_range, b2_range, beta_range));

    flint_printf("NEG SIGN AT %.10f: %d\n", TYPE_5_BETA_HI,
                 eval_low(TYPE_5_B1_HARD, TYPE_5_B2_HARD,TYPE_5_BETA_HI) < 0);
    flint_printf("POS SIGN AT %.10f: %d\n", TYPE_5_BETA_LO,
                 pos_check(Arb(TYPE_5_B1_HARD-2*TYPE_5_EPS, TYPE_5_B1_HARD+2*TYPE_5_EPS),
                           Arb(TYPE_5_B2_HARD-2*TYPE_5_EPS, TYPE_5_B2_HARD+2*TYPE_5_EPS), TYPE_5_BETA_LO));

    
    flint_cleanup_master();

    return 0;
}
