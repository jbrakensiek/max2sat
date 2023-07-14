/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include <cstdio>
#include <cassert>
#include "max2sat.hpp"
#include "constants.hpp"
#define NUM_THREADS 1

Arb vol_est(0);
Arb excl_est(0);
Arb tri_est(0);

Arb vol(const Arb &a, const Arb &b, const Arb &c) {
    return a.rad() * b.rad() * c.rad();
}

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


int check(const Arb &b1, const Arb &b2, const Arb &rho, const Arb &beta) {
    //int t = Config::tri_check_rel_rho(b1, b2, rho);

    Arb b12 = Config::b12_from_rel_rho(b1, b2, rho);

    if (b12 < -1 + Arb::abs(b1 + b2)) {
        // invalid region, so "good" by default
        tri_est = tri_est + vol(b1, b2, rho);
        return 1;
    }

    if (obj(b1, b2, rho, beta) > TYPE_5_LOWER_BOUND) {
        vol_est = vol_est + vol(b1, b2, rho);
        return 1;
    }

    // derivative checks

    if (b12 < 1 - Arb::abs(b1 - b2)) {
        Arb d_b1 = obj_d_b1(b1, b2, rho, beta);
        Arb d_b2 = obj_d_b2(b1, b2, rho, beta);
        Arb d_rho = obj_d_rho(b1, b2, rho, beta);

#ifdef DEBUG
        flint_printf("PARTIALS\n");
        d_b1.pretty_println();
        d_b2.pretty_println();
        d_rho.pretty_println();
#endif

        if (!b1.is_nan() && (d_b1 > 0 || d_b1 < 0)) {
            excl_est = excl_est + vol(b1, b2, rho);
            return 1;
        }
        else if (!b2.is_nan() && (d_b2 > 0 || d_b2 < 0)) {
            excl_est = excl_est + vol(b1, b2, rho);
            return 1;
        }
        else if (!rho.is_nan() && (d_rho > 0 || d_rho < 0)) {
            excl_est = excl_est + vol(b1, b2, rho);
            return 1;
        }
    }

#ifdef DEBUG
    flint_printf("AT :\n");
    b1.println();
    b2.println();
    rho.println();
    beta.println();
    obj(b1,b2,rho,beta).println();
    flint_printf("VOL : ");
    vol_est.println();
    flint_printf("TRI : ");
    tri_est.println();
    flint_printf("EXCL: ");
    excl_est.println();
#endif

    // otherwise we need to split
    Arb rb1 = b1.rad();
    Arb rb2 = b2.rad();
    Arb rrho = rho.rad();
    Arb rbeta = beta.rad();

    if (rb1 >= rb2 && rb1 >= rrho && rb1 >= rbeta) {
        return check(b1.left_half(), b2, rho, beta) && check(b1.right_half(), b2, rho, beta);
    }
    else if (rb2 >= rrho && rb2 >= rbeta) {
        return check(b1, b2.left_half(), rho, beta) &&
                        check(b1, b2.right_half(), rho, beta);
    }
    else if (rrho >= rbeta){
        return check(b1, b2, rho.left_half(), beta) &&
                        check(b1, b2, rho.right_half(), beta);
    }
    else {
        return check(b1, b2, rho, beta.left_half()) &&
            check(b1, b2, rho, beta.right_half());

    }
}

int main(int argc, char* argv[]) {
    
    flint_set_num_threads(NUM_THREADS);
    
    /*    Arb b1(-0.1);
    Arb b2(-0.2);
    Arb rho(0.3);
    Arb eps(0.0000001);

    Arb b12 = Config::b12_from_rel_rho(b1, b2, rho);
    b12.println();
    Max2Sat c(b1, b2, b12);
    flint_printf("%d\n", c.triangle_check());

    Arb ob = obj(b1, b2, rho);
    ob.println();
    (c.prob(beta) - beta* c.value()).println();
    Arb ob_d_b1_est = (obj(b1+eps, b2, rho) - obj(b1, b2, rho)) / eps;
    Arb ob_d_b2_est = (obj(b1, b2+eps, rho) - obj(b1, b2, rho)) / eps;
    Arb ob_d_rho_est = (obj(b1, b2, rho+eps) - obj(b1, b2, rho)) / eps;

    ob_d_b1_est.println();
    obj_d_b1(b1, b2, rho).println();

    ob_d_b2_est.println();
    obj_d_b2(b1, b2, rho).println();

    ob_d_rho_est.println();
    obj_d_rho(b1, b2, rho).println(); */

    Arb b1_range(-1, 1);
    Arb b2_range(-1, 1);
    Arb rho_range(-1, 1);
    Arb beta_range(TYPE_5_COARSE_BETA_LO, TYPE_5_COARSE_BETA_HI);

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

    flint_printf("Step 1: gradient nonzero everywhere (or easy to approx)\n");

    flint_printf("RESULT: %d\n", check(b1_range, b2_range, rho_range, beta_range));

    flint_printf("VOL : ");
    vol_est.println();
    flint_printf("TRI : ");
    tri_est.println();
    flint_printf("EXCL: ");
    excl_est.println();
    
    flint_cleanup_master();

    return 0;
}
