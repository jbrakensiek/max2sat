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

Arb obj(const Arb &b1, const Arb &b2, const Arb &rho) {
 //   Max2Sat::prob_from_rel(b1, b2, rho, 1).pretty_println();
 //   Max2Sat::value_from_rel(b1, b2, rho).pretty_println();
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


Arb refine_alpha(const Arb &b1, const Arb &b2, const Arb &rho, const Arb &alpha, int depth) {
    int t = Config::tri_check_rel_rho(b1, b2, rho);

    if (t == -1) {
        // invalid region, so "good" by default
        tri_est = tri_est + vol(b1, b2, rho);
        return 1;
    }

    Arb ob = obj(b1, b2, rho);
    assert(! ob.is_nan());
#ifdef DEBUG
    flint_printf("obj: ");
    ob.println();
#endif
    
    // no progress to be made in this part of the tree
    if (depth == 0 || ob >= alpha) { 
        /*if (ob.left_edge() < -0.5) {
            flint_printf("WHAT %d %d\n", depth, t);
            b1.println();
            b2.println();
            rho.println();
            ob.pretty_println();
        }*/
        vol_est = vol_est + vol(b1, b2, rho);
        return ob;
    }

    if (t == 1) {
        // derivative checks -- only valid inside an entirely good region

        Arb d_b1 = obj_d_b1(b1, b2, rho);
        Arb d_b2 = obj_d_b2(b1, b2, rho);
        Arb d_rho = obj_d_rho(b1, b2, rho);

#ifdef DEBUG
        flint_printf("PARTIALS\n");
        d_b1.pretty_println();
        d_b2.pretty_println();
        d_rho.pretty_println();
#endif

        if (b1.rad() > 0 && d_b1 > 0) {
            excl_est = excl_est + vol(b1, b2, rho);
#ifdef DEBUG
        flint_printf("PATH B1-\n");
#endif
            return refine_alpha(b1.left_edge(), b2, rho, alpha, depth-1);
        }
        else if (b1.rad() > 0 && d_b1 < 0) {
#ifdef DEBUG
        flint_printf("PATH B1+\n");
#endif
            excl_est = excl_est + vol(b1, b2, rho);
            return refine_alpha(b1.right_edge(), b2, rho, alpha, depth-1);
        }
        else if (b2.rad() > 0 && d_b2 > 0) {
#ifdef DEBUG
        flint_printf("PATH B2-\n");
#endif
            excl_est = excl_est + vol(b1, b2, rho);
            return refine_alpha(b1, b2.left_edge(), rho, alpha, depth-1);
        }
        else if (b2.rad() > 0 && d_b2 < 0) {
#ifdef DEBUG
        flint_printf("PATH B2+\n");
#endif
            excl_est = excl_est + vol(b1, b2, rho);
            return refine_alpha(b1, b2.right_edge(), rho, alpha, depth-1);
        }
        else if (rho.rad() > 0 && d_rho > 0) {
#ifdef DEBUG
        flint_printf("PATH rho-\n");
#endif
            excl_est = excl_est + vol(b1, b2, rho);
            return refine_alpha(b1, b2, rho.left_edge(), alpha, depth-1);
        }
        else if (rho.rad() > 0 && d_rho < 0) {
#ifdef DEBUG
        flint_printf("PATH rho+\n");
#endif
            excl_est = excl_est + vol(b1, b2, rho);
            return refine_alpha(b1, b2, rho.right_edge(), alpha, depth-1);
        }
#ifdef DEBUG
        flint_printf("STAY\n");
#endif

    }

#ifdef DEBUG
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

    if (rb1 >= rb2 && rb1 >= rrho) {
#ifdef DEBUG
        flint_printf("SPLIT: b1\n");
#endif
        return Arb::min(refine_alpha(b1.left_half(), b2, rho, alpha, depth-1),
                        refine_alpha(b1.right_half(), b2, rho, alpha, depth-1));
    }
    else if (rb2 >= rrho) {
#ifdef DEBUG
        flint_printf("SPLIT: b2\n");
#endif
        return Arb::min(refine_alpha(b1, b2.left_half(), rho, alpha, depth-1),
                        refine_alpha(b1, b2.right_half(), rho, alpha, depth-1));
    }
    else {
#ifdef DEBUG
        flint_printf("SPLIT: rho\n");
#endif
        return Arb::min(refine_alpha(b1, b2, rho.left_half(), alpha, depth-1),
                        refine_alpha(b1, b2, rho.right_half(), alpha, depth-1));
    }
}

int main(int argc, char* argv[]) {
    
    flint_set_num_threads(NUM_THREADS);
    
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
    Arb rho_range(-1, 1);
    
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

    Arb alpha_range(-5, 5);

    for (int i = 0; i < 40; i++) {
        vol_est = 0;
        tri_est = 0;
        excl_est = 0;
        alpha_range = alpha_range.intersect(refine_alpha(b1_range, b2_range, rho_range, alpha_range, i));
        flint_printf("DEPTH %d: ", i);
        alpha_range.pretty_println();
        (1/(1-alpha_range)).pretty_println();
        flint_printf("VOL : ");
        vol_est.println();
        flint_printf("TRI : ");
        tri_est.println();
        flint_printf("EXCL: ");
        excl_est.println();
    }
    
    //flint_printf("Result: %d\n", res);
    
    flint_cleanup_master();

    return 0;
}
