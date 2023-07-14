/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include <cstdio>
#include "bivariate_normal.hpp"
#define NUM_THREADS 1

int main(int argc, char* argv[]) {
    
    flint_set_num_threads(NUM_THREADS);

    flint_printf("%d\n", RHO_THRESH < 1);

    biv_norm_cdf(0,0,0).println();
    biv_norm_cdf_unsafe(0,0,0).println();

    flint_printf("\n");

    biv_norm_cdf(0,0,1).println();
    biv_norm_cdf_unsafe(0,0,0.99999).println();

    flint_printf("\n");

    biv_norm_cdf(0,0,-1).println();
    biv_norm_cdf_unsafe(0,0,-0.99999).println();

    flint_printf("\n");
    
    Arb r = Arb(-0.5,0.5);

    biv_norm_cdf(0,0,r).println();
    biv_norm_cdf_unsafe(0,0,r).println();

    biv_norm_cdf(r,r,r).println();
    biv_norm_cdf_unsafe(r,r,r).println();

    flint_printf("\n");

    Arb r2 = Arb(-1, 0.5);
    Arb r3 = Arb(-0.99999, 0.5);
    biv_norm_cdf(0,0,r2).println();
    biv_norm_cdf(r,r,r2).println();

    biv_norm_cdf(0,0,r3).println();
    biv_norm_cdf(r,r,r3).println();

    biv_norm_cdf_unsafe(0,0,r3).println();
    biv_norm_cdf_unsafe(r,r,r3).println();

    flint_printf("\n");

    Arb r4 = Arb(-1, 1);
    Arb r5 = Arb(-0.9, 0.9);
    biv_norm_cdf(0.5,0.5,r4).println();
    biv_norm_cdf_unsafe(0.5,0.5,r5).println();

    flint_printf("\n");

    Arb a = biv_norm_cdf(0,0,0);
    Arb b = biv_norm_cdf(0.1,0.1,0.1);

    a.println();
    b.println();
    biv_norm_cdf(0,0,0.5).println();
    biv_norm_cdf(Arb(0,0.1),Arb(0,0.1),Arb(0,0.1)).println();
    Arb::join(a,b).println();
    biv_norm_cdf(0,0,Arb(-0.5,0.5)).println();

    biv_norm_cdf(0,0,0.999).println();
    biv_norm_cdf(0,0,0.99).println();
    biv_norm_cdf(0,0,Arb(0.99,0.999)).println();
    biv_norm_cdf(0,0,0.99999).println();
    biv_norm_cdf(0,10,0.1).println();
    biv_norm_cdf(10,10,0.1).println();

    // need more tests for partials

    flint_cleanup_master();

    return 0;
}
