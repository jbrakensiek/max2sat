/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include <cstdio>
#include "max2sat.hpp"
#include "constants.hpp"
#define NUM_THREADS 1

const Arb b12_range(HARD_LO_B12 - STEP_3_RAD, HARD_LO_B12 + STEP_3_RAD);

Arb g(const Arb &t, const Arb &beta, const Arb &b12) {
    Arb c12 = (1+b12)/2;
    Max2Sat c1(t-c12,t-c12,b12);
    Max2Sat c2(-c12,-c12,b12);
    return c1.prob(beta) - (1 - t/2) * c2.prob(beta);
}

Arb g_prime(const Arb &t, const Arb &beta, const Arb &b12) {
    Arb c12 = (1 + b12) / 2;
    Arb rho = (b12 - 1) / (1 - (t - c12).sqr()) + 1;
    Arb drho_dt = 2 * (b12 - 1) * (t - c12) / (1 - (t - c12).sqr()).sqr();
    Arb T = Arb::norm_cdf_inv((1 + beta * (t - c12)) / 2);
    Arb dPhi_drho = biv_norm_cdf_d_rho(T, T, rho);
    Arb dPhi_dx = biv_norm_cdf_d_t1(T, T, rho);
    Arb dT_dt = beta / (2 * Arb::norm_pdf(T));

    Arb z = 2 * dT_dt * dPhi_dx;
    //Arb z = beta * Arb::norm_cdf(T * ((1-rho)/(1+rho)).sqrt());

    Max2Sat conf(-c12, -c12, b12);
    Arb p = conf.prob(beta);

    return - (drho_dt * dPhi_drho + z) + (p/2);
}

int g_check(const Arb &t, const Arb &beta, const Arb &b12) {
    Arb res = g_prime(t, beta, b12);

#ifdef DEBUG
    t.print();
    flint_printf(" ");
    res.print();
    flint_printf(" %d\n", res > 0);
#endif

    if (res > 0) return 1;
    if (res < 0) return 0;

    return g_check(t.left_half(), beta, b12.left_half()) &&
        g_check(t.right_half(), beta, b12.left_half()) &&
        g_check(t.left_half(), beta, b12.right_half()) &&
        g_check(t.right_half(), beta, b12.right_half());
}

Arb h(const Arb &t, const Arb &beta, const Arb &b12) {
    Arb c12 = (1+b12)/2;
    Max2Sat c1(c12-t,c12-t,b12);
    Max2Sat c2(c12,c12,b12);
    return c1.prob(beta)*((1-b12)/2) - ((1-b12+t)/2) * c2.prob(beta);
}

Arb h_prime(const Arb &t, const Arb &beta, const Arb &b12) {
    Arb c12 = (1 + b12) / 2;
    Arb rho = (b12 - 1) / (1 - (t - c12).sqr()) + 1;
    Arb drho_dt = 2 * (b12 - 1) * (t - c12) / (1 - (t - c12).sqr()).sqr();
    Arb T = Arb::norm_cdf_inv((1 + beta * (c12 - t)) / 2);
    Arb dPhi_drho = biv_norm_cdf_d_rho(T, T, rho);
    Arb dPhi_dx = biv_norm_cdf_d_t1(T, T, rho);
    Arb dT_dt = -beta / (2 * Arb::norm_pdf(T));

    //Arb z = 2 * dT_dt * dPhi_dx;
    Arb z = -beta * Arb::norm_cdf(T * ((1-rho)/(1+rho)).sqrt());

    Max2Sat conf(c12, c12, b12);
    Arb p = conf.prob(beta);

    return - ((1-b12)/2) * (drho_dt * dPhi_drho + z) - (p/2);
}

int h_check(const Arb &t, const Arb &beta, const Arb &b12) {
    Arb res = h_prime(t, beta, b12);

#ifdef DEBUG
    t.print();
    flint_printf(" ");
    res.print();
    flint_printf(" %d\n", res > 0);
#endif

    if (res > 0) return 1;
    if (res < 0) return 0;

    return h_check(t.left_half(), beta, b12.left_half()) &&
        h_check(t.right_half(), beta, b12.left_half()) &&
        h_check(t.left_half(), beta, b12.right_half()) &&
        h_check(t.right_half(), beta, b12.right_half());
}

int main(int argc, char* argv[]) {
    
    flint_set_num_threads(NUM_THREADS);

    Arb beta_range(BETA_LO, BETA_HI);
    Arb t_range(0, T_MAX);

#ifdef DEBUG
   flint_printf("DEBUG G\n");
   g_prime(0, 0.94016567, -0.67504335).println();

   Arb t = 0.00000001;
   Arb g0 = g(0, 0.94016567, -0.67504335);
   g0.println();
   g(0.1, 0.94016567, -0.67504335).println();
   Arb gt = g(t, 0.94016567, -0.67504335);
   ((gt-g0)/t).println();
   flint_printf("\n");

   flint_printf("DEBUG H\n");
   h_prime(0, 0.94016567, -0.67504335).println();

   Arb h0 = h(0, 0.94016567, -0.67504335);
   h0.println();
   h(0.1, 0.94016567, -0.67504335).println();
   Arb ht = h(t, 0.94016567, -0.67504335);
   ((ht-h0)/t).println();
   flint_printf("\n");

#endif
    
 
    flint_printf("Step 3a: prove g_prime(t) >= 0, where\n");
    flint_printf("t   : ");
    t_range.println();
    flint_printf("beta: ");
    beta_range.println();
    flint_printf("b12 : ");
    b12_range.println();

//    Arb res = g_prime(t_range, beta_range, hard_lo_b12_range);
    
//    flint_printf("Result: ");
//    res.println();
    
    flint_printf("Positive? %d\n", g_check(t_range, beta_range, b12_range));

    flint_printf("\nStep 3b: prove h_prime(t) >= 0, where\n");
    flint_printf("t   : ");
    t_range.println();
    flint_printf("beta: ");
    beta_range.println();
    flint_printf("b12 : ");
    b12_range.println();

//    Arb res = g_prime(t_range, beta_range, hard_lo_b12_range);
    
//    flint_printf("Result: ");
//    res.println();
    
    flint_printf("Positive? %d\n", h_check(t_range, beta_range, b12_range));


    //int res = check(beta_range, b1_range, b2_range, b12_range);
    
    //flint_printf("Result: %d\n", res);
    
    flint_cleanup_master();

    return 0;
}
