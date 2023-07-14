/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include <cstdio>
#include "arb_wrapper.hpp"
#define NUM_THREADS 1


int main(int argc, char* argv[]) {
    
    flint_set_num_threads(NUM_THREADS);

    Arb x(2.0), y(3.0);

    x.println();
    y.println();

    Arb z = x + y;

    z.println();

    x.join(y).println();
    (Arb::join(x, y)).println();

    (z+2.0).println();

    Arb w(2.0, 3.0);

    w.println();

    (2.0+w).println();
    (w+2).println();

    (w+w).println();
    (w*w).println();

    (7.0-w).println();
    (1 / w).println();

    (w.exp()).println();
    (Arb::exp(w)).println();

    (w.sqrt()).println();
    (Arb::sqrt(w)).println();

    (w.pow(w)).println();
    (Arb::pow(w,w)).println();
    (Arb::pow(2,w)).println();
    (Arb::pow(w,2)).println();

    (Arb::pi()).println();

    (Arb::norm_pdf(0)).println();
    (Arb::norm_pdf(1)).println();

    (Arb::norm_cdf(-1)).println();
    (Arb::norm_cdf(0)).println();
    (Arb::norm_cdf(10)).println();

    (Arb::norm_cdf_inv(0)).println();
    (Arb::norm_cdf_inv(0.025)).println();
    (Arb::norm_cdf_inv(0.5)).println();
    (Arb::norm_cdf_inv(0.975)).println();
    (Arb::norm_cdf_inv(1)).println();
    (Arb::norm_cdf_inv(Arb::join(0.98,0.99))).println();

    (w.mid()).println();
    w.intersect(w+0.5).println();
    w.intersect(w+1.5).println(); // should be NaN
    w.join(w+0.5).println();

    w.left_half().println();
    w.right_half().println();

    flint_printf("%d\n", w.is_nan());
    flint_printf("%d\n", w.intersect(w+1.5).is_nan());
    flint_printf("%d\n", w == w);
    flint_printf("%d\n", w != w);
    flint_printf("%d\n", w != (w+1));
    flint_printf("%d\n", w != (w+1.01));
    flint_printf("%d\n", w < (w+1));
    flint_printf("%d\n", w < (w+1.000001));
    flint_printf("%d\n", w <= (w+1));
    flint_printf("%d\n", w <= (w+1.000001));
    flint_printf("%d\n", (w+1.000001) > w);
    flint_printf("%d\n", (w+1.000001) >= w);
    flint_printf("%d\n", ((Arb) 1) == ((Arb) 1));
    flint_printf("%d\n", ((Arb) 1) > ((Arb) 1));
    flint_printf("%d\n", ((Arb) 1) >= ((Arb) 1));

    flint_printf("%d\n", sizeof(Arb));
    flint_printf("%d\n", sizeof(arb_t));

    w.println();
    w.left_edge().println();
    w.right_edge().println();

    flint_printf("%d\n", w >= w.left_edge());
    flint_printf("%d\n", w <= w.right_edge());

    w.left_half().println();
    w.right_half().println();
    (2*w-2.5).println();
    w.min(2*w-2.5).println();
    w.max(2*w-2.5).println();

    Arb::abs(w-2.5).println();
    Arb::abs(w-3.5).println();

    flint_cleanup_master();

    return 0;
}
