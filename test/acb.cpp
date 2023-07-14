/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include <cstdio>
#include "acb_wrapper.hpp"
#define NUM_THREADS 1


int main(int argc, char* argv[]) {
    
    flint_set_num_threads(NUM_THREADS);

    Acb x(1.0), y(0, 1.0), z(2.0, 3.0);

    (x+y).println();
    x.real().println();
    z.imag().println();

    Acb w(x.real(),y.imag());
    Acb v = Acb::join(w,z);

    v.println();
    (v+v).println();
    (v-v).println();
    (v*v).println();
    (v/v).println();
    (x.real()*v).println();

    Acb::pi().println();
    
    flint_printf("%d\n", Acb::pi().is_real());
    flint_printf("%d\n", v.is_real());

    x.sqrt().println();
    v.sqrt().println();

    Acb iv(Arb(-1,1));

    Acb::sqrt_analytic(iv, 1).println();
    Acb::sqrt_analytic(iv, 0).println();

    Acb::exp(iv).println();

    Acb::pow(v, v).println();
    Acb::pow_analytic(iv, v, 0).println();
    Acb::pow_analytic(iv, v, 1).println();

    Acb::norm_pdf(v).println();
    Acb::norm_cdf(v).println();

    Acb(0.0).norm_pdf().println();
    Acb(1.96).norm_cdf().println();

    flint_printf("%d\n", v.is_nan());
    flint_printf("%d\n", Acb::nan().is_nan());
    flint_printf("%d\n", sizeof(Acb));
    flint_printf("%d\n", sizeof(acb_t));

    flint_cleanup_master();

    return 0;
}
