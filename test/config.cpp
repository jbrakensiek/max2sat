/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include <cstdio>
#include "config.hpp"
#include "max2sat.hpp"
#define NUM_THREADS 1

int main(int argc, char* argv[]) {
    
    flint_set_num_threads(NUM_THREADS);

    Config w(0,0,0);
    Config w1 = Config::from_relative(0,0,0.1);
    Config w2 = Config::from_relative(0,0,0.9);

    w.println();
    w1.println();
    w2.println();

    w.lower_b12().println();
    w.upper_b12().println();

    flint_printf("%d\n", w.triangle_check());

    Config w3(0.5, 0.4, -0.3); 

    w3.println();
    flint_printf("%d\n", w3.triangle_check());

    Config w4(Arb(0,0.5), Arb(0,0.4), Arb(0,-0.3));
    w4.println();
    flint_printf("%d\n", w4.triangle_check());

    w.rho().println();
    w1.rho().println();
    w2.rho().println();
    w3.rho().println();
    w4.rho().println();

    Max2Sat w5(0, 0, 0);
    w5.value().println();

    flint_cleanup_master();

    return 0;
}
