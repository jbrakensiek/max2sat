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

    Max2Sat w(0, 0, 0);
    w.value().println();
    w.prob(0).println();
    w.ratio(0).println();
    w.obj(0, 0).println();
    w.obj(1, 0).println();

    Arb beta(0.94016567);

    /* Data from Uri:
       c =   0.90000  ratio =  0.94285733  b =  0.09489274  0.09489281 -0.78978555  0.02042889
       c =   0.91000  ratio =  0.94348622  b =  0.06370651  0.06370647 -0.76741298  0.10517404
       c =   0.92000  ratio =  0.94382378  b =  0.03615154  0.03615155 -0.75230309  0.17539381
       c =   0.93000  ratio =  0.94394670  b =  0.01008972  0.01008970 -0.74017942  0.23964116
       c =   0.94000  ratio =  0.94388932  b = -0.01518071 -0.01518070 -0.72963860  0.24000000
       c =   0.95000  ratio =  0.94366882  b = -0.04000018 -0.04000018 -0.71999964  0.20000000
       c =   0.96000  ratio =  0.94329275  b = -0.06456802 -0.06456807 -0.71086391  0.16000000
       c =   0.97000  ratio =  0.94276187  b = -0.08901563 -0.08901564 -0.70196873  0.12000000
       c =   0.98000  ratio =  0.94207139  b = -0.11343717 -0.11343722 -0.69312561  0.08000000
       c =   0.99000  ratio =  0.94121122  b = -0.13790514 -0.13790511 -0.68418975  0.04000000
       c =   1.00000  ratio =  0.94016567  b = -0.16247834 -0.16247830 -0.67504335  0.00000000
    */

    Max2Sat c(0.09489274, 0.09489281, -0.78978555);
    c.value().println();
    c.prob(beta).println();
    c.ratio(beta).println();
    c.obj(beta, beta).println();

    Max2Sat d(-0.16247834, -0.16247830, -0.67504335);
    d.value().println();
    d.prob(beta).println();
    d.ratio(beta).println();
    d.obj(beta, beta).println();

    Arb b1(-0.04000018, -0.06456802);
    Arb b12(-0.71999964, -0.71086391);

    Max2Sat e(b1, b1, b12);
    e.value().println();
    e.prob(beta).println();
    e.ratio(beta).println();
    e.obj(beta, beta).println();

    // TODO: tests with b1 != b2

    flint_cleanup_master();

    return 0;
}
