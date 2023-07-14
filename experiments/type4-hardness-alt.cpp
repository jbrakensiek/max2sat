/*
  Copyright (c) 2022-23 Joshua Brakensiek, Neng Huang, Aaron Potechin and Uri Zwick

  This code is licensed under the MIT License.
*/

#include <cstdio>
#include <cassert>
#include "max2sat.hpp"
#include "constants.hpp"
#define NUM_THREADS 1

Arb prob (const Arb &t1, const Arb &t2, const Arb &b) {
    Arb rho = - (1 - b) / (1 + b);
    Arb t = (1-b) / 2;
    Arb r12 = 1 - Arb::norm_cdf(Arb::safe_sqrt((1 - rho) / (1 + rho)) * Arb::norm_cdf_inv(t));
    Arb s1 = r12 * (1 - biv_norm_cdf_norm_thresh(t, t, rho)) +
        (1 - r12) * (1 - biv_norm_cdf_norm_thresh(1-t, 1-t, rho));
    Arb s2 = 1;
    Arb diff = s2 - s1;
    Arb p3 = (diff * 2) / (1 + diff * 2);
    Arb p1 = r12 * (1 - p3);
    Arb p2 = 1 - p1 - p3;
    /* flint_printf("p1: "); p1.println();
    flint_printf("p2: "); p2.println();
    flint_printf("p3: "); p3.println(); */

    Arb q1 = p3 / 2;
    Arb q2 = q1;
    Arb q3 = p2 - p3/2;
    Arb q4 = p1 - p3/2;
    Arb q5 = q1;
    Arb q6 = q1;

    /*flint_printf("q1: "); q1.println();
    flint_printf("q2: "); q2.println();
    flint_printf("q3: "); q3.println();
    flint_printf("q4: "); q4.println();
    flint_printf("q5: "); q5.println();
    flint_printf("q6: "); q6.println();*/

    Arb ans_1 = q1 * (1 - biv_norm_cdf_norm_thresh(t1, t1, rho));
    Arb ans_2 = q2 * (1 - biv_norm_cdf_norm_thresh(t2, t2, rho));
    Arb ans_3 = q3 * (1 - biv_norm_cdf_norm_thresh(1-t1, t2, rho));
    Arb ans_4 = q4 * (1 - biv_norm_cdf_norm_thresh(1-t2, t1, rho));
    Arb ans_5 = q5 * t2;
    Arb ans_6 = q6 * t1;
    
    Arb num = ans_1 + ans_2 + ans_3 + ans_4 + ans_5 + ans_6;
    assert(! (num < -1e-9));

    return num;
}

Arb prob_d_t1(const Arb &t1, const Arb &t2, const Arb &b) {
    Arb rho = - (1 - b) / (1 + b);
    Arb t = (1-b) / 2;
    Arb r12 = 1 - Arb::norm_cdf(Arb::safe_sqrt((1 - rho) / (1 + rho)) * Arb::norm_cdf_inv(t));
    Arb s1 = r12 * (1 - biv_norm_cdf_norm_thresh(t, t, rho)) +
        (1 - r12) * (1 - biv_norm_cdf_norm_thresh(1-t, 1-t, rho));
    Arb s2 = 1;
    Arb diff = s2 - s1;
    Arb p3 = (diff * 2) / (1 + diff * 2);
    Arb p1 = r12 * (1 - p3);
    Arb p2 = 1 - p1 - p3;

    Arb q1 = p3 / 2;
    Arb q2 = q1;
    Arb q3 = p2 - p3/2;
    Arb q4 = p1 - p3/2;
    Arb q5 = q1;
    Arb q6 = q1;

    Arb t1i = Arb::norm_cdf_inv(t1);
    Arb t2i = Arb::norm_cdf_inv(t2);

    Arb pdf_t1i = Arb::norm_pdf(t1i);
    Arb y = pdf_t1i;

    Arb coef_1 = -2 * Arb::norm_cdf(Arb::safe_sqrt((1 - rho) / (1 + rho)) * t1i);
    Arb coef_2 = 0;
    Arb coef_3 = Arb::norm_cdf((t2i + rho * t1i) / Arb::safe_sqrt(1 - rho.sqr()));
    Arb coef_4 = -Arb::norm_cdf((-t2i - rho * t1i) / Arb::safe_sqrt(1 - rho.sqr()));
    Arb coef_5 = 0;
    Arb coef_6 = 1;

    return y * (coef_1 * q1 + coef_2 * q2 + coef_3 * q3 + coef_4 * q4 + coef_5 * q5 + coef_6 * q6);
}

Arb prob_d_t2(const Arb &t1, const Arb &t2, const Arb &b) {
        Arb rho = - (1 - b) / (1 + b);
    Arb t = (1-b) / 2;
    Arb r12 = 1 - Arb::norm_cdf(Arb::safe_sqrt((1 - rho) / (1 + rho)) * Arb::norm_cdf_inv(t));
    Arb s1 = r12 * (1 - biv_norm_cdf_norm_thresh(t, t, rho)) +
        (1 - r12) * (1 - biv_norm_cdf_norm_thresh(1-t, 1-t, rho));
    Arb s2 = 1;
    Arb diff = s2 - s1;
    Arb p3 = (diff * 2) / (1 + diff * 2);
    Arb p1 = r12 * (1 - p3);
    Arb p2 = 1 - p1 - p3;

    Arb q1 = p3 / 2;
    Arb q2 = q1;
    Arb q3 = p2 - p3/2;
    Arb q4 = p1 - p3/2;
    Arb q5 = q1;
    Arb q6 = q1;

    Arb t1i = Arb::norm_cdf_inv(t1);
    Arb t2i = Arb::norm_cdf_inv(t2);

    Arb pdf_t2i = Arb::norm_pdf(t2i);
    Arb y = pdf_t2i;

    Arb coef_1 = 0;
    Arb coef_2 = -2 * Arb::norm_cdf(Arb::safe_sqrt((1 - rho) / (1 + rho)) * t2i);
    Arb coef_3 = -Arb::norm_cdf((-t1i - rho * t2i) / Arb::safe_sqrt(1 - rho.sqr()));
    Arb coef_4 = Arb::norm_cdf((t1i + rho * t2i) / Arb::safe_sqrt(1 - rho.sqr()));
    Arb coef_5 = 1;
    Arb coef_6 = 0;

    return y * (coef_1 * q1 + coef_2 * q2 + coef_3 * q3 + coef_4 * q4 + coef_5 * q5 + coef_6 * q6);
}

Arb prob_d_t1_d_t1(const Arb &t1, const Arb &t2, const Arb &b) {
    Arb rho = - (1 - b) / (1 + b);
    Arb t = (1-b) / 2;
    Arb r12 = 1 - Arb::norm_cdf(Arb::safe_sqrt((1 - rho) / (1 + rho)) * Arb::norm_cdf_inv(t));
    Arb s1 = r12 * (1 - biv_norm_cdf_norm_thresh(t, t, rho)) +
        (1 - r12) * (1 - biv_norm_cdf_norm_thresh(1-t, 1-t, rho));
    Arb s2 = 1;
    Arb diff = s2 - s1;
    Arb p3 = (diff * 2) / (1 + diff * 2);
    Arb p1 = r12 * (1 - p3);
    Arb p2 = 1 - p1 - p3;

    Arb q1 = p3 / 2;
    Arb q2 = q1;
    Arb q3 = p2 - p3/2;
    Arb q4 = p1 - p3/2;
    Arb q5 = q1;
    Arb q6 = q1;

    Arb t1i = Arb::norm_cdf_inv(t1);
    Arb t2i = Arb::norm_cdf_inv(t2);

    Arb pdf_t1i = Arb::norm_pdf(t1i);
    Arb x = -t1i * pdf_t1i;
    Arb y = pdf_t1i;

    Arb z = Arb::safe_sqrt((1 - rho) / (1 + rho)) * t1i;
    Arb zz = (t2i + rho * t1i) / Arb::safe_sqrt(1 - rho.sqr());

    Arb coef_1a = -2 * x * Arb::norm_cdf(z);
    Arb coef_1b = -2 * y * Arb::safe_sqrt((1 - rho) / (1 + rho)) * Arb::norm_pdf(z);
    Arb coef_2 = 0;
    Arb coef_3 = x * Arb::norm_cdf(zz);
    Arb coef_4 = -x * Arb::norm_cdf(-zz);
    Arb coef_34 = y * rho / Arb::safe_sqrt(1 - rho.sqr()) * Arb::norm_pdf(zz);
    Arb coef_5 = 0;
    Arb coef_6 = x;

    return (coef_1a + coef_1b) * q1 + coef_2 * q2 + (coef_3 + coef_34) * q3
        + (coef_4 + coef_34)* q4 + coef_5 * q5 + coef_6 * q6;
}

Arb prob_d_t2_d_t2(const Arb &t1, const Arb &t2, const Arb &b) {
    Arb rho = - (1 - b) / (1 + b);
    Arb t = (1-b) / 2;
    Arb r12 = 1 - Arb::norm_cdf(Arb::safe_sqrt((1 - rho) / (1 + rho)) * Arb::norm_cdf_inv(t));
    Arb s1 = r12 * (1 - biv_norm_cdf_norm_thresh(t, t, rho)) +
        (1 - r12) * (1 - biv_norm_cdf_norm_thresh(1-t, 1-t, rho));
    Arb s2 = 1;
    Arb diff = s2 - s1;
    Arb p3 = (diff * 2) / (1 + diff * 2);
    Arb p1 = r12 * (1 - p3);
    Arb p2 = 1 - p1 - p3;

    Arb q1 = p3 / 2;
    Arb q2 = q1;
    Arb q3 = p2 - p3/2;
    Arb q4 = p1 - p3/2;
    Arb q5 = q1;
    Arb q6 = q1;

    Arb t1i = Arb::norm_cdf_inv(t1);
    Arb t2i = Arb::norm_cdf_inv(t2);

    Arb pdf_t2i = Arb::norm_pdf(t2i);
    Arb x = -t2i * pdf_t2i;
    Arb y = pdf_t2i;

    Arb z = Arb::safe_sqrt((1 - rho) / (1 + rho)) * t2i;
    Arb zz = (t1i + rho * t2i) / Arb::safe_sqrt(1 - rho.sqr());

    Arb coef_1 = 0;
    Arb coef_2a = -2 * x * Arb::norm_cdf(z);
    Arb coef_2b = -2 * y * Arb::safe_sqrt((1 - rho) / (1 + rho)) * Arb::norm_pdf(z);
    Arb coef_3 = -x * Arb::norm_cdf(-zz);
    Arb coef_4 = x * Arb::norm_cdf(zz);
    Arb coef_34 = y * rho / Arb::safe_sqrt(1 - rho.sqr()) * Arb::norm_pdf(zz);
    Arb coef_5 = x;
    Arb coef_6 = 0;

    return coef_1 * q1 + (coef_2a + coef_2b) * q2 + (coef_3 + coef_34) * q3
        + (coef_4 + coef_34)* q4 + coef_5 * q5 + coef_6 * q6;
}


Arb prob_d_t1_d_t2(const Arb &t1, const Arb &t2, const Arb &b) {
    Arb rho = - (1 - b) / (1 + b);
    Arb t = (1-b) / 2;
    Arb r12 = 1 - Arb::norm_cdf(Arb::safe_sqrt((1 - rho) / (1 + rho)) * Arb::norm_cdf_inv(t));
    Arb s1 = r12 * (1 - biv_norm_cdf_norm_thresh(t, t, rho)) +
        (1 - r12) * (1 - biv_norm_cdf_norm_thresh(1-t, 1-t, rho));
    Arb s2 = 1;
    Arb diff = s2 - s1;
    Arb p3 = (diff * 2) / (1 + diff * 2);
    Arb p1 = r12 * (1 - p3);
    Arb p2 = 1 - p1 - p3;

    Arb q1 = p3 / 2;
    Arb q2 = q1;
    Arb q3 = p2 - p3/2;
    Arb q4 = p1 - p3/2;
    Arb q5 = q1;
    Arb q6 = q1;

    Arb t1i = Arb::norm_cdf_inv(t1);
    Arb t2i = Arb::norm_cdf_inv(t2);

    Arb pdf_t1i = Arb::norm_pdf(t1i);
    Arb x = -t1i * pdf_t1i;
    Arb y = pdf_t1i;
    Arb zz = (t2i + rho * t1i) / Arb::safe_sqrt(1 - rho.sqr());

    return ((q3 + q4) / Arb::safe_sqrt(1-rho.sqr())) * y * Arb::norm_pdf(zz);
}

int check(const Arb &t1, const Arb &t2, const Arb &b) {
    if (Arb::abs(t1 - (1-b)/2) < TYPE_4_HARD_EPS_ALT &&
        Arb::abs(t2 - (1+b)/2) < TYPE_4_HARD_EPS_ALT) {
        Arb d11 = prob_d_t1_d_t1(t1, t2, b);
        Arb d22 = prob_d_t2_d_t2(t1, t2, b);
        Arb d12 = prob_d_t1_d_t2(t1, t2, b);

        d11.println();
        (d11 * d22 - d12 * d12).println();
        
        assert(!d11.is_nan() && !d22.is_nan() && !d12.is_nan());
        if (d11 < 0 && (d11 * d22 - d12 * d12) > 0) {
            return 1;
        }
    }

    if (t1 < TYPE_4_HARD_EPS_ALT2 && t2 < TYPE_4_HARD_EPS_ALT2) {
        return 1;
    }

    if (1-t1 < TYPE_4_HARD_EPS_ALT2 && 1-t2 < TYPE_4_HARD_EPS_ALT2) {
        return 1;
    }

    if (prob(t1, t2, b) < prob((1-b)/2, (1+b)/2, b)) {
        return 1;
    }
    
    if (t1 > 0 && t1 < 1 && t2 > 0 && t2 < 1) {
        // derivative checks
        Arb d_t1 = prob_d_t1(t1, t2, b);
        Arb d_t2 = prob_d_t2(t1, t2, b);

        if (!d_t1.is_nan() && (d_t1 > 0 || d_t1 < 0)) {
            return 1;
        }
        else if (!d_t2.is_nan() && (d_t2 > 0 || d_t2 < 0)) {
            return 1;
        }
    }

    // otherwise we need to split
    Arb rt1 = t1.rad();
    Arb rt2 = t2.rad();
    Arb rb = b.rad();

    /*flint_printf("AT: \n");
    t1.println();
    t2.println();
    b.println();
    prob(t1, t2, b).println();
    prob_d_t1(t1, t2, b).println();
    prob_d_t2(t1, t2, b).println();*/
    

    if (rt1 >= rt2 && rt1 >= rb) {
#ifdef DEBUG
        flint_printf("SPLIT: t1\n");
#endif
        return check(t1.left_half(), t2, b) && check(t1.right_half(), t2, b);
    }
    else if (rt2 >= rb) {
#ifdef DEBUG
        flint_printf("SPLIT: t2\n");
#endif
        return check(t1, t2.left_half(), b) &&
                        check(t1, t2.right_half(), b);
    }
    else {
        return check(t1, t2, b.left_half()) &&
                        check(t1, t2, b.right_half());

    }
}

int main(int argc, char* argv[]) {
    
    flint_set_num_threads(NUM_THREADS);
    
    //(1/(1-eval_low(0.1489, 0.1489))).println();
    //(1/(1-eval_low(-0.1489, -0.1489))).println();
    /* Arb t1(.2);
    Arb t2(.3);
    eval_norm(0, 0).println();
    eval_norm(TYPE_4_T1, TYPE_4_T2).println();
    eval_norm(1, 1).println();*/

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

    Arb t1_range(0, 1);
    Arb t2_range(0, 1);
    Arb b_range(TYPE_4_B1_HARD-2*TYPE_4_HARD_EPS,
                TYPE_4_B1_HARD+2*TYPE_4_HARD_EPS);

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

    flint_printf("Type 4, hardness\n");

    /*Arb t1 =  Arb::norm_cdf(0.9);
    Arb t2 = Arb::norm_cdf(-2.3);

    t1.println();
    t2.println();
    prob(t1, t2, TYPE_4_B1_HARD).println();
    prob_d_t1(t1, t2, TYPE_4_B1_HARD).println();
    prob_d_t2(t1, t2, TYPE_4_B1_HARD).println();
    prob_d_t1_d_t1(t1, t2, TYPE_4_B1_HARD).println();
    prob_d_t2_d_t2(t1, t2, TYPE_4_B1_HARD).println();
    prob_d_t1_d_t2(t1, t2, TYPE_4_B1_HARD).println();*/

    //prob(0,0,TYPE_4_B1_HARD).println();
    flint_printf("RESULT: %d\n", check(t2_range, t2_range, b_range));
    
    flint_cleanup_master();

    return 0;
}
