/* bs_pow_diff.c — C side of the targeted differential that explains the
 * one-ULP divergence of the BS integrator's proposed timestep.
 *
 * The BS step-size controller (integrator_bs.c) computes
 *     exp   = 1.0 / (2*k + 1)
 *     fac   = stepControl2 / pow(error / stepControl1, exp)
 *     power = pow(stepControl3, exp)
 * i.e. pow() with the odd-reciprocal exponents 1/3, 1/5, ... 1/17.
 * This harness evaluates exactly those pow() calls over a wide range of
 * `error` values and dumps raw bit patterns, so the Rust twin can be
 * compared bit-for-bit and the divergence attributed precisely.
 *
 * Part of the rebound_rs port verification. GPL-3.0-or-later.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>

static unsigned long long bits(double x){
    unsigned long long u; memcpy(&u,&x,8); return u;
}

int main(void){
    const double stepControl1 = 0.65;
    const double stepControl2 = 0.94;
    const double stepControl3 = 0.02;
    FILE* f = fopen("bs_pow_c.txt","wb");
    long n = 0;
    for (int k = 1; k <= 8; k++){
        const double e = 1.0 / (2*k + 1);
        /* the constant call: pow(stepControl3, exp) */
        fprintf(f, "P %d %016llx\n", k, bits(pow(stepControl3, e)));
        /* the error-dependent call over 8 decades, 25000 samples per k */
        for (long i = 0; i < 25000; i++){
            double error = 1e-8 * pow(10.0, 8.0 * (double)i / 25000.0);
            double v = pow(error / stepControl1, e);
            fprintf(f, "E %d %ld %016llx\n", k, i, bits(v));
            n++;
        }
    }
    fclose(f);
    printf("bs_pow_c: %ld samples\n", n);
    return 0;
}
