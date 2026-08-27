/* frequency_test.c — C reference for the frequency_analysis cross-check.
 * Runs reb_frequency_analysis (MFT, FMFT, FMFT2) on a deterministic
 * three-frequency synthetic signal and dumps the raw IEEE-754 bit
 * patterns of the outputs to frequency_c.txt.
 * Part of the rebound_rs port verification. GPL-3.0-or-later. */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static unsigned long long bits(double x){
    unsigned long long u; memcpy(&u,&x,8); return u;
}

int main(void){
    const size_t ndata = 256;
    const size_t nfreq = 3;
    double* input = malloc(sizeof(double)*2*ndata);
    /* Quasi-periodic signal with three frequencies (rad per sample). */
    const double f1 = 0.30, a1 = 1.00, p1 = 0.40;
    const double f2 = 0.55, a2 = 0.35, p2 = 1.90;
    const double f3 = 0.11, a3 = 0.10, p3 = 5.10;
    for (size_t i=0;i<ndata;i++){
        double t = (double)i;
        input[2*i]   = a1*cos(f1*t+p1) + a2*cos(f2*t+p2) + a3*cos(f3*t+p3);
        input[2*i+1] = a1*sin(f1*t+p1) + a2*sin(f2*t+p2) + a3*sin(f3*t+p3);
    }

    FILE* f = fopen("frequency_c.txt","wb");
    int types[3] = {REB_FREQUENCY_ANALYSIS_MFT, REB_FREQUENCY_ANALYSIS_FMFT, REB_FREQUENCY_ANALYSIS_FMFT2};
    const char* names[3] = {"MFT","FMFT","FMFT2"};
    for (int ti=0; ti<3; ti++){
        double output[9] = {0};
        int ret = reb_frequency_analysis(output, nfreq, 0.05, 1.0, types[ti], input, ndata);
        fprintf(f, "%s ret %d\n", names[ti], ret);
        for (size_t k=0;k<3*nfreq;k++){
            fprintf(f, "%zu %016llx\n", k, bits(output[k]));
        }
    }
    fclose(f);
    free(input);
    printf("frequency_test done\n");
    return 0;
}
