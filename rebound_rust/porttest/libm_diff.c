/* libm_diff.c — dump f(x) bit patterns for a deterministic corpus, to
 * be compared against the Rust side (libm_diff.rs). Part of the
 * rebound_rs port verification. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

static uint64_t bits(double x){ uint64_t u; memcpy(&u,&x,8); return u; }

/* xorshift64 for a portable deterministic corpus */
static uint64_t s = 88172645463325252ULL;
static uint64_t xs(void){ s^=s<<13; s^=s>>7; s^=s<<17; return s; }

int main(void){
    FILE* f = fopen("libm_c.txt","wb");
    for (int i=0;i<200000;i++){
        /* doubles in (-1000, 1000) and small magnitudes */
        double x = ((double)(xs()%2000000000ULL)/1e6)-1000.0;
        double y = ((double)(xs()%2000000000ULL)/1e6)-1000.0;
        double xp = fabs(x)+1e-9;
        fprintf(f, "%016llx %016llx %016llx %016llx %016llx %016llx %016llx %016llx %016llx\n",
            (unsigned long long)bits(sin(x)),
            (unsigned long long)bits(cos(x)),
            (unsigned long long)bits(tan(x)),
            (unsigned long long)bits(atan2(y,x)),
            (unsigned long long)bits(pow(xp, -0.234)),
            (unsigned long long)bits(sqrt(xp)),
            (unsigned long long)bits(fmod(y, 3.7)),
            (unsigned long long)bits(exp(x/100.)),
            (unsigned long long)bits(log(xp)));
    }
    fclose(f);
    printf("done\n");
    return 0;
}
/* appended: exp/log differential */
