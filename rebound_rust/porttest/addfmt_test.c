/* addfmt_test.c — C reference for reb_simulation_add_fmt cross-check.
 * Adds the built-in solar system, then a particle from orbital
 * elements, then one from Pal coordinates, then one from a period.
 * Dumps all particle bit patterns. GPL-3.0-or-later. */
#include "rebound.h"
#include <stdio.h>
#include <string.h>

static unsigned long long bits(double x){
    unsigned long long u; memcpy(&u,&x,8); return u;
}

int main(void){
    struct reb_simulation* r = reb_simulation_create();
    r->G = 1.0;
    reb_simulation_add_fmt(r, "solar system");
    reb_simulation_add_fmt(r, "m a e inc Omega omega f", 1e-9, 12.5, 0.3, 0.2, 0.6, 1.1, 2.5);
    reb_simulation_add_fmt(r, "m a l h k ix iy", 2e-9, 15.5, 0.7, 0.05, -0.03, 0.01, 0.02);
    reb_simulation_add_fmt(r, "m P e M", 3e-9, 100.0, 0.1, 0.5);

    FILE* f = fopen("addfmt_c.txt","wb");
    for (size_t i=0;i<r->N;i++){
        struct reb_particle p = r->particles[i];
        fprintf(f, "%llu %016llx %016llx %016llx %016llx %016llx %016llx %016llx\n",
                (unsigned long long)i, bits(p.m),
                bits(p.x), bits(p.y), bits(p.z),
                bits(p.vx), bits(p.vy), bits(p.vz));
    }
    fclose(f);
    printf("addfmt done N=%llu\n", (unsigned long long)r->N);
    return 0;
}
