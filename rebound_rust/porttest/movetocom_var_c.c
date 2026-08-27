/* movetocom_var_c.c — C reference for reb_simulation_move_to_com with
 * variational particles.
 *
 * The port audit found that the Rust translation summed `particles_var`
 * where the C sums `particles` when building the first-order `dm`
 * accumulator, which silently dropped a whole term of the variational
 * centre-of-mass shift (and therefore changed every MEGNO / Lyapunov
 * result). This program produces the C's answer for a Sun+Jupiter system
 * with MEGNO enabled, as raw IEEE-754 bits, so the Rust twin can be
 * compared against it exactly.
 *
 * Part of the rebound_rs port verification. GPL-3.0-or-later.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rebound.h"

static unsigned long long bits(double x){
    unsigned long long u; memcpy(&u,&x,8); return u;
}

int main(void){
    struct reb_simulation* r = reb_simulation_create();
    r->G = 1.0;
    r->dt = 0.01;
    reb_simulation_set_integrator(r, "ias15");

    struct reb_particle sun = {0};
    sun.m = 1.0;
    reb_simulation_add(r, sun);
    reb_simulation_add_fmt(r, "m a e", 9.54579e-4, 5.2, 0.0489);

    reb_simulation_init_megno_seed(r, 12345);

    struct reb_particle com = reb_simulation_com(r);
    printf("com.m %016llx com.x %016llx\n", bits(com.m), bits(com.x));

    reb_simulation_move_to_com(r);

    FILE* f = fopen("movetocom_var_c.txt","wb");
    fprintf(f, "N %llu N_var %llu\n",
            (unsigned long long)r->N, (unsigned long long)r->N_var);
    for (size_t i=0;i<r->N;i++){
        struct reb_particle p = r->particles[i];
        fprintf(f, "p %llu %016llx %016llx %016llx %016llx %016llx %016llx\n",
                (unsigned long long)i, bits(p.x), bits(p.y), bits(p.z),
                bits(p.vx), bits(p.vy), bits(p.vz));
    }
    for (size_t i=0;i<r->N_var;i++){
        struct reb_particle p = r->particles_var[i];
        fprintf(f, "v %llu %016llx %016llx %016llx %016llx %016llx %016llx\n",
                (unsigned long long)i, bits(p.x), bits(p.y), bits(p.z),
                bits(p.vx), bits(p.vy), bits(p.vz));
    }
    fclose(f);
    printf("movetocom_var_c done\n");
    reb_simulation_free(r);
    return 0;
}
