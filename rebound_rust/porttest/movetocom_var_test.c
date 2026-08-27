/* movetocom_var_test.c — probe reb_simulation_move_to_com's 1st-order
 * variational block (dm accumulator array question). */
#include "rebound.h"
#include <stdio.h>
#include <string.h>

static unsigned long long bits(double x){
    unsigned long long u; memcpy(&u,&x,8); return u;
}

int main(void){
    struct reb_simulation* r = reb_simulation_create();
    r->G = 1.;
    struct reb_particle sun = {0};
    sun.m = 1.;
    reb_simulation_add(r, sun);
    struct reb_particle jup = {0};
    jup.m = 0.000954588;
    jup.x = 5.2;
    jup.vy = 0.4396;
    reb_simulation_add(r, jup);

    reb_simulation_init_megno_seed(r, 12345);

    struct reb_particle com = reb_simulation_com(r);
    printf("com.m = %.17e (bits %016llx)\n", com.m, bits(com.m));
    printf("com.x = %.17e\n", com.x);
    printf("N=%zu N_var=%zu N_var_config=%zu index=%zu\n",
           r->N, r->N_var, r->N_var_config, r->var_config[0].index);
    for (size_t i=0;i<r->N_var;i++){
        printf("BEFORE var[%zu] m=%.17e x=%.17e vx=%.17e\n", i,
               r->particles_var[i].m, r->particles_var[i].x, r->particles_var[i].vx);
    }
    /* what dm would be under each reading */
    double dm_real=0., dm_var=0.;
    for (size_t i=0;i<r->N;i++){ dm_real += r->particles[i].m; dm_var += r->particles_var[i].m; }
    printf("dm(particles)=%.17e  dm(particles_var)=%.17e\n", dm_real, dm_var);

    reb_simulation_move_to_com(r);

    for (size_t i=0;i<r->N_var;i++){
        printf("AFTER  var[%zu] x=%.17e (bits %016llx) vx=%.17e (bits %016llx)\n", i,
               r->particles_var[i].x, bits(r->particles_var[i].x),
               r->particles_var[i].vx, bits(r->particles_var[i].vx));
    }
    reb_simulation_free(r);
    return 0;
}
