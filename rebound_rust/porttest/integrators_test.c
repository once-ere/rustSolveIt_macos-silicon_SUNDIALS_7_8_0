/* integrators_test.c — C reference for the integrator cross-checks.
 * Two-body problem with explicit Cartesian initial conditions, fixed
 * particle data, no randomness, no pow() anywhere on the path.
 * Usage: integrators_test <integrator> [order] [steps]
 * Dumps the final state as raw bit patterns to state_c_final.txt.
 * Part of the rebound_rs port verification. GPL-3.0-or-later. */
#include "rebound.h"
#include "integrator_leapfrog.h"
#include "integrator_whfast.h"
#include "integrator_saba.h"
#include "integrator_janus.h"
#include "integrator_eos.h"
#include "integrator_mercurius.h"
#include "integrator_bs.h"
#include "integrator_trace.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned long long bits(double x){
    unsigned long long u; memcpy(&u,&x,8); return u;
}

int main(int argc, char* argv[]){
    const char* integrator = argc>1 ? argv[1] : "ias15";
    unsigned int order = argc>2 ? (unsigned int)atoi(argv[2]) : 2;
    unsigned long long nsteps = argc>3 ? strtoull(argv[3],NULL,10) : 1000;

    struct reb_simulation* r = reb_simulation_create();
    /* whfast configurations are encoded as pseudo names:
     *   whfast        default (jacobi, safe_mode)
     *   whfast-c11    corrector 11        whfast-c17    corrector 17
     *   whfast-dh     democratic heliocentric
     *   whfast-whds   WHDS coordinates
     *   whfast-bary   barycentric coordinates
     *   whfast-mk     modified kick kernel
     *   whfast-comp   composition kernel
     *   whfast-lazy   lazy implementer's kernel
     *   whfast-usafe  safe_mode = 0 */
    const char* real_integrator = integrator;
    if (strncmp(integrator, "whfast", 6)==0) real_integrator = "whfast";
    if (strncmp(integrator, "saba", 4)==0) real_integrator = "saba";
    if (strncmp(integrator, "janus", 5)==0) real_integrator = "janus";
    if (strncmp(integrator, "eos", 3)==0) real_integrator = "eos";
    if (strncmp(integrator, "mercurius", 9)==0) real_integrator = "mercurius";
    if (strncmp(integrator, "bs", 2)==0) real_integrator = "bs";
    if (strncmp(integrator, "trace", 5)==0) real_integrator = "trace";
    void* state = reb_simulation_set_integrator(r, real_integrator);
    if (strcmp(integrator,"leapfrog")==0){
        struct reb_integrator_leapfrog_state* lf = state;
        lf->order = order;
    }
    if (strncmp(integrator, "whfast", 6)==0){
        struct reb_integrator_whfast_state* wh = state;
        if (strcmp(integrator,"whfast-c11")==0)  wh->corrector = 11;
        if (strcmp(integrator,"whfast-c17")==0){ wh->corrector = 17; wh->corrector2 = 1; }
        if (strcmp(integrator,"whfast-dh")==0)   wh->coordinates = REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;
        if (strcmp(integrator,"whfast-whds")==0) wh->coordinates = REB_INTEGRATOR_WHFAST_COORDINATES_WHDS;
        if (strcmp(integrator,"whfast-bary")==0) wh->coordinates = REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC;
        if (strcmp(integrator,"whfast-mk")==0)   wh->kernel = REB_INTEGRATOR_WHFAST_KERNEL_MODIFIEDKICK;
        if (strcmp(integrator,"whfast-comp")==0) wh->kernel = REB_INTEGRATOR_WHFAST_KERNEL_COMPOSITION;
        if (strcmp(integrator,"whfast-lazy")==0) wh->kernel = REB_INTEGRATOR_WHFAST_KERNEL_LAZY;
        if (strcmp(integrator,"whfast-usafe")==0) wh->safe_mode = 0;
    }
    if (strncmp(integrator, "saba", 4)==0){
        struct reb_integrator_saba_state* sb = state;
        if (strcmp(integrator,"saba-1")==0)     sb->type = REB_INTEGRATOR_SABA_TYPE_1;
        if (strcmp(integrator,"saba-2")==0)     sb->type = REB_INTEGRATOR_SABA_TYPE_2;
        if (strcmp(integrator,"saba-3")==0)     sb->type = REB_INTEGRATOR_SABA_TYPE_3;
        if (strcmp(integrator,"saba-4")==0)     sb->type = REB_INTEGRATOR_SABA_TYPE_4;
        if (strcmp(integrator,"saba-cm2")==0)   sb->type = REB_INTEGRATOR_SABA_TYPE_CM_2;
        if (strcmp(integrator,"saba-cl2")==0)   sb->type = REB_INTEGRATOR_SABA_TYPE_CL_2;
        if (strcmp(integrator,"saba-104")==0)   sb->type = REB_INTEGRATOR_SABA_TYPE_10_4;
        if (strcmp(integrator,"saba-864")==0)   sb->type = REB_INTEGRATOR_SABA_TYPE_8_6_4;
        if (strcmp(integrator,"saba-h844")==0)  sb->type = REB_INTEGRATOR_SABA_TYPE_H_8_4_4;
        if (strcmp(integrator,"saba-h864")==0)  sb->type = REB_INTEGRATOR_SABA_TYPE_H_8_6_4;
        if (strcmp(integrator,"saba-h1064")==0) sb->type = REB_INTEGRATOR_SABA_TYPE_H_10_6_4;
        if (strcmp(integrator,"saba-usafe")==0) sb->safe_mode = 0;
    }
    if (strncmp(integrator, "janus", 5)==0){
        struct reb_integrator_janus_state* jn = state;
        if (strcmp(integrator,"janus-2")==0)  jn->order = 2;
        if (strcmp(integrator,"janus-4")==0)  jn->order = 4;
        if (strcmp(integrator,"janus-8")==0)  jn->order = 8;
        if (strcmp(integrator,"janus-10")==0) jn->order = 10;
    }
    if (strncmp(integrator, "eos", 3)==0){
        struct reb_integrator_eos_state* es = state;
        /* eos-<phi0>-<phi1> with numeric type ids 0-8, and eos-usafe */
        if (strcmp(integrator,"eos-usafe")==0){ es->safe_mode = 0; }
        else if (strlen(integrator)==7 && integrator[3]=='-' && integrator[5]=='-'){
            es->phi0 = integrator[4]-'0';
            es->phi1 = integrator[6]-'0';
        }
    }
    if (strncmp(integrator, "mercurius", 9)==0){
        struct reb_integrator_mercurius_state* mc = state;
        /* mercurius          default (L_mercury, r_crit_hill=3, safe_mode)
         * mercurius-usafe    safe_mode = 0
         * mercurius-c4       Hernandez C4 switching function
         * mercurius-c5       Hernandez C5 switching function
         * mercurius-inf      infinitely differentiable switching function
         * mercurius-hill01   r_crit_hill = 0.1 (no encounters: pure WH path) */
        if (strcmp(integrator,"mercurius-usafe")==0)  mc->safe_mode = 0;
        if (strcmp(integrator,"mercurius-c4")==0)     mc->L = reb_integrator_mercurius_L_C4;
        if (strcmp(integrator,"mercurius-c5")==0)     mc->L = reb_integrator_mercurius_L_C5;
        if (strcmp(integrator,"mercurius-inf")==0)    mc->L = reb_integrator_mercurius_L_infinity;
        if (strcmp(integrator,"mercurius-hill01")==0) mc->r_crit_hill = 0.1;
    }
    if (strncmp(integrator, "bs", 2)==0){
        struct reb_integrator_bs_state* b = state;
        /* bs         default (eps_abs = eps_rel = 1e-8)
         * bs-tight   eps_abs = eps_rel = 1e-11
         * bs-loose   eps_abs = eps_rel = 1e-6
         * bs-maxdt   max_dt = 0.02 */
        if (strcmp(integrator,"bs-tight")==0){ b->eps_abs = 1e-11; b->eps_rel = 1e-11; }
        if (strcmp(integrator,"bs-loose")==0){ b->eps_abs = 1e-6;  b->eps_rel = 1e-6; }
        if (strcmp(integrator,"bs-maxdt")==0){ b->max_dt = 0.02; }
    }
    if (strncmp(integrator, "trace", 5)==0){
        struct reb_integrator_trace_state* tr = state;
        /* trace           default (FULL_BS peri mode, r_crit_hill=3)
         * trace-pbs       peri_mode = PARTIAL_BS
         * trace-ias15     peri_mode = FULL_IAS15
         * trace-hill1     r_crit_hill = 1
         * trace-perinone  S_peri = switch_peri_none
         * trace-eta001    peri_crit_eta = 0.01 (forces pericenter flags) */
        if (strcmp(integrator,"trace-pbs")==0)      tr->peri_mode = REB_INTEGRATOR_TRACE_PERIMODE_PARTIAL_BS;
        if (strcmp(integrator,"trace-ias15")==0)    tr->peri_mode = REB_INTEGRATOR_TRACE_PERIMODE_FULL_IAS15;
        if (strcmp(integrator,"trace-hill1")==0)    tr->r_crit_hill = 1;
        if (strcmp(integrator,"trace-perinone")==0) tr->S_peri = reb_integrator_trace_switch_peri_none;
        if (strcmp(integrator,"trace-eta001")==0)   tr->peri_crit_eta = 0.01;
    }
    r->G = 1.0;
    r->dt = 0.01;

    struct reb_particle star = {0};
    star.m = 1.0;
    reb_simulation_add(r, star);

    struct reb_particle planet = {0};
    planet.m = 1e-3;
    planet.x = 1.6;             /* apocenter of a=1, e=0.6 orbit */
    planet.vy = 0.5;            /* roughly the apocenter speed   */
    reb_simulation_add(r, planet);

    struct reb_particle moon = {0};
    moon.m = 1e-7;
    moon.x = 1.7;
    moon.vy = 0.6;
    moon.z = 0.01;
    moon.vz = 0.001;
    reb_simulation_add(r, moon);

    reb_simulation_steps(r, nsteps);

    FILE* f = fopen("state_c_final.txt","wb");
    fprintf(f, "integrator %s order %u steps %llu\n", integrator, order, nsteps);
    fprintf(f, "t %016llx\n", bits(r->t));
    fprintf(f, "dt %016llx\n", bits(r->dt));
    fprintf(f, "steps_done %llu\n", (unsigned long long)r->steps_done);
    for (size_t i=0;i<r->N;i++){
        struct reb_particle p = r->particles[i];
        fprintf(f, "%llu %016llx %016llx %016llx %016llx %016llx %016llx\n",
            (unsigned long long)i,
            bits(p.x), bits(p.y), bits(p.z),
            bits(p.vx), bits(p.vy), bits(p.vz));
    }
    fclose(f);
    printf("%s done: t=%.17e steps=%llu\n", integrator, r->t, (unsigned long long)r->steps_done);
    reb_simulation_free(r);
    return 0;
}
