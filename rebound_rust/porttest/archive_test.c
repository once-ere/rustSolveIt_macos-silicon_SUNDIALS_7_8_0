/* archive_test.c — C side of the Simulationarchive cross-language
 * round-trip verification.
 * Usage: archive_test <integrator> write     — run 3x100 steps, saving
 *          a snapshot after each 100 to archive_c_<integrator>.bin;
 *          dump the final (300 step) state bits to archive_state_c.txt.
 *        archive_test <integrator> continue  — load snapshot 1 (the
 *          200-step state) from archive_rust_<integrator>.bin, run 100
 *          more steps, dump state bits to archive_state_c.txt.
 * Part of the rebound_rs port verification. GPL-3.0-or-later. */
#include "rebound.h"
#include "integrator_whfast.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned long long bits(double x){
    unsigned long long u; memcpy(&u,&x,8); return u;
}

static void setup_particles(struct reb_simulation* r){
    struct reb_particle star = {0};
    star.m = 1.0;
    reb_simulation_add(r, star);
    struct reb_particle planet = {0};
    planet.m = 1e-3; planet.x = 1.6; planet.vy = 0.5;
    reb_simulation_add(r, planet);
    struct reb_particle moon = {0};
    moon.m = 1e-7; moon.x = 1.7; moon.vy = 0.6; moon.z = 0.01; moon.vz = 0.001;
    reb_simulation_add(r, moon);
}

static void configure(struct reb_simulation* r, const char* integrator, void* state){
    if (strcmp(integrator,"whfast-usafe")==0){
        struct reb_integrator_whfast_state* wh = state;
        wh->safe_mode = 0;
    }
    (void)r;
}

static void dump_state(struct reb_simulation* r, const char* integrator){
    FILE* f = fopen("archive_state_c.txt","wb");
    fprintf(f, "integrator %s\n", integrator);
    fprintf(f, "t %016llx\n", bits(r->t));
    fprintf(f, "dt %016llx\n", bits(r->dt));
    for (size_t i=0;i<r->N;i++){
        struct reb_particle p = r->particles[i];
        fprintf(f, "%llu %016llx %016llx %016llx %016llx %016llx %016llx\n",
                (unsigned long long)i,
                bits(p.x), bits(p.y), bits(p.z),
                bits(p.vx), bits(p.vy), bits(p.vz));
    }
    fclose(f);
}

int main(int argc, char* argv[]){
    const char* integrator = argc>1 ? argv[1] : "whfast-usafe";
    const char* mode = argc>2 ? argv[2] : "write";

    const char* real_integrator = integrator;
    if (strncmp(integrator, "whfast", 6)==0) real_integrator = "whfast";

    char fname[256];

    if (strcmp(mode,"load")==0){
        /* Load snapshot 0 from the file given in argv[3], dump state. */
        const char* lfname = argc>3 ? argv[3] : "served.bin";
        struct reb_simulation* r = reb_simulation_create_from_file((char*)lfname, 0);
        if (!r){
            printf("Failed to load %s\n", lfname);
            return 1;
        }
        dump_state(r, integrator);
        printf("load done: t=%e\n", r->t);
        return 0;
    }
    if (strcmp(mode,"write")==0){
        struct reb_simulation* r = reb_simulation_create();
        void* state = reb_simulation_set_integrator(r, real_integrator);
        configure(r, integrator, state);
        r->G = 1.0;
        r->dt = 0.01;
        setup_particles(r);
        sprintf(fname, "archive_c_%s.bin", integrator);
        remove(fname);
        for (int s=0;s<3;s++){
            reb_simulation_steps(r, 100);
            reb_simulation_save_to_file(r, fname);
        }
        dump_state(r, integrator);
        printf("write done: t=%e\n", r->t);
    }else{
        sprintf(fname, "archive_rust_%s.bin", integrator);
        struct reb_simulation* r = reb_simulation_create_from_file(fname, 1);
        if (!r){
            printf("Failed to load %s\n", fname);
            return 1;
        }
        reb_simulation_steps(r, 100);
        dump_state(r, integrator);
        printf("continue done: t=%e\n", r->t);
    }
    return 0;
}
