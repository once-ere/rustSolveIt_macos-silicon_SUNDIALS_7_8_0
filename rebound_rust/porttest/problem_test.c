/**
 * shearing_sheet port test — C reference side.
 *
 * The stock shearing_sheet/problem.c with three controlled changes so
 * that the run is exactly reproducible and comparable:
 *   1. r->rand_seed = 42 (stock: time+pid seed)
 *   2. no web server, no heartbeat (wall-clock output excluded)
 *   3. run exactly 400 timesteps via reb_simulation_steps(), then dump
 *      every particle's state as raw IEEE-754 bit patterns.
 *
 * Everything else (constants, particle initialization loop, modules,
 * tolerances) is byte-for-byte the stock example.
 *
 * Part of the rebound_rs port verification. GPL-3.0-or-later, based on
 * REBOUND (c) Hanno Rein et al.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* The stock example computes eps = 0.32*pow(fabs(v)*100., -0.234).
 * `pow` is the ONE libm function where Rust's runtime does not defer to
 * the UCRT (differential: 60/200000 inputs differ by <=2 ulp; sin, cos,
 * tan, atan2, sqrt, fmod, exp, log are all bit-identical). For the
 * bit-identity harness the same law is therefore written via exp/log —
 * IDENTICALLY on both sides — so that any remaining divergence would
 * expose a real port defect rather than the documented pow difference.
 * The `pow` form is exercised by the separate stock-run comparison. */
double coefficient_of_restitution_bridges(const struct reb_simulation* const r, double v){
    (void)r;
    // assumes v in units of [m/s]
    double eps = 0.32*exp(-0.234*log(fabs(v)*100.));
    if (eps>1) eps=1;
    if (eps<0) eps=0;
    return eps;
}

static unsigned long long bits(double x){
    unsigned long long u;
    memcpy(&u, &x, 8);
    return u;
}

int main(int argc, char* argv[]) {
    unsigned long long nsteps = 400;
    if (argc>1) nsteps = strtoull(argv[1], NULL, 10);
    struct reb_simulation* r = reb_simulation_create();
    r->rand_seed         = 42;                  // CONTROLLED SEED
    r->opening_angle2    = .5;
    reb_simulation_set_integrator(r, "sei");
    r->boundary          = REB_BOUNDARY_SHEAR;
    r->gravity           = REB_GRAVITY_TREE;
    r->collision         = REB_COLLISION_TREE;
    r->collision_resolve = reb_collision_resolve_hardsphere;
    r->OMEGA             = 0.00013143527;       // 1/s
    r->G                 = 6.67428e-11;         // N / (1e-5 kg)^2 m^2
    r->softening         = 0.1;                 // m
    r->dt                = 1e-3*2.*M_PI/r->OMEGA;  // s
    double surfacedensity          = 400;     // kg/m^2
    double particle_density        = 400;     // kg/m^3
    double particle_radius_min     = 1;       // m
    double particle_radius_max     = 4;       // m
    double particle_radius_slope   = -3;
    double root_size             = 100;         // m
    r->root_size = root_size;
    r->N_root_x = 2;
    r->N_root_y = 2;
    r->N_ghost_x = 2;
    r->N_ghost_y = 2;
    r->N_ghost_z = 0;
    struct reb_vec3d boxsize = {
        .x = r->root_size*(double)r->N_root_x,
        .y = r->root_size*(double)r->N_root_y,
        .z = r->root_size*(double)r->N_root_z,
    };

    printf("Toomre wavelength: %f\n",4.*M_PI*M_PI*surfacedensity/r->OMEGA/r->OMEGA*r->G);
    r->coefficient_of_restitution = coefficient_of_restitution_bridges;
    r->minimum_collision_velocity = particle_radius_min*r->OMEGA*0.001;

    // Add all ring particles
    double total_mass = surfacedensity*boxsize.x*boxsize.y;
    double mass = 0;
    while(mass<total_mass){
        struct reb_particle pt = {0};
        pt.x         = reb_random_uniform(r, -boxsize.x/2.,boxsize.x/2.);
        pt.y         = reb_random_uniform(r, -boxsize.y/2.,boxsize.y/2.);
        pt.z         = reb_random_normal(r, 1.);                    // m
        pt.vx        = 0;
        pt.vy        = -1.5*pt.x*r->OMEGA;
        pt.vz        = 0;
        pt.ax        = 0;
        pt.ay        = 0;
        pt.az        = 0;
        double radius     = reb_random_powerlaw(r, particle_radius_min,particle_radius_max,particle_radius_slope);
        pt.r         = radius;                        // m
        double        particle_mass = particle_density*4./3.*M_PI*radius*radius*radius;
        pt.m         = particle_mass;     // kg
        reb_simulation_add(r, pt);
        mass += particle_mass;
    }

    printf("N after init: %llu\n", (unsigned long long)r->N);

    // Dump the initial conditions too (verifies the RNG stream match).
    FILE* f0 = fopen("state_c_init.txt","wb");
    fprintf(f0, "N %llu\n", (unsigned long long)r->N);
    for (size_t i=0;i<r->N;i++){
        struct reb_particle p = r->particles[i];
        fprintf(f0, "%llu %016llx %016llx %016llx %016llx %016llx %016llx %016llx %016llx\n",
            (unsigned long long)i,
            bits(p.x), bits(p.y), bits(p.z),
            bits(p.vx), bits(p.vy), bits(p.vz),
            bits(p.m), bits(p.r));
    }
    fclose(f0);

    reb_simulation_steps(r, nsteps);

    FILE* f = fopen("state_c_final.txt","wb");
    fprintf(f, "N %llu\n", (unsigned long long)r->N);
    fprintf(f, "t %016llx %.17e\n", bits(r->t), r->t);
    fprintf(f, "steps_done %llu\n", (unsigned long long)r->steps_done);
    fprintf(f, "collisions_log_n %lld\n", (long long)r->collisions_log_n);
    fprintf(f, "collisions_plog %016llx %.17e\n", bits(r->collisions_plog), r->collisions_plog);
    fprintf(f, "rand_seed %u\n", r->rand_seed);
    for (size_t i=0;i<r->N;i++){
        struct reb_particle p = r->particles[i];
        fprintf(f, "%llu %016llx %016llx %016llx %016llx %016llx %016llx %016llx %016llx\n",
            (unsigned long long)i,
            bits(p.x), bits(p.y), bits(p.z),
            bits(p.vx), bits(p.vy), bits(p.vz),
            bits(p.m), bits(p.r));
    }
    fclose(f);
    printf("final: t=%.17e steps=%llu collisions=%lld\n",
        r->t, (unsigned long long)r->steps_done, (long long)r->collisions_log_n);
    reb_simulation_free(r);
    return 0;
}
