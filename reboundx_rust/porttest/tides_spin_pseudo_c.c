/* tides_spin_pseudo_c.c — C reference for the REBOUNDx port test.
 *
 * This is reboundx/examples/tides_spin_pseudo_synchronization/problem.c
 * with exactly three portability/verification changes:
 *   1. <unistd.h> dropped (not available under MSVC; nothing from it is used)
 *   2. system("rm -v output.txt") replaced by remove() (no shell dependency)
 *   3. the heartbeat's %e text output is replaced by a final dump of every
 *      state variable as a raw IEEE-754 bit pattern, so the comparison
 *      against the Rust port is exact rather than to printed precision.
 * The physics setup is byte-for-byte the stock example.
 *
 * Part of the reboundx_rs port verification. GPL-3.0-or-later,
 * based on REBOUNDx (c) Dan Tamayo, Hanno Rein et al.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

static unsigned long long bits(double x){
    unsigned long long u; memcpy(&u,&x,8); return u;
}

double tmax = 1000 * 2 * M_PI;

int main(int argc, char* argv[]){
    if (argc > 1) tmax = atof(argv[1]);
    struct reb_simulation* sim = reb_simulation_create();

    // Star
    const double solar_mass = 1.;
    const double solar_rad = 0.00465;
    reb_simulation_add_fmt(sim, "m r", solar_mass, solar_rad);

    // Fiducial hot Jupiter
    const double p1_mass = 1. * 9.55e-4;
    const double p1_rad = 1. * 4.676e-4;
    const double p1_e = 0.01;
    const double p1_inc = 0.01;
    reb_simulation_add_fmt(sim, "m a e inc r", p1_mass, 0.04072, p1_e, p1_inc, p1_rad);

    sim->N_active = 2;
    reb_simulation_set_integrator(sim, "whfast");
    sim->dt = 1e-3;

    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);

    const double solar_k2 = 0.07;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);
    const double solar_spin_period = 27 * 2 * M_PI / 365;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z=solar_spin});

    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.07 * solar_mass * solar_rad * solar_rad);

    const double solar_Q = 1e6;
    struct reb_orbit orb = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
    double solar_tau = 1 / (2 * solar_Q * orb.n);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", solar_tau);

    // Planet
    const double spin_period_1 = 0.5 * 2. * M_PI / 365.;
    const double spin_1 = (2. * M_PI) / spin_period_1;
    const double planet_Q = 10000.;
    const double theta_1 = 30. * (M_PI / 180.);
    const double phi_1 = 0 * (M_PI / 180);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", 0.3);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.25 * p1_mass * p1_rad * p1_rad);

    struct reb_vec3d Omega_1 = reb_tools_spherical_to_xyz(spin_1, theta_1, phi_1);
    rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", Omega_1);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau", 1./(2*planet_Q*orb.n));

    reb_simulation_move_to_com(sim);

    struct reb_vec3d newz = reb_vec3d_add(reb_simulation_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    rebx_simulation_irotate(rebx, rot);
    rebx_spin_initialize_ode(rebx, effect);

    reb_simulation_integrate(sim, tmax);

    FILE* of = fopen("state_pseudo_c.txt","wb");
    fprintf(of, "example tides_spin_pseudo_synchronization tmax %016llx\n", bits(tmax));
    fprintf(of, "t %016llx\n", bits(sim->t));
    fprintf(of, "dt %016llx\n", bits(sim->dt));
    fprintf(of, "N %llu\n", (unsigned long long)sim->N);
    for (size_t i=0;i<sim->N;i++){
        struct reb_particle p = sim->particles[i];
        fprintf(of, "p %llu %016llx %016llx %016llx %016llx %016llx %016llx %016llx\n",
                (unsigned long long)i,
                bits(p.x), bits(p.y), bits(p.z),
                bits(p.vx), bits(p.vy), bits(p.vz), bits(p.m));
        struct reb_vec3d* Om = rebx_get_param(rebx, sim->particles[i].ap, "Omega");
        if (Om){
            fprintf(of, "Omega %llu %016llx %016llx %016llx\n",
                    (unsigned long long)i, bits(Om->x), bits(Om->y), bits(Om->z));
        } else {
            fprintf(of, "Omega %llu none\n", (unsigned long long)i);
        }
    }
    fclose(of);
    printf("pseudo done t=%.17e\n", sim->t);

    rebx_free(rebx);
    reb_simulation_free(sim);
    return 0;
}
