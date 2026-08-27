/* tides_spin_kozai_c.c — C reference for the REBOUNDx port test.
 *
 * reboundx/examples/tides_spin_kozai/problem.c with the same three
 * verification changes as the other harnesses: no <unistd.h>, no
 * system("rm"), and a final raw IEEE-754 bit dump instead of %e text.
 * Exercises tides_spin + gr_potential under the ADAPTIVE IAS15 integrator,
 * so it also checks that the REBOUNDx force feeds IAS15's error estimator
 * identically (the adaptive timestep sequence must match bit-for-bit).
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

double tmax = 1e5;

int main(int argc, char* argv[]){
    if (argc > 1) tmax = atof(argv[1]);
    struct reb_simulation* sim = reb_simulation_create();
    sim->dt             = M_PI*1e-1;
    reb_simulation_set_integrator(sim, "ias15");

    struct reb_particle star = {0};
    star.m  = 1;
    star.r = 0.00465;
    reb_simulation_add(sim, star);

    double planet_m  = 0.054 * 9.55e-4;
    double planet_r = 0.3 * 4.676e-4;
    double planet_a = 2.;
    double planet_e = 0.001;
    reb_simulation_add_fmt(sim, "m r a e", planet_m, planet_r, planet_a, planet_e);

    double perturber_m  = 1;
    double perturber_a = 50.;
    double perturber_e = 0.7 * M_PI / 180.;
    double perturber_inc = 80. * M_PI / 180.;
    reb_simulation_add_fmt(sim, "m a e inc", perturber_m, perturber_a, perturber_e, perturber_inc);

    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);

    const double solar_k2 = 0.07;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", solar_k2);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.07 * star.m * star.r * star.r);

    const double solar_spin_period = 4.6 * 2. * M_PI / 365.;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z=solar_spin});

    const double solar_Q = 1e6;
    struct reb_orbit orb = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
    double solar_tau = 1 / (2 * solar_Q * orb.n);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", solar_tau);

    const double planet_k2 = 0.4;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", planet_k2);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.25 * planet_m * planet_r * planet_r);

    const double spin_period_p = 1. * 2. * M_PI / 365.;
    const double spin_p = (2. * M_PI) / spin_period_p;
    const double theta_p = 0. * M_PI / 180.;
    const double phi_p = 0. * M_PI / 180;
    struct reb_vec3d Omega_sv = reb_tools_spherical_to_xyz(spin_p, theta_p, phi_p);
    rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", Omega_sv);

    const double planet_Q = 3e5;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau", 1./(2.*planet_Q*orb.n));

    struct rebx_force* gr = rebx_load_force(rebx, "gr_potential");
    rebx_add_force(rebx, gr);
    rebx_set_param_double(rebx, &gr->ap, "c", 10065.32);

    reb_simulation_move_to_com(sim);

    struct reb_vec3d newz = reb_vec3d_add(reb_simulation_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    rebx_simulation_irotate(rebx, rot);
    rebx_spin_initialize_ode(rebx, effect);

    reb_simulation_integrate(sim, tmax);

    FILE* of = fopen("state_kozai_c.txt","wb");
    fprintf(of, "example tides_spin_kozai tmax %016llx\n", bits(tmax));
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
    printf("kozai done t=%.17e\n", sim->t);

    rebx_free(rebx);
    reb_simulation_free(sim);
    return 0;
}
