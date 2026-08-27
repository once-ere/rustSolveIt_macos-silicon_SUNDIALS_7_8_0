/* tides_spin_migration_c.c — C reference for the REBOUNDx port test.
 *
 * reboundx/examples/tides_spin_migration_driven_obliquity_tides/problem.c
 * with the same three verification changes as the other harnesses: no
 * <unistd.h>, no system("rm"), and a final raw IEEE-754 bit dump.
 * Exercises tides_spin + modify_orbits_forces (TWO simultaneous REBOUNDx
 * forces), and mid-run parameter mutation: migration is switched off at
 * tmax/2 by setting tau_a = INFINITY, so it also checks that changing a
 * parameter between integrate() calls behaves identically.
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

double tmax = 100 * 2 * M_PI;

int main(int argc, char* argv[]){
    if (argc > 1) tmax = atof(argv[1]);
    struct reb_simulation* sim = reb_simulation_create();

    const double solar_mass = 1.;
    const double solar_rad = 0.00465;
    reb_simulation_add_fmt(sim, "m r", solar_mass, solar_rad);

    const double p1_mass = 5. * 3.0e-6;
    const double p1_rad = 2.5 * 4.26e-5;
    reb_simulation_add_fmt(sim, "m a e r inc Omega pomega M", p1_mass, 0.17308688, 0.01, p1_rad, 0.5 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.));

    const double p2_mass = 5. * 3.0e-6;
    const double p2_rad = 2.5 * 4.26e-5;
    reb_simulation_add_fmt(sim, "m a e r inc Omega pomega M", p2_mass, 0.23290608, 0.01, p2_rad, -0.431 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.));

    sim->N_active = 3;
    reb_simulation_set_integrator(sim, "whfast");
    sim->dt = 1e-3;

    struct rebx_extras* rebx = rebx_attach(sim);
    struct rebx_force* effect = rebx_load_force(rebx, "tides_spin");
    rebx_add_force(rebx, effect);

    const double solar_spin_period = 20 * 2 * M_PI / 365;
    const double solar_spin = (2 * M_PI) / solar_spin_period;
    const double solar_Q = 1000000.;
    rebx_set_param_double(rebx, &sim->particles[0].ap, "k2", 0.07);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "I", 0.07 * solar_mass * solar_rad * solar_rad);
    rebx_set_param_vec3d(rebx, &sim->particles[0].ap, "Omega", (struct reb_vec3d){.z = solar_spin});

    struct reb_orbit orb = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
    rebx_set_param_double(rebx, &sim->particles[0].ap, "tau", 1./(2.*orb.n*solar_Q));

    const double spin_period_1 = 5. * 2. * M_PI / 365.;
    const double spin_1 = (2. * M_PI) / spin_period_1;
    const double planet_Q = 10000.;
    rebx_set_param_double(rebx, &sim->particles[1].ap, "k2", 0.4);
    rebx_set_param_double(rebx, &sim->particles[1].ap, "I", 0.25 * p1_mass * p1_rad * p1_rad);
    rebx_set_param_vec3d(rebx, &sim->particles[1].ap, "Omega", (struct reb_vec3d){.y=spin_1 * -0.0261769, .z=spin_1 * 0.99965732});
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau", 1./(2.*orb.n*planet_Q));

    double spin_period_2 = 3. * 2. * M_PI / 365.;
    double spin_2 = (2. * M_PI) / spin_period_2;
    rebx_set_param_double(rebx, &sim->particles[2].ap, "k2", 0.4);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "I", 0.25 * p2_mass * p2_rad * p2_rad);
    rebx_set_param_vec3d(rebx, &sim->particles[2].ap, "Omega", (struct reb_vec3d){.y=spin_2 * 0.0249736, .z=spin_2 * 0.99968811});

    struct reb_orbit orb2 = reb_orbit_from_particle(sim->G, sim->particles[2], sim->particles[0]);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau", 1./(2.*orb2.n*planet_Q));

    struct rebx_force* mo = rebx_load_force(rebx, "modify_orbits_forces");
    rebx_add_force(rebx, mo);

    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", -5e6 * 2 * M_PI);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", (-5e6 * 2 * M_PI) / 1.1);

    reb_simulation_move_to_com(sim);

    struct reb_vec3d newz = reb_vec3d_add(reb_simulation_angular_momentum(sim), rebx_tools_spin_angular_momentum(rebx));
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    rebx_simulation_irotate(rebx, rot);
    rebx_spin_initialize_ode(rebx, effect);

    reb_simulation_integrate(sim, tmax/2);

    /* Migration switching off */
    rebx_set_param_double(rebx, &sim->particles[1].ap, "tau_a", INFINITY);
    rebx_set_param_double(rebx, &sim->particles[2].ap, "tau_a", INFINITY);

    reb_simulation_integrate(sim, tmax);

    FILE* of = fopen("state_migration_c.txt","wb");
    fprintf(of, "example tides_spin_migration_driven_obliquity_tides tmax %016llx\n", bits(tmax));
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
    printf("migration done t=%.17e\n", sim->t);

    rebx_free(rebx);
    reb_simulation_free(sim);
    return 0;
}
