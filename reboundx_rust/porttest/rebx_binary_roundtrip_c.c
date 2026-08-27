/* rebx_binary_roundtrip_c.c — C reference for the REBOUNDx binary
 * format. Builds exactly the state that
 * reboundx_rust/examples/rebx_binary_roundtrip.rs builds and writes it
 * with the C's own rebx_output_binary, so the two files can be diffed
 * byte for byte. Verification scaffolding only.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"

int main(int argc, char* argv[]){
    const char* out = (argc > 1) ? argv[1] : "rebx_c_reference.bin";

    struct reb_simulation* sim = reb_simulation_create();

    struct reb_particle star = {0};
    star.m = 1.; star.r = 0.005;
    reb_simulation_add(sim, star);
    struct reb_particle p1 = {0};
    p1.m = 1e-3; p1.x = 1.; p1.vy = 1.;
    reb_simulation_add(sim, p1);
    struct reb_particle p2 = {0};
    p2.m = 5e-4; p2.x = 2.; p2.vy = 0.65;
    reb_simulation_add(sim, p2);

    struct rebx_extras* rebx = rebx_attach(sim);

    struct rebx_force* gr = rebx_load_force(rebx, "gr_potential");
    rebx_add_force(rebx, gr);
    struct rebx_force* cf = rebx_load_force(rebx, "central_force");
    rebx_add_force(rebx, cf);

    struct rebx_operator* mm = rebx_load_operator(rebx, "modify_mass");
    struct rebx_operator* dr = rebx_load_operator(rebx, "drift");
    rebx_add_operator_step(rebx, mm, 0.5, REBX_TIMING_PRE);
    rebx_add_operator_step(rebx, mm, 0.5, REBX_TIMING_POST);
    rebx_add_operator_step(rebx, dr, 1.0, REBX_TIMING_POST);

    rebx_set_param_double(rebx, &gr->ap, "c", 10065.320005560323);
    rebx_set_param_int(rebx, &gr->ap, "gr_source", 7);
    rebx_set_param_double(rebx, &cf->ap, "Acentral", 1.2345678901234567e-8);
    rebx_set_param_double(rebx, &mm->ap, "c", -3.14159265358979e12);

    struct rebx_node** ap1 = (struct rebx_node**)&sim->particles[1].ap;
    rebx_set_param_double(rebx, ap1, "tau_mass", 1.7976931348623157e30);
    rebx_set_param_int(rebx, ap1, "primary", -12345);
    struct reb_vec3d Om;
    Om.x = 1.5e-7; Om.y = -2.5e13; Om.z = 0.30000000000000004;
    rebx_set_param_vec3d(rebx, ap1, "Omega", Om);
    rebx_set_param_pointer(rebx, ap1, "force", gr);

    struct rebx_node** ap2 = (struct rebx_node**)&sim->particles[2].ap;
    rebx_set_param_double(rebx, ap2, "beta", 0.1 + 0.2);
    rebx_set_param_int(rebx, ap2, "coordinates", 2);

    rebx_output_binary(rebx, (char*)out);
    printf("C wrote %s\n", out);
    return 0;
}
