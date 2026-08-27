/* kepler_rectilinear_c.c — does the ORIGINAL C Kepler solver terminate for
 * (near-)rectilinear, i.e. (near-)zero angular momentum, hyperbolic motion?
 *
 * The Rust port's new test suite flagged this regime as one where the
 * quartic/Newton iteration and its bisection fallback are most strained,
 * and claimed that h == 0 makes the bisection bounds NaN so the loop
 * cannot terminate. Before recording that as "inherited from the C" we
 * have to check what the C actually does. This program calls the C solver
 * with exactly that input and prints markers before and after, so a hang
 * is unambiguous, plus the final state as raw bits so the Rust twin can be
 * compared against it exactly.
 *
 *   argv[1] (optional) = vy, the transverse velocity.
 *     vy = 0     -> h = 0 exactly (purely radial)
 *     vy = 1e-12 -> h tiny but non-zero (the "near-rectilinear" case)
 *
 * Part of the rebound_rs port verification. GPL-3.0-or-later.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rebound.h"
#include "integrator_whfast.h"

static unsigned long long bits(double x){
    unsigned long long u; memcpy(&u,&x,8); return u;
}

int main(int argc, char* argv[]){
    double vy = (argc > 1) ? atof(argv[1]) : 0.0;
    int nsteps  = (argc > 2) ? atoi(argv[2]) : 20;
    double dt   = (argc > 3) ? atof(argv[3]) : 0.1;

    /* r = (1,0,0), v = (3,vy,0), mu = 1.
       h = r x v = (0,0,vy).  v^2 = 9 > 2*mu/r = 2, so hyperbolic. */
    struct reb_particle p = {0};
    p.x = 1.; p.y = 0.; p.z = 0.;
    p.vx = 3.; p.vy = vy; p.vz = 0.;

    printf("BEFORE: 20 kepler steps, vy=%.17e\n", vy);
    fflush(stdout);

    for (int k = 0; k < nsteps; k++){
        reb_integrator_whfast_kepler_solver(&p, 1.0 /* mu */, dt, NULL);
    }

    printf("AFTER: x=%016llx y=%016llx vx=%016llx vy=%016llx\n",
           bits(p.x), bits(p.y), bits(p.vx), bits(p.vy));
    fflush(stdout);
    return 0;
}
