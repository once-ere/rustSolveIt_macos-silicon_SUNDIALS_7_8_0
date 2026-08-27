/* rebx_binary_read_c.c — loads a REBOUNDx binary with the C library and
 * dumps everything the Rust round-trip example checks, as raw IEEE-754
 * bit patterns, so the two readers can be compared. Verification only.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rebound.h"
#include "reboundx.h"

static unsigned long long bits(double x){
    unsigned long long u; memcpy(&u,&x,8); return u;
}

static int len(struct rebx_node* n){ int c=0; while(n){c++;n=n->next;} return c; }

int main(int argc, char* argv[]){
    const char* in = (argc > 1) ? argv[1] : "rebx_c_reference.bin";
    struct reb_simulation* sim = reb_simulation_create();
    struct reb_particle star = {0}; star.m = 1.; star.r = 0.005;
    reb_simulation_add(sim, star);
    struct reb_particle p1 = {0}; p1.m = 1e-3; p1.x = 1.; p1.vy = 1.;
    reb_simulation_add(sim, p1);
    struct reb_particle p2 = {0}; p2.m = 5e-4; p2.x = 2.; p2.vy = 0.7071067811865476;
    reb_simulation_add(sim, p2);

    struct rebx_extras* rebx = rebx_create_extras_from_binary(sim, in);
    if (rebx == NULL){ printf("LOAD FAILED\n"); return 1; }

    printf("registered_params %d\n", len(rebx->registered_params));
    printf("allocated_forces %d\n", len(rebx->allocated_forces));
    printf("allocated_operators %d\n", len(rebx->allocated_operators));
    printf("additional_forces %d\n", len(rebx->additional_forces));
    printf("pre %d post %d\n", len(rebx->pre_timestep_modifications), len(rebx->post_timestep_modifications));

    /* order of additional_forces, head first */
    for (struct rebx_node* n = rebx->additional_forces; n; n = n->next){
        struct rebx_force* f = n->object;
        printf("additional_force %s\n", f->name);
    }
    for (struct rebx_node* n = rebx->post_timestep_modifications; n; n = n->next){
        struct rebx_step* s = n->object;
        printf("post_step %s %016llx\n", s->operator->name, bits(s->dt_fraction));
    }
    for (struct rebx_node* n = rebx->pre_timestep_modifications; n; n = n->next){
        struct rebx_step* s = n->object;
        printf("pre_step %s %016llx\n", s->operator->name, bits(s->dt_fraction));
    }

    struct rebx_force* gr = rebx_get_force(rebx, "gr_potential");
    struct rebx_force* cf = rebx_get_force(rebx, "central_force");
    struct rebx_operator* mm = rebx_get_operator(rebx, "modify_mass");
    printf("gr %s cf %s mm %s\n", gr?"yes":"no", cf?"yes":"no", mm?"yes":"no");

    double* c = rebx_get_param(rebx, gr->ap, "c");
    int* src   = rebx_get_param(rebx, gr->ap, "gr_source");
    double* ac = rebx_get_param(rebx, cf->ap, "Acentral");
    double* cop= rebx_get_param(rebx, mm->ap, "c");
    printf("gr.c %016llx\n", c?bits(*c):0ULL);
    printf("gr.gr_source %d\n", src?*src:-999999);
    printf("cf.Acentral %016llx\n", ac?bits(*ac):0ULL);
    printf("mm.c %016llx\n", cop?bits(*cop):0ULL);

    /* parameter name order for gr, head first */
    for (struct rebx_node* n = gr->ap; n; n = n->next){
        struct rebx_param* p = n->object;
        printf("gr.param %s\n", p->name);
    }

    double* tm = rebx_get_param(rebx, sim->particles[1].ap, "tau_mass");
    int* pr    = rebx_get_param(rebx, sim->particles[1].ap, "primary");
    struct reb_vec3d* om = rebx_get_param(rebx, sim->particles[1].ap, "Omega");
    struct rebx_force* fp = rebx_get_param(rebx, sim->particles[1].ap, "force");
    double* be  = rebx_get_param(rebx, sim->particles[2].ap, "beta");
    int* co     = rebx_get_param(rebx, sim->particles[2].ap, "coordinates");
    printf("p1.tau_mass %016llx\n", tm?bits(*tm):0ULL);
    printf("p1.primary %d\n", pr?*pr:-999999);
    printf("p1.Omega %016llx %016llx %016llx\n", om?bits(om->x):0ULL, om?bits(om->y):0ULL, om?bits(om->z):0ULL);
    printf("p1.force %s\n", fp?fp->name:"NULL");
    printf("p2.beta %016llx\n", be?bits(*be):0ULL);
    printf("p2.coordinates %d\n", co?*co:-999999);
    for (struct rebx_node* n = sim->particles[1].ap; n; n = n->next){
        struct rebx_param* p = n->object;
        printf("p1.param %s\n", p->name);
    }
    printf("p0.params %d\n", len(sim->particles[0].ap));
    return 0;
}
