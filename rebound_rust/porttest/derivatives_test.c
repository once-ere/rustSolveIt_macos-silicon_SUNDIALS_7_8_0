/* derivatives_test.c — C reference for the derivatives.c cross-checks.
 * Calls all 65 reb_particle_derivative_* functions for two fixed
 * configurations and dumps the raw IEEE-754 bit patterns of the
 * resulting particle (x,y,z,vx,vy,vz,m) to derivatives_c.txt.
 * Part of the rebound_rs port verification. GPL-3.0-or-later. */
#include "rebound.h"
#include "derivatives.h"
#include <stdio.h>
#include <string.h>

static unsigned long long bits(double x){
    unsigned long long u; memcpy(&u,&x,8); return u;
}

typedef struct reb_particle (*derivfn)(double, struct reb_particle, struct reb_particle);

struct entry { const char* name; derivfn fn; };

static const struct entry fns[] = {
    {"lambda",        reb_particle_derivative_lambda},
    {"h",             reb_particle_derivative_h},
    {"k",             reb_particle_derivative_k},
    {"k_k",           reb_particle_derivative_k_k},
    {"h_h",           reb_particle_derivative_h_h},
    {"lambda_lambda", reb_particle_derivative_lambda_lambda},
    {"k_lambda",      reb_particle_derivative_k_lambda},
    {"h_lambda",      reb_particle_derivative_h_lambda},
    {"k_h",           reb_particle_derivative_k_h},
    {"a",             reb_particle_derivative_a},
    {"a_a",           reb_particle_derivative_a_a},
    {"ix",            reb_particle_derivative_ix},
    {"ix_ix",         reb_particle_derivative_ix_ix},
    {"iy",            reb_particle_derivative_iy},
    {"iy_iy",         reb_particle_derivative_iy_iy},
    {"k_ix",          reb_particle_derivative_k_ix},
    {"h_ix",          reb_particle_derivative_h_ix},
    {"lambda_ix",     reb_particle_derivative_lambda_ix},
    {"lambda_iy",     reb_particle_derivative_lambda_iy},
    {"h_iy",          reb_particle_derivative_h_iy},
    {"k_iy",          reb_particle_derivative_k_iy},
    {"ix_iy",         reb_particle_derivative_ix_iy},
    {"a_ix",          reb_particle_derivative_a_ix},
    {"a_iy",          reb_particle_derivative_a_iy},
    {"a_lambda",      reb_particle_derivative_a_lambda},
    {"a_h",           reb_particle_derivative_a_h},
    {"a_k",           reb_particle_derivative_a_k},
    {"m",             reb_particle_derivative_m},
    {"m_a",           reb_particle_derivative_m_a},
    {"m_lambda",      reb_particle_derivative_m_lambda},
    {"m_h",           reb_particle_derivative_m_h},
    {"m_k",           reb_particle_derivative_m_k},
    {"m_ix",          reb_particle_derivative_m_ix},
    {"m_iy",          reb_particle_derivative_m_iy},
    {"m_m",           reb_particle_derivative_m_m},
    {"e",             reb_particle_derivative_e},
    {"e_e",           reb_particle_derivative_e_e},
    {"inc",           reb_particle_derivative_inc},
    {"inc_inc",       reb_particle_derivative_inc_inc},
    {"Omega",         reb_particle_derivative_Omega},
    {"Omega_Omega",   reb_particle_derivative_Omega_Omega},
    {"omega",         reb_particle_derivative_omega},
    {"omega_omega",   reb_particle_derivative_omega_omega},
    {"f",             reb_particle_derivative_f},
    {"f_f",           reb_particle_derivative_f_f},
    {"a_e",           reb_particle_derivative_a_e},
    {"a_inc",         reb_particle_derivative_a_inc},
    {"a_Omega",       reb_particle_derivative_a_Omega},
    {"a_omega",       reb_particle_derivative_a_omega},
    {"a_f",           reb_particle_derivative_a_f},
    {"e_inc",         reb_particle_derivative_e_inc},
    {"e_Omega",       reb_particle_derivative_e_Omega},
    {"e_omega",       reb_particle_derivative_e_omega},
    {"e_f",           reb_particle_derivative_e_f},
    {"m_e",           reb_particle_derivative_m_e},
    {"inc_Omega",     reb_particle_derivative_inc_Omega},
    {"inc_omega",     reb_particle_derivative_inc_omega},
    {"inc_f",         reb_particle_derivative_inc_f},
    {"m_inc",         reb_particle_derivative_m_inc},
    {"omega_Omega",   reb_particle_derivative_omega_Omega},
    {"Omega_f",       reb_particle_derivative_Omega_f},
    {"m_Omega",       reb_particle_derivative_m_Omega},
    {"omega_f",       reb_particle_derivative_omega_f},
    {"m_omega",       reb_particle_derivative_m_omega},
    {"m_f",           reb_particle_derivative_m_f},
};

int main(void){
    const double G = 1.0;

    struct reb_particle primary1 = {0};
    primary1.m = 1.0;
    primary1.x = 0.1;  primary1.y = -0.2;  primary1.z = 0.05;
    primary1.vx = 0.01; primary1.vy = -0.03; primary1.vz = 0.002;

    struct reb_particle po1 = {0};
    po1.m = 1e-3;
    po1.x = 1.3;  po1.y = 0.4;  po1.z = 0.1;
    po1.vx = -0.2; po1.vy = 0.9; po1.vz = 0.03;

    struct reb_particle primary2 = {0};
    primary2.m = 2.3;

    struct reb_particle po2 = {0};
    po2.m = 1e-5;
    po2.x = 0.7;  po2.y = -0.5;  po2.z = 0.2;
    po2.vx = 0.4; po2.vy = 1.1;  po2.vz = -0.05;

    FILE* f = fopen("derivatives_c.txt","w");
    if (!f) return 1;
    const int n = (int)(sizeof(fns)/sizeof(fns[0]));
    for (int cfg=1; cfg<=2; cfg++){
        struct reb_particle primary = cfg==1 ? primary1 : primary2;
        struct reb_particle po      = cfg==1 ? po1      : po2;
        for (int i=0; i<n; i++){
            struct reb_particle np = fns[i].fn(G, primary, po);
            fprintf(f, "%s cfg%d %016llx %016llx %016llx %016llx %016llx %016llx %016llx\n",
                    fns[i].name, cfg,
                    bits(np.x), bits(np.y), bits(np.z),
                    bits(np.vx), bits(np.vy), bits(np.vz), bits(np.m));
        }
    }
    fclose(f);
    printf("derivatives_c.txt written (%d functions x 2 configs)\n", n);
    return 0;
}
