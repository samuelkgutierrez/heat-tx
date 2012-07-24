/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* author: samuel k. gutierrez */

/* what we are solving
 *
 * u_t = c * (u_xx * u_yy), 0 <= x,y <= NX, t >= 0
 */

/* just for fun. */

/* http://thecodedecanter.wordpress.com/2010/04/30/modelling-the-2d-heat-equation-in-f-using-100-lines-of-code/ */
/* http://www.cosy.sbg.ac.at/events/parnum05/book/horak1.pdf */

/*
 * double[n+2][n+2] u_old, u_new;
 * double c, delta_t, delta_s;
 * initialize u_old, u_new with initial values and boundary conditions;
 * while (still time points to compute) {
 *     for (int i = 1; i <= n; i++) {
 *         for (int j = 1; j <= n; j++) {
 *             u_new[i, j] = u_old[i, j] + c * delta_t / delta_s^2 *
 *             (u_old[i + 1, j] + u_old[i - 1, j] - 4 * u_old[i, j] +
 *             u_old[i, j + 1] + u_old[i, j - 1]);
 *         }
 *     }
 *     u_old = u_new;
 * }
 */

/*
 * NOTES
 * see also: Crank–Nicolson method
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

/* max simulation time */
#define T_MAX 4096
/* nx and ny */
#define N 256
/* thermal conductivity */
#define THERM_COND 1.6
/* some constant */
#define K 0.1

/* return codes */
enum {
    SUCCESS = 0,
    FAILURE,
    FAILURE_OOR,
    FAILURE_IO,
    FAILURE_INVALID_ARG
};

typedef struct mesh_t {
    int nx;
    int ny;
    double **vals;
} mesh_t;

/* simulation parameters */
typedef struct simulation_params_t {
    /* mesh size (xy) */
    int nx;
    int ny;
    /* thermal conductivity */
    double c;
    double delta_s;
    /* time interval */
    double delta_t;
    /* max simulation time */
    int max_t;
} simulation_params_t;

typedef struct simulation_t {
    simulation_params_t *params;
    mesh_t *u_old;
    mesh_t *u_new;
} simulation_t;

static char *app_name = "heat-tx";
static char *app_ver = "0.1";

/* static forward declarations */
static int
mesh_construct(mesh_t **new_mesh,
               int x /* number of rows */,
               int y /* number of columns */);

static int
mesh_destruct(mesh_t *mesh);

static int
params_construct(simulation_params_t **params);

static int
params_destruct(simulation_params_t **params);

static int
dump_mesh(const mesh_t *mesh);

static double
max_val(const mesh_t *mesh);

/* ////////////////////////////////////////////////////////////////////////// */
static int
sim_param_cp(const simulation_params_t *from,
             simulation_params_t *to)
{
    (void)memcpy(to, from, sizeof(*from));
    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
params_construct(simulation_params_t **params)
{
    simulation_params_t *tmp = NULL;

    if (NULL == params) return FAILURE_INVALID_ARG;

    tmp = (simulation_params_t *)calloc(1, sizeof(simulation_params_t));
    if (NULL == tmp) {
        fprintf(stderr, "out of resources @ %s:%d\n", __FILE__, __LINE__);
        return FAILURE_OOR;
    }
    *params = tmp;

    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
params_destruct(simulation_params_t **params)
{
    if (NULL == params) return FAILURE_INVALID_ARG;

    if (NULL != *params) {
        free(*params);
        *params = NULL;
    }

    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
mesh_construct(mesh_t **new_mesh,
               int x /* number of rows */,
               int y /* number of columns */)
{
    mesh_t *tmp_mesh = NULL;
    int i;

    if (NULL == new_mesh) return FAILURE_INVALID_ARG;

    tmp_mesh = (mesh_t *)calloc(1, sizeof(*tmp_mesh));
    if (NULL == tmp_mesh) {
        fprintf(stderr, "out of resources @ %s:%d\n", __FILE__, __LINE__);
        /* just bail */
        return FAILURE_OOR;
    }
    tmp_mesh->vals = (double **)calloc(x, sizeof(double *));
    if (NULL == tmp_mesh->vals) {
        fprintf(stderr, "out of resources @ %s:%d\n", __FILE__, __LINE__);
        goto error;
    }
    for (i = 0; i < y; ++i) {
        tmp_mesh->vals[i] = (double *)calloc(y, sizeof(double));
        if (NULL == tmp_mesh->vals[i]) {
            fprintf(stderr, "out of resources @ %s:%d\n", __FILE__, __LINE__);
            goto error;
        }
    }
    tmp_mesh->nx = x;
    tmp_mesh->ny = y;

    *new_mesh = tmp_mesh;
    return SUCCESS;

error:
    mesh_destruct(tmp_mesh);
    return FAILURE_OOR;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
mesh_destruct(mesh_t *mesh)
{
    int i;

    if (NULL == mesh) return FAILURE_INVALID_ARG;

    if (NULL != mesh->vals) {
        for (i = 0; i < mesh->ny; ++i) {
            if (NULL != mesh->vals[i]) {
                free(mesh->vals[i]);
                mesh->vals[i] = NULL;
            }
        }
        free(mesh->vals);
        mesh->vals = NULL;
    }
    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
gen_meshes(simulation_t *sim)
{
    int rc = FAILURE;

    printf("    o constructing %d x %d mesh...", sim->params->nx,
                                                 sim->params->ny);

    if (SUCCESS != (rc = mesh_construct(&sim->u_old,
                                        sim->params->nx,
                                        sim->params->ny))) {
        fprintf(stderr, "\nmesh_construct failure @ %s:%d\n", __FILE__,
                __LINE__);
        goto out;
    }
    if (SUCCESS != (rc = mesh_construct(&sim->u_new,
                                        sim->params->nx,
                                        sim->params->ny))) {
        fprintf(stderr, "\nmesh_construct failure @ %s:%d\n", __FILE__,
                __LINE__);
        goto out;
    }

    printf("done!\n");

out:
    if (SUCCESS != rc) {
        mesh_destruct(sim->u_new);
        mesh_destruct(sim->u_old);
    }
    return rc;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
new_simulation(simulation_t **new_sim,
               const simulation_params_t *params)
{
    simulation_t *sim = NULL;
    int rc = FAILURE;

    if (NULL == new_sim) return FAILURE_INVALID_ARG;
    if (NULL == params) return FAILURE_INVALID_ARG;

    printf("::: initializing simulation... \n");

    sim = (simulation_t *)calloc(1, sizeof(*sim));
    if (NULL == sim) {
        fprintf(stderr, "out of resources @ %s:%d\n", __FILE__, __LINE__);
        return FAILURE_OOR;
    }
    if (SUCCESS != (rc = params_construct(&sim->params))) {
        fprintf(stderr, "params_construct failure @ %s:%d: rc = %d\n", __FILE__,
                __LINE__, rc);
        goto out;
    }
    if (SUCCESS != (rc = sim_param_cp(params, sim->params))) {
        fprintf(stderr, "sim_param_cp failure @ %s:%d: rc = %d\n", __FILE__,
                __LINE__, rc);
        goto out;
    }
    if (SUCCESS != (rc = gen_meshes(sim))) {
        fprintf(stderr, "gen_meshes failure @ %s:%d: rc = %d\n", __FILE__,
                __LINE__, rc);
        /* on failure, gen_meshes cleans up after itself */
        goto out;
    }

    *new_sim = sim;

out:
    return rc;

}

/* ////////////////////////////////////////////////////////////////////////// */
static int
init_params(simulation_params_t *params,
            int n,
            double c,
            int max_t)
{

    if (NULL == params) return FAILURE_INVALID_ARG;

    printf("    o initializing simulation parameters...\n");

    params->nx = n;
    params->ny = n;
    params->c = c;
    params->max_t = max_t;
    params->delta_s = 1.0 / (double)(n + 1);
    /* we know from theory that we have to obey the restriction:
     * delta_t <= (delta_s)^2/2c. so just make them equal.
     */
#if 0
    As different choices of n lead to different values for ∆s, it is also
    sensible to choose different values for ∆t. In our case, we will always set
    ∆t = (∆s)^2/4c
    order to guarantee stability.
    params->delta_t = ((1.0) / 2.0 * params->c) / 2.0;
    /* seems to be okay VVVV */
    params->delta_t = (1.0 / (double)(pow(n, 2) + (2 * n) + 1.0)) / (2.0 * c);
#endif
    params->delta_t = pow(params->delta_s, 2.0) / (4.0 * params->c);

    printf("      . max_t: %d\n", params->max_t);
    printf("      . nx: %d\n", params->nx);
    printf("      . ny: %d\n", params->ny);
    printf("      . nx x ny: %d\n", (params->nx * params->ny));
    printf("      . c: %lf\n", params->c);
    printf("      . delta_s: %lf\n", params->delta_s);
    printf("      . delta_t: %lf\n", params->delta_t);
    printf("    o done!\n");

    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
dump_image(const simulation_t *sim,
           const char *where)
{
    FILE *imgfp = NULL;
    const char *header = "P2\n#\n";
    int i, j;

    if (NULL == (imgfp = fopen("test.ppm", "w+"))) {
        fprintf(stderr, "fopen failure @ %s:%d\n", __FILE__, __LINE__);
        return FAILURE_IO;
    }
    /* write image header */
    fprintf(imgfp, "%s", header);
    /* write columns rows */
    fprintf(imgfp, "%d %d\n", sim->params->ny - 2, sim->params->nx - 2);
    /* write the max mesh value */
    fprintf(imgfp, "%d\n", (int)max_val(sim->u_new));
    /* write the matrix */
    for (i = 1; i < sim->u_new->nx - 1; ++i) {
        for (j = 1; j < sim->u_new->ny - 1; ++j) {
            fprintf(imgfp, "%d%s", (int)sim->u_new->vals[i][j],
                    (j == sim->u_new->ny - 2) ? "" : " ");
        }
        fprintf(imgfp, "\n");
    }

    fflush(imgfp);
    if (NULL != imgfp) fclose(imgfp);
    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
mesh_cp(const mesh_t *from,
        mesh_t *to)
{
    int i;
    int j;

    for (i = 0; i < from->nx; ++i) {
        for (j = 0; j < from->ny; ++j) {
            to->vals[i][j] = from->vals[i][j];
        }
    }

    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
run_simulation(simulation_t *sim)
{
    int t, i, j;

    printf("::: starting simulation...\n");
    for (t = 0; t < sim->params->max_t; ++t) {
        printf("      starting loop %d of %d\n", t, sim->params->max_t);
        for (i = 1; i < sim->params->nx - 1; ++i) {
            for (j = 1; j < sim->params->ny - 1; ++j) {
                sim->u_new->vals[i][j] =
                    /* boundary conditions */
                    (i <= 1 || j <= 1 || i >= sim->params->nx - 1 ||
                     j >= sim->params->ny - 1) ? 0.0 :
                    (sim->u_old->vals[i][j] +
                    (sim->params->c *
                    sim->params->delta_t /
                    pow(sim->params->delta_s, 2)) *
                    (sim->u_old->vals[i + 1][j] +
                     sim->u_old->vals[i - 1][j] -
                     4.0 * sim->u_old->vals[i][j] +
                     sim->u_old->vals[i][j + 1] +
                     sim->u_old->vals[i][j - 1]));
            }
        }
        /*      from        to        */
        mesh_cp(sim->u_new, sim->u_old);
    }
    printf("::: done!\n");
    printf("::: starting visualization dump...");
    dump_image(sim, NULL);
    printf("done!\n");
    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static double
max_val(const mesh_t *mesh)
{
    int i, j;
    bool first = true;
    double max = 0.0;

    for (i = 1; i < mesh->nx - 1; ++i) {
        for (j = 1; j < mesh->ny - 1; ++j) {
            if (first) {
                max = mesh->vals[i][j];
                first = false;
                continue;
            }
            if (max < mesh->vals[i][j]) {
                max = mesh->vals[i][j];
            }
        }
    }
    return max;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
dump_mesh(const mesh_t *mesh)
{
    int i, j;

    printf("///////////////////////////////////////////////////////////////\n");
    for (i = 0; i < mesh->nx; ++i) {
        for (j = 0; j < mesh->ny; ++j) {
            printf("%e ", mesh->vals[i][j]);
        }
        printf("\n");
    }
    printf("///////////////////////////////////////////////////////////////\n");

    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
dump_sim(const simulation_t *sim)
{
    int i;
    int j;
    for (i = 0; i < sim->params->nx; ++i) {
        for (j = 0; j < sim->params->ny; ++j) {
            printf("%.2e ", sim->u_new->vals[i][j]);
        }
        printf("\n");
    }

    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
static int
set_initial_conds(simulation_t *sim)
{
    int i;

    if (NULL == sim) return FAILURE_INVALID_ARG;

    printf("    o setting initial conditions...");

#if 0
    for (i = 1; i < sim->params->nx - 1; ++i) {
        sim->u_old->vals[i][1] = K;
    }
#endif
#if 0
    sim->u_old->vals[sim->params->nx / 2][sim->params->ny / 2] = K;
#endif
    sim->u_old->vals[2][2] = 100000000.0;
    sim->u_old->vals[sim->params->nx - 2][2] = 100000000.0;
    sim->u_old->vals[sim->params->nx - 2][sim->params->ny - 2] = 100000000.0;
    sim->u_old->vals[sim->params->nx - 2][2] = 100000000.0;
    mesh_cp(sim->u_old, sim->u_new);

    printf("done\n");

    return SUCCESS;
}

/* ////////////////////////////////////////////////////////////////////////// */
int
main(int argc, char **argv)
{
    int rc = FAILURE;
    int erc = EXIT_FAILURE;
    simulation_params_t *params = NULL;
    simulation_t *sim = NULL;

    /* print application banner */
    printf("\n::: %s %s :::\n", app_name, app_ver);

    if (SUCCESS != (rc = params_construct(&params))) {
        fprintf(stderr, "params_construct failure @ %s:%d: rc = %d\n", __FILE__,
                __LINE__, rc);
        goto cleanup;
    }
    if (SUCCESS != (rc = init_params(params, N, THERM_COND, T_MAX))) {
        fprintf(stderr, "init_params failure @ %s:%d: rc = %d\n", __FILE__,
                __LINE__, rc);
        goto cleanup;
    }
    if (SUCCESS != (rc = new_simulation(&sim, params))) {
        fprintf(stderr, "new_simulation failure @ %s:%d: rc = %d\n", __FILE__,
                __LINE__, rc);
        goto cleanup;
    }
    if (SUCCESS != (rc = set_initial_conds(sim))) {
        fprintf(stderr, "set_initial_conds failure @ %s:%d: rc = %d\n",
                __FILE__, __LINE__, rc);
        goto cleanup;
    }
    if (SUCCESS != (rc = run_simulation(sim))) {
        fprintf(stderr, "run_simulation failure @ %s:%d: rc = %d\n",
                __FILE__, __LINE__, rc);
        goto cleanup;
    }

cleanup:
    params_destruct(&params);
    return erc;
}
