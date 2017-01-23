/* Standard library includes. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* The ENT API. */
#include "ent.h"

/* Handle the Physics. */
static struct ent_physics * physics = NULL;

/* Error handler: dump any error message and exit to the OS. */
static void handle_error(enum ent_return rc, ent_function_t * caller)
{
        /* Dump an error message. */
        ent_error_print(stderr, rc, caller, "\t", "\n");
        fprintf(stderr, "\n");

        /* Finalise and exit to the OS. */
        ent_physics_destroy(&physics);
        exit(EXIT_FAILURE);
}

/* Uniform distribution over [0,1]. */
static double random(struct ent_context * context)
{
        return rand() / (double)RAND_MAX;
}

int main(int nargc, char * argv[])
{
        /* Parse the input arguments. */
        double energy = (nargc > 1) ? atof(argv[1]) : 1E+12;
        enum ent_pid projectile = (nargc > 2) ? atoi(argv[2]) : ENT_PID_NU_TAU;
        enum ent_process process =
            (nargc > 3) ? atoi(argv[3]) : ENT_PROCESS_DIS_CC;
        enum ent_pid target = (nargc > 4) ? atoi(argv[4]) : ENT_PID_PROTON;
        int events = (nargc > 5) ? atoi(argv[5]) : 100000;

        /* Register the error handler for ENT library functions. */
        ent_error_handler_set(&handle_error);

        /* Create a new Physics environment. */
        ent_physics_create(&physics, "data/pdf/CT14nnlo_0000.dat");

        /* Compute the DCS by numeric integration. */
        const double Z = (target == ENT_PID_PROTON) ? 1. : 0.;
        const double A = 1.;
        const double ymin = 1E-07;
        const double xmin = 1E-12;
        const int ny = 141;
        const int nx = 241;
        const double dly = -log(ymin) / (ny - 1);
        const double dlx = -log(xmin) / (nx - 1);
        int i;
        FILE * stream = fopen("dcs.dat", "w+");
        for (i = 0; i < ny; i++) {
                double y = (i < ny - 1) ? ymin * exp(i * dly) : 1. - 1E-08;
                double dcs = 0.;
                int j;
                double x0 = xmin;
                double d0;
                ent_physics_dcs(
                    physics, projectile, energy, Z, A, process, x0, y, &d0);
                for (j = 1; j < nx; j++) {
                        const double x1 = xmin * exp(j * dlx);
                        double d1;
                        ent_physics_dcs(physics, projectile, energy, Z, A,
                            process, x1, y, &d1);
                        dcs += 0.5 * (x1 - x0) * (d1 + d0);
                        x0 = x1;
                        d0 = d1;
                }
                fprintf(stream, "%12.5E %12.5E\n", y, dcs);
        }
        fclose(stream);

        /* Instanciate a new simulation context. */
        struct ent_context context = { NULL, &random, 1 };

        /* Run a batch of Monte-Carlo vertices. */
        stream = fopen("y.dat", "w+");
        for (i = 0; i < events; i++) {
                struct ent_state neutrino = { projectile, energy, 0., 0., 1.,
                        { 0., 0., 0. }, { 0., 0., 1. } };
                struct ent_state product;
                ent_vertex(
                    physics, &context, process, target, &neutrino, &product);
                if (process == ENT_PROCESS_DIS_CC)
                        fprintf(stream, "%12.5lE\n", neutrino.energy / energy);
                else
                        fprintf(stream, "%12.5lE\n", product.energy / energy);
        }
        fclose(stream);

        /* Finalise and exit to the OS. */
        ent_physics_destroy(&physics);
        exit(EXIT_SUCCESS);
}
