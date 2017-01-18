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

int main()
{
        /* Register the error handler for ENT library functions. */
        ent_error_handler_set(&handle_error);

        /* Create a new Physics environment. */
        ent_physics_create(&physics, "data/pdf/CT14nnlo_0000.dat");

        /* Compute the DCS. */
        double energy = 1E+05;
        const double ymin = 1E-07;
        const double xmin = 1E-12;
        const int ny = 141;
        const int nx = 241;
        const double dly = -log(ymin) / (ny - 1);
        const double dlx = -log(xmin) / (nx - 1);
        int i;
        FILE * stream = fopen("dcs.dat", "w+");
        for (i = 0; i < ny - 1; i++) {
                const double y = ymin * exp(i * dly);
                double dcs = 0.;
                int j;
                double x0 = xmin;
                double d0;
                ent_physics_dcs(physics, ENT_PID_NU_TAU, energy, 1., 1.,
                    ENT_PROCESS_DIS_CC, x0, y, &d0);
                for (j = 1; j < nx; j++) {
                        const double x1 = xmin * exp(j * dlx);
                        double d1;
                        ent_physics_dcs(physics, ENT_PID_NU_TAU, energy, 1., 1.,
                            ENT_PROCESS_DIS_NC, x1, y, &d1);
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
        for (i = 0; i < 100000; i++) {
                struct ent_state neutrino = { ENT_PID_NU_TAU, energy, 0., 0.,
                        1., { 0., 0., 0. }, { 0., 0., 1. } };
                struct ent_state product;
                ent_vertex(physics, &context, ENT_PROCESS_DIS_NC,
                    ENT_PID_PROTON, &neutrino, &product);
                fprintf(stream, "%12.5lE\n", product.energy / energy);
        }
        fclose(stream);

        /* Finalise and exit to the OS. */
        ent_physics_destroy(&physics);
        exit(EXIT_SUCCESS);
}
