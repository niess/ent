/*
 * Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC
 * Author: Valentin NIESS (niess@in2p3.fr)
 *
 * an Engine for high energy Neutrinos Transport (ENT)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */

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
        double energy = (nargc > 1) ? atof(argv[1]) : 1E+09;
        enum ent_pid projectile = (nargc > 2) ? atoi(argv[2]) : ENT_PID_NU_TAU;
        enum ent_process process =
            (nargc > 3) ? atoi(argv[3]) : ENT_PROCESS_DIS_CC;
        double Z = (nargc > 4) ? atof(argv[4]) : 0.5;
        double A = (nargc > 5) ? atof(argv[5]) : 1.;
        int events = (nargc > 6) ? atoi(argv[6]) : 1000000;

        /* Register the error handler for ENT library functions. */
        ent_error_handler_set(&handle_error);

        /* Create a new Physics environment. */
        ent_physics_create(&physics, "data/pdf/CT14nnlo_0000.dat");

        /* Compute the DCS by numeric integration. */
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

        /* Instanciate a new simulation context and target medium. */
        struct ent_context context = { NULL, &random, 1 };
        struct ent_medium medium = { Z, A };

        /* Run a batch of Monte-Carlo vertices. */
        stream = fopen("y.dat", "w+");
        for (i = 0; i < events; i++) {
                struct ent_state neutrino = { projectile, energy, 0., 0., 1.,
                        { 0., 0., 0. }, { 0., 0., 1. } };
                struct ent_state product;
                ent_vertex(
                    physics, &context, &neutrino, &medium, process, &product);
                if (process == ENT_PROCESS_DIS_CC)
                        fprintf(stream, "%12.5lE\n", product.energy / energy);
                else
                        fprintf(stream, "%12.5lE\n", neutrino.energy / energy);
        }
        fclose(stream);

        /* Finalise and exit to the OS. */
        ent_physics_destroy(&physics);
        exit(EXIT_SUCCESS);
}
