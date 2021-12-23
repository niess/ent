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

/* Density callback with a uniform density. */
static double density(
    struct ent_medium * medium, struct ent_state * state, double * density)
{
        *density = 2.65E+03;
        return 0.;
}

/* Medium callback with a single infinite medium. */
static double medium(struct ent_context * context, struct ent_state * state,
    struct ent_medium ** medium)
{
        static struct ent_medium medium_ = { 13., 26., &density };
        *medium = &medium_;
        return 0.;
}

/* Uniform distribution over [0,1]. */
static double random(struct ent_context * context)
{
        return rand() / (double)RAND_MAX;
}

int main(int nargc, char * argv[])
{
        /* Parse the input arguments. */
        double energy = (nargc > 1) ? atof(argv[1]) : 5.7E+06;
        enum ent_pid projectile =
            (nargc > 2) ? atoi(argv[2]) : ENT_PID_NU_BAR_E;
        double depth = (nargc > 3) ? atof(argv[3]) : 6400E+03;
        int events = (nargc > 4) ? atoi(argv[4]) : 10000;

        /* Register the error handler for ENT library functions. */
        ent_error_handler_set(&handle_error);

        /* Create a new Physics environment. */
        ent_physics_create(&physics,
            "share/pdf/NNPDF31sx_nlo_as_0118_LHCb_nf_6_0000.dat",
            "share/cs/BGR18.txt");

        /* Instanciate a new simulation context. */
        struct ent_context context = { &medium, &random, NULL, NULL, depth };

        /* Run a batch of Monte-Carlo transports. */
        FILE * stream = fopen("transport.dat", "w+");
        int i;
        for (i = 0; i < events; i++) {
                struct ent_state neutrino = { projectile, energy, 0., 0., 1.,
                        { 0., 0., 0. }, { 0., 0., 1. } };
                enum ent_event event;
                for (;;) {
                        ent_transport(
                            physics, &context, &neutrino, NULL, &event);
                        if (neutrino.energy <= 0.) break;
                        if (event == ENT_EVENT_LIMIT_DISTANCE) {
                                fprintf(stream, "%3d %12.5lE\n", neutrino.pid,
                                    neutrino.energy);
                                break;
                        }
                        int aid = abs(neutrino.pid);
                        if ((aid != ENT_PID_NU_E) && (aid != ENT_PID_NU_MU) &&
                            (aid != ENT_PID_NU_TAU))
                                break;
                }
        }
        fclose(stream);

        /* Finalise and exit to the OS. */
        ent_physics_destroy(&physics);
        exit(EXIT_SUCCESS);
}
