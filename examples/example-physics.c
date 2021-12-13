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

int main(int nargc, char * argv[])
{
        /* Parse the input arguments. */
        enum ent_pid projectile = (nargc > 1) ? atoi(argv[1]) : ENT_PID_NU_TAU;
        enum ent_process process =
            (nargc > 2) ? atoi(argv[2]) : ENT_PROCESS_DIS_CC;
        double Z = (nargc > 3) ? atof(argv[3]) : 0.5;
        double A = (nargc > 4) ? atof(argv[4]) : 1.;

        /* Register the error handler for ENT library functions. */
        ent_error_handler_set(&handle_error);

        /* Create a new Physics environment. */
        ent_physics_create(
            &physics, "share/pdf/CT14nlo_0000.dat", "share/cs/CSMS.txt");

        /* Get the total cross-section. */
        const double Emin = 1E+00;
        const double Emax = 1E+12;
        const int nE = 241;
        const double lnE = log(Emax / Emin) / (nE - 1);
        int i;
        for (i = 0; i < nE; i++) {
                const double energy = Emin * exp(i * lnE);
                double cs;
                ent_physics_cross_section(
                    physics, projectile, energy, Z, A, process, &cs);
                printf("%.5lE %.5lE\n", energy, cs);
        }

        /* Finalise and exit to the OS. */
        ent_physics_destroy(&physics);
        exit(EXIT_SUCCESS);
}
