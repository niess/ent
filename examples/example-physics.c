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

/* The ENT library. */
#include "ent.h"


int main(int nargc, char * argv[])
{
        /* Parse the input arguments. */
        enum ent_pid projectile = (nargc > 1) ? atoi(argv[1]) : ENT_PID_NU_TAU;
        enum ent_process process =
            (nargc > 2) ? atoi(argv[2]) : ENT_PROCESS_DIS_CC;
        double energy = (nargc > 3) ? atof(argv[4]) : 1E+09; /* GeV */
        double Z = (nargc > 4) ? atof(argv[4]) : 0.5;
        double A = (nargc > 5) ? atof(argv[5]) : 1.;

        /* Create the physics using DIS SFs data computed with APFEL.
         *
         * Alternatively, a PDF file can be provided as well, in LHA format.
         * However, in this case ENT uses LO expressions for DIS structure
         * functions.
         */
        struct ent_physics * physics;
        ent_physics_create(&physics, "share/ent/BGR18-sf.ent");

        /* Print physics meta-data. */
        puts(ent_physics_metadata(physics));

        /* Get the total cross-section. */
        double cs0;
        ent_physics_cross_section(
            physics, projectile, energy, Z, A, process, &cs0);
        printf("sigma0(%G) = %.3E\n", energy, cs0);

        /* Rescale the physics and print the new cross-section. */
        ent_physics_rescale(physics, "share/ent/BGR18-cross-section.txt");

        double cs1;
        ent_physics_cross_section(
            physics, projectile, energy, Z, A, process, &cs1);
        printf("sigma1(%G) = %.3E (%+.1f%%)\n", energy, cs1,
            100. * (cs1 / cs0 - 1));

        /* N.B.: Physics data can be dumped to a binary file (using
         * `ent_physics_dump`) and directly re-loaded, for faster
         * initialisation..
         */

        /* Finalise and exit to the OS. */
        ent_physics_destroy(&physics);
        exit(EXIT_SUCCESS);
}
