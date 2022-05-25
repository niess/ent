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
        double energy = (nargc > 1) ? atof(argv[1]) : 1E+09; /* GeV */
        enum ent_pid projectile = (nargc > 2) ? atoi(argv[2]) : ENT_PID_NU_TAU;
        enum ent_process process =
            (nargc > 3) ? atoi(argv[3]) : ENT_PROCESS_DIS_CC;
        double Z = (nargc > 4) ? atof(argv[4]) : 0.5;
        double A = (nargc > 5) ? atof(argv[5]) : 1.;

        /* Create the physics using a data dump. */
        struct ent_physics * physics;
        ent_physics_create(&physics, "share/ent/CMS11-physics.ent");

        /* Create a new simulation context. */
        struct ent_context * context;
        ent_context_create(&context);

        /* Initialise the projectile Monte Carlo state. */
        struct ent_state lepton = {
            .pid = projectile,
            .energy = energy,
            .weight = 1.,
            .direction = { 0., 0., 1. } /* Must be a unit vector! */
        };

        /* Instanciate the target medium. */
        struct ent_medium medium = { Z, A };

        /* Generate a Monte-Carlo collision. */
        struct ent_state products[ENT_PRODUCTS_SIZE];
        ent_collide(physics, context, &lepton, &medium, process, products);

        /* Print collision products. */
        puts(" PID     energy     theta");
        puts("          (GeV)     (rad)");

        int i;
        for (i = 0; i < ENT_PRODUCTS_SIZE + 1; i++) {
                struct ent_state * product;
                if (i == 0) {
                        /* After the collision, the initial ``projectile'' state
                         * contains the lepton product.
                         .*/
                        product = &lepton;
                } else {
                        /* Other collision products are provided separately.
                         * Note that a `NULL` pointer can be provided, if this
                         * is not of interest.
                         */
                        product = products + (i - 1);
                }

                if (product->weight <= 0.) continue;

                const double theta = acos(product->direction[2]);
                printf("%4d   %.3E  %.3E\n",
                    product->pid, product->energy, theta);

                /* N.B.: ENT does not hadronize hadronic products. However,
                 * semi-leptonic decays of top quarks are simulated. A pid of
                 * `100` is a special value used by ENT in order to designate a
                 * generic hadronic product, i.e. not hadronized.
                 */
        }

        /* Finalise and exit to the OS. */
        ent_context_destroy(&context);
        ent_physics_destroy(&physics);

        exit(EXIT_SUCCESS);
}
