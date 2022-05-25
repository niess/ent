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

/* The ENT library. */
#include "ent.h"


/* Density callback with a uniform density. */
static double density(
    struct ent_medium * medium, struct ent_state * state, double * density)
{
        *density = 2.65E+03; /* kg / m^3. */
        return 0.;
}

/* Medium callback with a single infinite medium. */
static double medium(struct ent_context * context, struct ent_state * state,
    struct ent_medium ** medium)
{
        /* Standard rock, see e.g.
         * https://pdg.lbl.gov/2021/AtomicNuclearProperties/HTML/standard_rock.html
         */
        static struct ent_medium medium_ = { 11., 22., &density };
        *medium = &medium_;
        return 0.;
}

int main(int nargc, char * argv[])
{
        /* Parse the input arguments. */
        double energy = (nargc > 1) ? atof(argv[1]) : 1E+06;
        enum ent_pid projectile =
            (nargc > 2) ? atoi(argv[2]) : ENT_PID_NU_E;
        double depth = (nargc > 3) ? atof(argv[3]) : 6400E+03;

        /* Create the physics using a data dump. */
        struct ent_physics * physics;
        ent_physics_create(&physics, "share/ent/BGR18-physics.ent");

        /* Create an configure the simulation context. */
        struct ent_context * context;
        ent_context_create(&context);

        context->medium = &medium;
        context->distance_max = depth;

        /* Initialise the neutrino state */
        struct ent_state lepton = {
            .pid = projectile,
            .energy = energy,
            .weight = 1.,
            .direction = { 0., 0., 1. } /* Must be a unit vector! */
        };

        /* Transport a Monte-Carlo event through a given depth of rock. */
        enum ent_event event;
        for (;;) {
                ent_transport(
                    physics, context, &lepton, NULL, &event);

                if ((lepton.energy <= 0.) || (lepton.pid != projectile)) {
                        puts("projectile did not exit.");
                        break;
                }

                if (event == ENT_EVENT_LIMIT_DISTANCE) {
                        printf("projectile exits with %G GeV left\n",
                            lepton.energy);
                        break;
                }
        }

        /* Finalise and exit to the OS. */
        ent_context_destroy(&context);
        ent_physics_destroy(&physics);

        exit(EXIT_SUCCESS);
}
