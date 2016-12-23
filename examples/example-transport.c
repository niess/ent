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
static double medium(struct ent_state * state, struct ent_medium ** medium)
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

int main()
{
        /* Register the error handler for GULL library functions. */
        ent_error_handler_set(&handle_error);

        /* Create a new Physics environment. */
        ent_physics_create(&physics, "data/pdf/CT14nnlo_0000.dat");

        /* Instanciate a new simulation context. */
        struct ent_context context = { &medium, &random, 1, 0., 6400E+03, 0. };

        /* Run a Monte-Carlo transport. */
        struct ent_state neutrino = { ENT_PID_NU_TAU, 1E+09, 0., 0., 1.,
                { 0., 0., 0. }, { 0., 0., 1. } };
        struct ent_state product;
        enum ent_event event;
        ent_transport(physics, &context, &neutrino, &product, &event);

        /* Finalise and exit to the OS. */
        ent_physics_destroy(&physics);
        exit(EXIT_SUCCESS);
}
