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
        ent_physics_create(&physics, "data/pdf/CT14nnlo_0000.dat");

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
