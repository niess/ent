/* Standard library includes. */
#include <stdio.h>
#include <stdlib.h>

/* The ENT API. */
#include "ent.h"

/* Handle the Physics. */
static struct ent_physics * physics = NULL;

/* Error handler: dump any error message and exit to the OS. */
void handle_error(enum ent_return rc, ent_function_t * caller)
{
        /* Dump an error message. */
        ent_error_print(stderr, rc, caller, "\t", "\n");
        fprintf(stderr, "\n");

        /* Finalise and exit to the OS. */
        ent_physics_destroy(&physics);
        exit(EXIT_FAILURE);
}

int main()
{
        /* Register the error handler for GULL library functions. */
        ent_error_handler_set(&handle_error);

        /* Create a new Physics environment. */
        ent_physics_create(&physics, "data/pdf/cteq6d.tbl");

        /* Test some DCS. */
        const double energy = 6.3E+06;
        const double x = 1E-05;
        const double y = 1E-02;
        double dcs;
        ent_physics_dcs(physics, ENT_PROJECTILE_NU_E_BAR, energy, 0.5, 1.,
            ENT_PROCESS_GLASHOW_HADRON, x, y, &dcs);
        printf("DCS(%.5lE, %.5lE, %.5lE) = %.5lE\n", energy, x, y, dcs);

        const double q2 = energy * x * y * 0.931;
        double pdf;
        ent_physics_pdf(physics, ENT_PARTON_U, 0.5, q2, &pdf);
        printf("PDF(%.5lE, %.5lE) = %.5lE\n", x, q2, pdf);

        /* Finalise and exit to the OS. */
        ent_physics_destroy(&physics);
        exit(EXIT_SUCCESS);
}
