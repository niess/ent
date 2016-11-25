/* Standard library includes. */
#include <stdio.h>
#include <stdlib.h>

/* The ENT API. */
#include "ent.h"

int main()
{
        struct ent_physics * physics;
        ent_physics_create(&physics, "data/pdf/cteq6d.tbl");

        const double energy = 1E+06;
        const double x = 1E-05;
        const double y = 1E-02;
        double value;
        ent_physics_dcs(physics, ENT_PROJECTILE_NU_TAU, energy, 0.5, 1.,
            ENT_PROCESS_DIS_CC, x, y, &value);
        printf("DCS(%.5lE, %.5lE, %.5lE) = %.5lE\n", energy, x, y, value);

        ent_physics_destroy(&physics);
        exit(EXIT_SUCCESS);
}
