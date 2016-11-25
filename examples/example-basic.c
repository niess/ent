/* Standard library includes. */
#include <stdio.h>
#include <stdlib.h>

/* The ENT API. */
#include "ent.h"

int main()
{
        struct ent_dcs * dcs;
        ent_dcs_create("data/pdf/cteq6d.tbl", &dcs);

        const double energy = 1E+06;
        const double x = 1E-05;
        const double y = 1E-02;
        double value;
        ent_dcs_compute(dcs, ENT_PROJECTILE_NU_TAU, energy, 0.5, 1.,
            ENT_PROCESS_INVERSE_TAU, x, y, &value);
        printf("DCS(%.5lE, %.5lE, %.5lE) = %.5lE\n", energy, x, y, value);

        ent_dcs_destroy(&dcs);
        exit(EXIT_SUCCESS);
}
