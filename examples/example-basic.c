/* Standard library includes. */
#include <stdio.h>
#include <stdlib.h>

/* The ANT API. */
#include "ant.h"

int main()
{
        struct ant_dcs * dcs;
        ant_dcs_create("data/pdf/cteq6d.tbl", &dcs);

        const double energy = 1E+06;
        const double x = 1E-06;
        const double y = 1E-02;
        double value;
        ant_dcs_compute(dcs, ANT_PROJECTILE_NU_TAU, energy, 0.5, 1.,
            ANT_PROCESS_CC, x, y, &value);
        printf("DCS(%.5lE, %.5lE, %.5lE) = %.5lE\n", energy, x, y, value);

        ant_dcs_destroy(&dcs);
        exit(EXIT_SUCCESS);
}
