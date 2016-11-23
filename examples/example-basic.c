/* Standard library includes. */
#include <stdlib.h>

/* The ANT API. */
#include "ant.h"

int main()
{
        struct ant_dcs * dcs;
        ant_dcs_create("data/pdf/cteq6d.tbl", &dcs);

        ant_dcs_destroy(&dcs);
        exit(EXIT_SUCCESS);
}
