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
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* The ENT API. */
#include "ent.h"

/* The electron mass, in GeV/c^2. */
#define ENT_MASS_ELECTRON 0.511E-03
/* The muon mass, in GeV/c^2. */
#define ENT_MASS_MUON 0.10566
/* The tau mass, in GeV/c^2. */
#define ENT_MASS_TAU 1.77682
/* The nucleon mass, in GeV/c^2. */
#define ENT_MASS_NUCLEON 0.931494
/* The W boson mass, in GeV/c^2. */
#define ENT_MASS_W 80.385
/* The Z boson mass, in GeV/c^2. */
#define ENT_MASS_Z 91.1876
/* The W boson width, in GeV/c^2. */
#define ENT_WIDTH_W 2.085
/* Branching ratio for W decay to electron. */
#define ENT_BR_W_TO_ELECTRON 0.1071
/* Branching ratio for W decay to muon. */
#define ENT_BR_W_TO_MUON 0.1063
/* Branching ratio for W decay to tau. */
#define ENT_BR_W_TO_TAU 0.1138
/* The muon decay length, in m. */
#define ENT_CTAU_MUON 659.09433
/* The tau decay length, in m. */
#define ENT_CTAU_TAU 8.709E-05
/* The Fermi coupling constant GF/(hbar*c)^3, in GeV^-2. */
#define ENT_PHYS_GF 1.1663787E-05
/* The Planck constant as hbar*c, in GeV * m. */
#define ENT_PHYS_HBC 1.97326978E-16
/* The Weinberg angle at MZ, as sin(theta_W)^2. */
#define ENT_PHYS_SIN_THETA_W_2 0.231295
/* Avogadro's number, in mol^-1. */
#define ENT_PHYS_NA 6.022E+23

/* Energy range for tabulations. */
#define ENERGY_MIN 1E+02
#define ENERGY_MAX 1E+12
#define ENERGY_N 301

/* Sampling for the tabulations of DIS cumulative cross-sections. */
#define DIS_Y_N 301

/* Exponent for DIS bias model, w.r.t. x. */
#define DIS_X_EXPONENT 2.5

#ifndef M_PI
/* Define pi, if unknown. */
#define M_PI 3.14159265358979323846
#endif

/* ENT data format tag. */
#define ENT_FORMAT_TAG "/ent/"

/* ENT data format version. */
#define ENT_FORMAT_VERSION 3

/* Indices for Physics processes with a specific projectile and target. */
enum proget_index {
        /* Backward decay from a tau. */
        PROGET_BACKWARD_DECAY_TAU = -2,
        /* Backward decay from a muon. */
        PROGET_BACKWARD_DECAY_MUON = -1,
        /* Charged current DIS of a neutrino on a neutron. */
        PROGET_CC_OTHER_NU_NEUTRON = 0,
        /* Neutral current DIS of a neutrino on a neutron. */
        PROGET_NC_NU_NEUTRON,
        /* Charged current DIS of an anti-neutrino on a neutron. */
        PROGET_CC_OTHER_NU_BAR_NEUTRON,
        /* Neutral current DIS of an anti-neutrino on a neutron. */
        PROGET_NC_NU_BAR_NEUTRON,
        /* Charged current DIS of a neutrino on a proton. */
        PROGET_CC_OTHER_NU_PROTON,
        /* Neutral current DIS of a neutrino on a proton. */
        PROGET_NC_NU_PROTON,
        /* Charged current DIS of an anti-neutrino on a proton. */
        PROGET_CC_OTHER_NU_BAR_PROTON,
        /* Neutral current DIS of an anti-neutrino on a proton. */
        PROGET_NC_NU_BAR_PROTON,
        /* Top production for a CC neutrino on a neutron. */
        PROGET_CC_TOP_NU_NEUTRON,
        /* Top production for a CC anti-neutrino on a neutron. */
        PROGET_CC_TOP_NU_BAR_NEUTRON,
        /* Top production for a CC neutrino on a proton. */
        PROGET_CC_TOP_NU_PROTON,
        /* Top production for a CC anti-neutrino on a proton. */
        PROGET_CC_TOP_NU_BAR_PROTON,
        /* Number of DIS PRoces-projectile-tarGET cases. */
        PROGET_N_DIS,
        /* Elastic scattering of a nu_e on a an electron. */
        PROGET_ELASTIC_NU_E = PROGET_N_DIS,
        /* Elastic scattering of an anti nu_e on a an electron. */
        PROGET_ELASTIC_NU_BAR_E,
        /* Elastic scattering of a nu_mu on a an electron. */
        PROGET_ELASTIC_NU_MU,
        /* Elastic scattering of an anti nu_mu on a an electron. */
        PROGET_ELASTIC_NU_BAR_MU,
        /* Elastic scattering of a nu_tau on a an electron. */
        PROGET_ELASTIC_NU_TAU,
        /* Elastic scattering of an anti nu_tau on a an electron. */
        PROGET_ELASTIC_NU_BAR_TAU,
        /* Inverse muon decay with a nu_mu projectile. */
        PROGET_INVERSE_NU_MU_MU,
        /* Inverse tau decay with a nu_tau projectile. */
        PROGET_INVERSE_NU_TAU_TAU,
        /* Inverse muon decay with an anti nu_e projectile. */
        PROGET_INVERSE_NU_BAR_E_MU,
        /* Inverse tau decay with an anti nu_e projectile. */
        PROGET_INVERSE_NU_BAR_E_TAU,
        /* Anti nu_e projectile on an electron with hadrons production. */
        PROGET_GLASHOW_HADRONS,
        /* Number of PROcess-projectile-tarGET cases. */
        PROGET_N
};

/* Compute the process-target index. */
static enum ent_return proget_compute(enum ent_pid projectile,
    enum ent_pid target, enum ent_process process, int * proget)
{
        *proget = -1;

        const int apid = abs(projectile);
        if ((apid != ENT_PID_NU_E) && (apid != ENT_PID_NU_MU) &&
            (apid != ENT_PID_NU_TAU))
                return ENT_RETURN_DOMAIN_ERROR;
        if (process == ENT_PROCESS_DIS_CC) {
                return ENT_RETURN_DOMAIN_ERROR;
        } else if (process == ENT_PROCESS_DIS_CC_OTHER) {
                if (target == ENT_PID_NEUTRON)
                        *proget = (projectile > 0) ?
                            PROGET_CC_OTHER_NU_NEUTRON :
                            PROGET_CC_OTHER_NU_BAR_NEUTRON;
                else if (target == ENT_PID_PROTON)
                        *proget = (projectile > 0) ?
                            PROGET_CC_OTHER_NU_PROTON :
                            PROGET_CC_OTHER_NU_BAR_PROTON;
                else
                        return ENT_RETURN_DOMAIN_ERROR;
        } else if (process == ENT_PROCESS_DIS_CC_TOP) {
                if (target == ENT_PID_NEUTRON)
                        *proget = (projectile > 0) ?
                            PROGET_CC_TOP_NU_NEUTRON :
                            PROGET_CC_TOP_NU_BAR_NEUTRON;
                else if (target == ENT_PID_PROTON)
                        *proget = (projectile > 0) ?
                            PROGET_CC_TOP_NU_PROTON :
                            PROGET_CC_TOP_NU_BAR_PROTON;
                else
                        return ENT_RETURN_DOMAIN_ERROR;
        } else if (process == ENT_PROCESS_DIS_NC) {
                if (target == ENT_PID_NEUTRON)
                        *proget = (projectile > 0) ? PROGET_NC_NU_NEUTRON :
                                                     PROGET_NC_NU_BAR_NEUTRON;
                else if (target == ENT_PID_PROTON)
                        *proget = (projectile > 0) ? PROGET_NC_NU_PROTON :
                                                     PROGET_NC_NU_BAR_PROTON;
                else
                        return ENT_RETURN_DOMAIN_ERROR;
        } else {
                if (target != ENT_PID_ELECTRON) return ENT_RETURN_DOMAIN_ERROR;
                if (process == ENT_PROCESS_ELASTIC) {
                        if (projectile == ENT_PID_NU_E)
                                *proget = PROGET_ELASTIC_NU_E;
                        else if (projectile == ENT_PID_NU_BAR_E)
                                *proget = PROGET_ELASTIC_NU_BAR_E;
                        else if (projectile == ENT_PID_NU_MU)
                                *proget = PROGET_ELASTIC_NU_MU;
                        else if (projectile == ENT_PID_NU_BAR_MU)
                                *proget = PROGET_ELASTIC_NU_BAR_MU;
                        else if (projectile == ENT_PID_NU_TAU)
                                *proget = PROGET_ELASTIC_NU_TAU;
                        else
                                *proget = PROGET_ELASTIC_NU_BAR_TAU;
                } else if (process == ENT_PROCESS_INVERSE_MUON) {
                        if (projectile == ENT_PID_NU_MU)
                                *proget = PROGET_INVERSE_NU_MU_MU;
                        else if (projectile == ENT_PID_NU_BAR_E)
                                *proget = PROGET_INVERSE_NU_BAR_E_MU;
                        else
                                return ENT_RETURN_DOMAIN_ERROR;
                } else if (process == ENT_PROCESS_INVERSE_TAU) {
                        if (projectile == ENT_PID_NU_TAU)
                                *proget = PROGET_INVERSE_NU_TAU_TAU;
                        else if (projectile == ENT_PID_NU_BAR_E)
                                *proget = PROGET_INVERSE_NU_BAR_E_TAU;
                        else
                                return ENT_RETURN_DOMAIN_ERROR;
                } else if (process == ENT_PROCESS_GLASHOW_HADRON) {
                        if (projectile == ENT_PID_NU_BAR_E)
                                *proget = PROGET_GLASHOW_HADRONS;
                        else
                                return ENT_RETURN_DOMAIN_ERROR;
                } else
                        return ENT_RETURN_DOMAIN_ERROR;
        }
        return ENT_RETURN_SUCCESS;
}

/* Opaque Physics object. */
struct ent_physics {
        /* Index to the PDF data. */
        struct grid * pdf;
        /* Index to the DIS SFs data. */
        struct grid * sf[PROGET_N_DIS];
        /* Effective thresholds for top production. */
        double mteff[5];
        /* Entry point for kinematic cross-sections, computed from DCSs. */
        double * cs_k;
        /* Entry point for transport cross-sections, which the user might
         * override.
         */
        double * cs_t;
        /* Entry point for the probability density function for DIS. */
        double * dis_pdf;
        /* Entry point for the cumulative density function for DIS. */
        double * dis_cdf;
        /* Entry point for DIS DDCS support in x. */
        double * dis_xlim;
        /* Entry point for forward DIS DCS support in y, as function of Ei. */
        double * dis_ylim_f;
        /* Entry point for backward DIS DCS support in y, as function of Ef. */
        double * dis_ylim_b;
        /* Entry point for backward DIS DCS support in y, as function of Eh. */
        double * dis_ylim_h;
        /* Bounds for x and Q2, in DIS. */
        double dis_Q2min, dis_Q2max;
        /* Placeholder for dynamic data. */
        char data[];
};

/* Compute the padded memory size.
 *
 * The padded memory size is the smallest integer multiple of *pad_size* and
 * greater or equal to *size*. It allows to align memory addresses on multiples
 * of pad_size.
 */
static int memory_padded_size(int size, int pad_size)
{
        int i = size / pad_size;
        if ((size % pad_size) != 0) i++;
        return i * pad_size;
}

/* Generic allocator for a Physics object. */
static void * physics_create(void)
{
        const int np = PROGET_N - 1;
        void * v = malloc(sizeof(struct ent_physics) +
            (2 * np * ENERGY_N +
             5 * PROGET_N_DIS * ENERGY_N * DIS_Y_N +
             6 * PROGET_N_DIS * ENERGY_N) * sizeof(double));
        if (v == NULL) return NULL;
        struct ent_physics * p = v;
        p->pdf = NULL;
        memset(p->sf, 0x0, sizeof(p->sf));
        p->cs_k = (double *)(p->data);
        p->cs_t = p->cs_k + np * ENERGY_N;
        p->dis_pdf = p->cs_t + np * ENERGY_N;
        p->dis_cdf = p->dis_pdf + PROGET_N_DIS * ENERGY_N * DIS_Y_N;
        p->dis_xlim = p->dis_cdf + PROGET_N_DIS * ENERGY_N * DIS_Y_N;
        p->dis_ylim_f = p->dis_xlim +
            3 * PROGET_N_DIS * DIS_Y_N * ENERGY_N;
        p->dis_ylim_b = p->dis_ylim_f + 2 * PROGET_N_DIS * ENERGY_N;
        p->dis_ylim_h = p->dis_ylim_b + 2 * PROGET_N_DIS * ENERGY_N;

        return v;
}

/* Get a return code as a string. */
const char * ent_error_string(enum ent_return code)
{
        static const char * msg[ENT_N_RETURNS] = { "Operation succeeded",
                "Bad address", "Value is out of validity range",
                "Invalid file format", "Couldn't read file",
                "Not enough memory", "No such file or directory" };

        if ((code < 0) || (code >= ENT_N_RETURNS))
                return NULL;
        else
                return msg[code];
}

/* Get a library function name as a string. */
const char * ent_error_function(ent_function_t * caller)
{
#define REGISTER_FUNCTION(function)                                            \
        if (caller == (ent_function_t *)function) return #function

        /* API functions with error codes. */
        REGISTER_FUNCTION(ent_physics_create);
        REGISTER_FUNCTION(ent_physics_cross_section);
        REGISTER_FUNCTION(ent_physics_dcs);
        REGISTER_FUNCTION(ent_physics_sf);
        REGISTER_FUNCTION(ent_physics_pdf);
        REGISTER_FUNCTION(ent_collide);
        REGISTER_FUNCTION(ent_transport);

        /* Other API functions. */
        REGISTER_FUNCTION(ent_physics_destroy);
        REGISTER_FUNCTION(ent_error_string);
        REGISTER_FUNCTION(ent_error_print);
        REGISTER_FUNCTION(ent_error_function);
        REGISTER_FUNCTION(ent_error_handler_get);
        REGISTER_FUNCTION(ent_error_handler_set);

        return NULL;
#undef REGISTER_FUNCTION
}

/* The user supplied error handler, if any. */
static ent_handler_cb * _handler = NULL;

/* Getter for the error handler. */
ent_handler_cb * ent_error_handler_get(void) { return _handler; }

/* Setter for the error handler. */
void ent_error_handler_set(ent_handler_cb * handler) { _handler = handler; }

/* Encapsulation of `return` with an error handler. */
static enum ent_return error_handle(
    enum ent_return code, ent_function_t * caller)
{
        if (code == ENT_RETURN_SUCCESS)
                return code;
        else if (_handler != NULL)
                _handler(code, caller);
        return code;
}

/* Helper macros for returning an encapsulated error code. */
#define ENT_ACKNOWLEDGE(caller)                                                \
        ent_function_t * _caller = (ent_function_t *)caller

#define ENT_RETURN(code) return error_handle(code, _caller)

/* Format an error message in JSON. */
void ent_error_print(FILE * stream, enum ent_return code,
    ent_function_t * function, const char * tabulation, const char * newline)
{
        const char * tab = (tabulation == NULL) ? "" : tabulation;
        const char * cr = (newline == NULL) ? "" : newline;

        fprintf(stream, "{%s%s\"code\" : %d,%s%s\"message\" : \"%s\"", cr, tab,
            code, cr, tab, ent_error_string(code));
        if (function != NULL)
                fprintf(stream, ",%s%s\"function\" : \"%s\"", cr, tab,
                    ent_error_function(function));
        fprintf(stream, "%s}", cr);
}

/* Buffer for reading data from a text file. */
struct file_buffer {
        jmp_buf env;
        FILE * stream;
        enum ent_return code;
        unsigned int size;
        char * cursor;
        char line[];
};

/* Memory blocks size for file buffers. */
#define FILE_BLOCK_SIZE 2048

/* Create a new file buffer */
static enum ent_return file_buffer_create(
    struct file_buffer ** buffer, FILE * stream, jmp_buf env)
{
        if ((*buffer = malloc(FILE_BLOCK_SIZE)) == NULL)
                return ENT_RETURN_MEMORY_ERROR;
        memcpy((*buffer)->env, env, sizeof(jmp_buf));
        (*buffer)->stream = stream;
        (*buffer)->code = ENT_RETURN_SUCCESS;
        (*buffer)->size = FILE_BLOCK_SIZE - sizeof(**buffer);
        (*buffer)->cursor = NULL;

        return ENT_RETURN_SUCCESS;
}

/* Raise a file parsing error. */
static void file_raise_error(struct file_buffer * buffer, enum ent_return code)
{
        buffer->code = code;
        longjmp(buffer->env, 1);
}

/* Load next line from file to buffer. */
static enum ent_return file_get_line_(struct file_buffer ** buffer, int skip)
{
        /* Get a new line. */
        char * ptr = (*buffer)->line;
        int size = (*buffer)->size;
        for (; skip >= 0; skip--) {
                for (;;) {
                        const char check = 0x1;
                        ptr[size - 1] = check;
                        if (fgets(ptr, size, (*buffer)->stream) == NULL)
                                return ENT_RETURN_IO_ERROR;

                        /* Check that no overflow occured. */
                        if (ptr[size - 1] == check) break;

                        /* Get more memory then ... */
                        void * tmp = realloc(*buffer,
                            sizeof(**buffer) + (*buffer)->size +
                                FILE_BLOCK_SIZE);
                        if (tmp == NULL) return ENT_RETURN_MEMORY_ERROR;
                        *buffer = tmp;
                        ptr = (*buffer)->line + (*buffer)->size - 1;
                        (*buffer)->size += FILE_BLOCK_SIZE;
                        size = FILE_BLOCK_SIZE + 1;
                }
        }

        /* Reset the line cursor and return. */
        (*buffer)->cursor = (*buffer)->line;
        return ENT_RETURN_SUCCESS;
}

#undef FILE_BLOCK_SIZE /* No more needed. */

static void file_get_line(struct file_buffer ** buffer, int skip)
{
        enum ent_return rc;
        if ((rc = file_get_line_(buffer, skip)) != ENT_RETURN_SUCCESS)
                file_raise_error(*buffer, rc);
}

/* Count the number of words in the file's buffer. */
static int file_count_words(struct file_buffer * buffer)
{
        /* Scan the current buffer for a valid word. */
        int count = 0;
        char * p = buffer->cursor;
        while (p != 0x0) {
                for (; isspace(*p); p++)
                        ;
                if (*p == 0x0) break;
                count++;
                for (; (*p != 0x0) && !isspace(*p); p++)
                        ;
        }
        return count;
}

/* Parse the next float in the buffer. */
static enum ent_return file_get_float_(struct file_buffer * buffer, float * f)
{
        char * endptr;
        *f = strtof(buffer->cursor, &endptr);
        if (buffer->cursor == endptr) return ENT_RETURN_FORMAT_ERROR;
        buffer->cursor = endptr;
        return ENT_RETURN_SUCCESS;
}

/* Parse a table of floats from an opened file using a buffer. */
static void file_get_table(
    struct file_buffer ** buffer, int size, float * table)
{
        enum ent_return rc;
        int failed = 0, n = 0;
        while (n < size) {
                if ((rc = file_get_float_(*buffer, table)) !=
                    ENT_RETURN_SUCCESS) {
                        if (failed) file_raise_error(*buffer, rc);
                        if ((rc = file_get_line_(buffer, 0)) !=
                            ENT_RETURN_SUCCESS)
                                file_raise_error(*buffer, rc);
                        failed = 1;
                } else {
                        failed = 0;
                        n++;
                        table++;
                }
        }
}

/* Containers for (x, Q2) grid data. */
struct grid {
        /* Link to the next table. */
        struct grid * next;
        /* The table size. */
        int nx, nQ2, nf;
        /* Links to the data tables. */
        float *x, *Q2, *lambda, *f;
        /* Placeholder for variable size data. */
        float data[];
};

/* Create a new grid instance. */
static struct grid * grid_create(int nx, int nQ2, int nf)
{
        const unsigned int size_x =
            memory_padded_size(nx * sizeof(float), sizeof(double));
        const unsigned int size_Q2 =
            memory_padded_size(nQ2 * sizeof(float), sizeof(double));
        const unsigned int size_lambda =
            memory_padded_size(nf * nQ2 * sizeof(float), sizeof(double));
        const unsigned int size_sf =
            memory_padded_size(nf * nx * nQ2 * sizeof(float), sizeof(double));
        unsigned int size =
            sizeof(struct grid) + size_x + size_Q2 + size_lambda + size_sf;

        struct grid * g;
        if ((g = malloc(size)) == NULL) {
                return NULL;
        }

        g->next = NULL;
        g->nx = nx;
        g->nQ2 = nQ2;
        g->nf = nf;
        g->x = g->data;
        g->Q2 = (float *)((char *)g->data + size_x);
        g->lambda = (float *)((char *)g->Q2 + size_Q2);
        g->f = (float *)((char *)g->lambda + size_lambda);

        return g;
}

/* Recursively destroy grid data. */
static void grid_destroy(struct grid ** grid)
{
        if (*grid == NULL) return;

        struct grid * g = *grid;
        while (g != NULL) {
                struct grid * next = g->next;
                free(g);
                g = next;
        }
        *grid = NULL;
}

/* Compute grid small x extrapolation. */
static enum ent_return grid_compute_lambda(struct grid * g, int * freeze)
{
        const unsigned int size_lambda =
            memory_padded_size(g->nf * g->nQ2 * sizeof(float), sizeof(double));
        memset(g->lambda, 0x0, size_lambda);

        if (g->x[0] <= 0.) {
                return ENT_RETURN_SUCCESS;
        }

        float lxi = log(g->x[1] / g->x[0]);
        if (lxi <= 0.) {
                return ENT_RETURN_FORMAT_ERROR;
        }
        lxi = 1. / lxi;

        int i;
        float * table;
        for (i = 0, table = g->lambda; i < g->nQ2; i++, table += g->nf) {
                const float * const f0 = g->f + i * g->nf;
                const float * const f1 = g->f + (g->nQ2 + i) * g->nf;
                int j;
                for (j = 0; j < g->nf; j++) {
                        if ((freeze != NULL) && (freeze[j])) {
                                continue;
                        } else {
                                const float f0j = f0[j];
                                const float f1j = f1[j];
                                if (f0j * f1j > 0.)
                                        table[j] = log(f0j / f1j) * lxi;
                        }
                }
        }

        return ENT_RETURN_SUCCESS;
}

/* Load a single grid in ENT format. */
static enum ent_return grid_load(
    struct grid ** g_ptr, FILE * stream, int * freeze)
{
        struct grid * g = NULL;
        enum ent_return rc;

        /* Read binary header. */
        int nx, nQ2, nf;
        if ((fread(&nx, sizeof(nx), 1, stream) != 1) ||
            (fread(&nQ2, sizeof(nQ2), 1, stream) != 1) ||
            (fread(&nf, sizeof(nf), 1, stream) != 1))
        {
                rc = ENT_RETURN_FORMAT_ERROR;
                goto error;
        }

        /* Allocate memory for tables and map it. */
        if ((g = grid_create(nx, nQ2, nf)) == NULL) {
                rc = ENT_RETURN_MEMORY_ERROR;
                goto error;
        }

        /* Load the grid data. */
        if (fread(g->x, sizeof(*g->x), nx, stream) != nx) {
                rc = ENT_RETURN_FORMAT_ERROR;
                goto error;
        }

        if (fread(g->Q2, sizeof(*g->Q2), nQ2, stream) != nQ2) {
                rc = ENT_RETURN_FORMAT_ERROR;
                goto error;
        }

        const int nr = nf * nx * nQ2;
        if (fread(g->f, sizeof(*g->f), nr, stream) != nr) {
                rc = ENT_RETURN_FORMAT_ERROR;
                goto error;
        }

        /* Compute lambda exponents. */
        if ((rc = grid_compute_lambda(g, freeze)) != ENT_RETURN_SUCCESS) {
                goto error;
        }

        /* Register the table and return. */
        *g_ptr = g;

        return ENT_RETURN_SUCCESS;
error:
        free(g);
        *g_ptr = NULL;

        return rc;
}

/* Recursive bracketing of a table value, using a dichotomy. */
static void table_bracket(const float * table, float value, int * p1, int * p2)
{
        int i3 = (*p1 + *p2) / 2;
        if (value >= table[i3])
                *p1 = i3;
        else
                *p2 = i3;
        if (*p2 - *p1 >= 2) table_bracket(table, value, p1, p2);
}

static void grid_compute(
    const struct grid * g, float x, float Q2, float * f)
{
        memset(f, 0x0, g->nf * sizeof(*f));

        /* Check the bounds. */
        while ((Q2 < g->Q2[0]) || (Q2 >= g->Q2[g->nQ2 - 1]) ||
            (x >= g->x[g->nx - 1])) {
                g = g->next;
                if (g == NULL) return;
        }

        if (x < g->x[0]) {
                /* Extrapolate with a power law for small x values, i.e.
                 * f(x) ~ 1 / x**lambda(Q2).
                 */
                int iQ0, iQ1;
                float hQ;
                if (Q2 >= g->Q2[g->nQ2 - 1]) {
                        iQ0 = iQ1 = g->nQ2 - 1;
                        hQ = 1.;
                } else {
                        iQ0 = 0;
                        iQ1 = g->nQ2 - 1;
                        table_bracket(g->Q2, Q2, &iQ0, &iQ1);
                        hQ = (Q2 - g->Q2[iQ0]) / (g->Q2[iQ1] - g->Q2[iQ0]);
                }

                int i;
                const float * const lambda0 = g->lambda + iQ0 * g->nf;
                const float * const lambda1 = g->lambda + iQ1 * g->nf;
                const float * const f0 = g->f + iQ0 * g->nf;
                const float * const f1 = g->f + iQ1 * g->nf;
                for (i = 0; i < g->nf; i++) {
                        const float y0 = lambda0[i];
                        const float y1 = lambda1[i];
                        if ((y0 == 0.) || (y1 == 0.)) {
                                f[i] = f0[i] * (1. - hQ) + f1[i] * hQ;
                        } else {
                                const float lambda = y0 * (1. - hQ) + y1 * hQ;
                                f[i] = (f0[i] * (1. - hQ) + f1[i] * hQ) *
                                        pow(x / g->x[0], -lambda);
                        }
                }
        } else {
                /* We are within the table bounds. Let's locate the bracketing
                 * table rows.
                 */
                int ix0 = 0, ix1 = g->nx - 1, iQ0 = 0, iQ1 = g->nQ2 - 1;
                table_bracket(g->x, x, &ix0, &ix1);
                table_bracket(g->Q2, Q2, &iQ0, &iQ1);
                const float hx = (x - g->x[ix0]) / (g->x[ix1] - g->x[ix0]);
                const float hQ = (Q2 - g->Q2[iQ0]) / (g->Q2[iQ1] - g->Q2[iQ0]);

                /* Interpolate the PDFs. */
                int i;
                const float * const f00 =
                    g->f + (ix0 * g->nQ2 + iQ0) * g->nf;
                const float * const f01 =
                    g->f + (ix0 * g->nQ2 + iQ1) * g->nf;
                const float * const f10 =
                    g->f + (ix1 * g->nQ2 + iQ0) * g->nf;
                const float * const f11 =
                    g->f + (ix1 * g->nQ2 + iQ1) * g->nf;
                const float r00 = (1. - hx) * (1. - hQ);
                const float r01 = (1. - hx) * hQ;
                const float r10 = hx * (1. - hQ);
                const float r11 = hx * hQ;
                for (i = 0; i < g->nf; i++)
                        f[i] = r00 * f00[i] + r01 * f01[i] +
                               r10 * f10[i] + r11 * f11[i];
        }
}

static double rest_mass(enum ent_pid pid);

/* Load DIS physics from an ENT file. */
static enum ent_return physics_load(
    struct ent_physics ** physics, FILE * stream)
{
        *physics = NULL;
        enum ent_return rc;

        /* Check the version number. */
        int version;
        if ((fscanf(stream, "%d", &version) != 1) ||
            (version != ENT_FORMAT_VERSION)) {
                rc = ENT_RETURN_FORMAT_ERROR;
                goto error;
        }

        /* Skip metadata. */
        while (fgetc(stream) != 0x0) {
                if (feof(stream)) {
                        rc = ENT_RETURN_FORMAT_ERROR;
                        goto error;
                }
        }

        /* Create the physics wrapper. */
        *physics = physics_create();
        if (*physics == NULL) {
                rc = ENT_RETURN_MEMORY_ERROR;
                goto error;
        }

        /* Load SFs tables. */
        int i;
        for (i = 0; i < PROGET_N_DIS; i++) {
                if ((i == PROGET_NC_NU_BAR_NEUTRON) ||
                    (i == PROGET_NC_NU_BAR_PROTON)) {
                        /* Anti-neutrino NC case. */
                        const int ii = i - 2;
                        const int nx = (*physics)->sf[ii]->nx;
                        const int nQ2 = (*physics)->sf[ii]->nQ2;
                        if (((*physics)->sf[i] = grid_create(nx, nQ2, 3)) ==
                            NULL) {
                                rc = ENT_RETURN_MEMORY_ERROR;
                                goto error;
                        }

                        /* Copy data but changing F3 sign. */
                        memcpy((*physics)->sf[i]->x, (*physics)->sf[ii]->x,
                            (*physics)->sf[i]->nx *
                            sizeof(*(*physics)->sf[i]->x));
                        memcpy((*physics)->sf[i]->Q2, (*physics)->sf[ii]->Q2,
                            (*physics)->sf[i]->nQ2 *
                            sizeof(*(*physics)->sf[i]->Q2));

                        const float * sf0 = (*physics)->sf[ii]->f;
                        float * sf1 = (*physics)->sf[i]->f;
                        int j;
                        for (j = 0; j < nx * nQ2; j++, sf0 += 3, sf1 += 3) {
                                sf1[0] = sf0[0];
                                sf1[1] = -sf0[1];
                                sf1[2] = sf0[2];
                        }

                        /* Compute lambda exponents. */
                        if ((rc = grid_compute_lambda(
                            (*physics)->sf[i], NULL)) != ENT_RETURN_SUCCESS) {
                                goto error;
                        }
                } else {
                        if ((rc = grid_load((*physics)->sf + i, stream,
                            NULL)) !=
                            ENT_RETURN_SUCCESS) {
                                goto error;
                        }
                }
        }

        /* Load top production threshold. */
        if (fread((*physics)->mteff, sizeof(*(*physics)->mteff), 4, stream)
            != 4) {
                rc = ENT_RETURN_FORMAT_ERROR;
                goto error;
        }
        const double mt = rest_mass(ENT_PID_TOP);
        double mtmax = 0.;
        for (i = 0; i < 4; i++) {
                if ((*physics)->mteff[i] < mt) (*physics)->mteff[i] = mt;
                if ((*physics)->mteff[i] > mtmax) mtmax = (*physics)->mteff[i];
        }
        (*physics)->mteff[4] = mtmax;

        return ENT_RETURN_SUCCESS;
error:
        free(*physics);
        *physics = NULL;
        return rc;
}

/* Load a PDF from a .dat file in lhagrid1 format. */
static enum ent_return lha_load_table(
    struct file_buffer ** buffer, struct grid ** pdf_ptr, int * eof)
{
#define LHAPDF_NF_MAX 13
        /* Locate the next data segment. */
        for (;;) {
                file_get_line(buffer, 0);
                if (strlen((*buffer)->cursor) < 3) continue;
                if (((*buffer)->cursor[0] == '-') &&
                    ((*buffer)->cursor[1] == '-') &&
                    ((*buffer)->cursor[2] == '-'))
                        break;
        }
        long int pos = ftell((*buffer)->stream);

        {
                /* Check for end of file. */
                *eof = 0;
                char tmp[64];
                if ((fgets(tmp, 63, (*buffer)->stream) == NULL) ||
                    (feof((*buffer)->stream))) {
                        *eof = 1;
                        return ENT_RETURN_SUCCESS;
                } else if (fseek((*buffer)->stream, pos, SEEK_SET) != 0) {
                        return ENT_RETURN_IO_ERROR;
                }
        }

        /* Parse the table format. */
        file_get_line(buffer, 0);
        const int nx = file_count_words(*buffer);
        file_get_line(buffer, 0);
        const int nQ = file_count_words(*buffer);
        file_get_line(buffer, 0);
        const int nf = file_count_words(*buffer);
        if ((nf != LHAPDF_NF_MAX - 2) && (nf != LHAPDF_NF_MAX)) {
                return ENT_RETURN_FORMAT_ERROR;
        }

        /* Allocate and map the memory for the tables. */
        struct grid * pdf;
        if ((*pdf_ptr = grid_create(nx, nQ, nf)) == NULL) {
                return ENT_RETURN_MEMORY_ERROR;
        } else {
                pdf = *pdf_ptr;
        }

        /* Roll back in the file and copy the data to memory. */
        if (fseek((*buffer)->stream, pos, SEEK_SET) != 0) {
                return ENT_RETURN_IO_ERROR;
        }
        file_get_line(buffer, 0);

        float row[nf];
        int index[nf];
        file_get_table(buffer, nx, pdf->x);
        file_get_table(buffer, nQ, pdf->Q2);
        file_get_table(buffer, nf, row);

        int i;
        for (i = 0; i < nQ; i++) pdf->Q2[i] *= pdf->Q2[i];
        for (i = 0; i < nf; i++) {
                if (row[i] == 21.)
                        index[i] = (pdf->nf - 1) / 2;
                else
                        index[i] = (int)row[i] + (nf - 1) / 2;
                if ((index[i] < 0) || (index[i] >= nf)) {
                        return ENT_RETURN_FORMAT_ERROR;
                }
        }

        float * table;
        for (i = 0, table = pdf->f; i < nx * nQ; i++, table += nf) {
                int j;
                file_get_table(buffer, nf, row);
                for (j = 0; j < nf; j++) table[index[j]] = row[j];
        }

        /* Tabulate the lambda exponents for the small x extrapolation. */
        return grid_compute_lambda(pdf, NULL);

#undef LHAPDF_NF_MAX
}

/* Load all PDFs from a .dat file in lhagrid1 format. */
static enum ent_return lha_load(struct ent_physics ** physics, FILE * stream)
{
        *physics = NULL;
        enum ent_return rc;
        struct file_buffer * buffer = NULL;
        struct grid * head = NULL;

        /* Redirect any subsequent file parsing error. */
        jmp_buf env;
        if (setjmp(env) != 0) {
                rc = buffer->code;
                goto exit;
        }

        /* Create the temporary file buffer. */
        if ((rc = file_buffer_create(&buffer, stream, env)) !=
            ENT_RETURN_SUCCESS)
                goto exit;

        /* Load the PDF set(s). */
        struct grid * pdf = NULL;
        int eof;
        for (eof = 0; !eof;) {
                struct grid * g;
                if ((rc = lha_load_table(&buffer, &g, &eof)) !=
                    ENT_RETURN_SUCCESS)
                        goto exit;
                if (!eof) {
                        if (head == NULL) {
                                head = pdf = g;
                        } else if (pdf != NULL) {
                                pdf->next = g;
                                pdf = g;
                        }
                }
        }
        if (head == NULL) {
                rc = ENT_RETURN_FORMAT_ERROR;
                goto exit;
        }

        /* Create the physics wrapper. */
        *physics = physics_create();
        if (*physics == NULL) {
                rc = ENT_RETURN_MEMORY_ERROR;
        } else {
                (*physics)->pdf = head;
        }

        /* Set cutoff for top production. */
        const double mt = rest_mass(ENT_PID_TOP);
        int i;
        for (i = 0; i < 5; i++) {
                (*physics)->mteff[i] = mt;
        }
exit:
        free(buffer);
        if (rc != ENT_RETURN_SUCCESS) {
                grid_destroy(&head);
                free(*physics);
                *physics = NULL;
        }
        return rc;
}

/* Compute SFs for DIS. */
static void dis_compute_sf(struct ent_physics * physics,
    enum ent_pid projectile, double Z, double A, enum ent_process process,
    float x, float Q2, double * sf)
{
        if (physics->pdf != NULL) {
                /* Get the PDFs. */
                const struct grid * const pdf = physics->pdf;
                float xfx[pdf->nf];
                grid_compute(pdf, x, Q2, xfx);

                /* Compute the relevant structure functions. */
                const int eps = (projectile > 0) ? 1 : -1;
                const int nf = (pdf->nf - 1) / 2;
                const double N = A - Z; /* Number of neutrons. */
                if (process == ENT_PROCESS_DIS_CC_OTHER) {
                        /* Charged current DIS process (w/o bottom). */
                        const double d = (Z <= 0.) ? 0. : xfx[1 * eps + nf];
                        const double u = (N <= 0.) ? 0. : xfx[2 * eps + nf];
                        const double s = xfx[3 * eps + nf];
                        const double F1 = Z * d + N * u + A * s;

                        const double dbar = (N <= 0.) ? 0. : xfx[-1 * eps + nf];
                        const double ubar = (Z <= 0.) ? 0. : xfx[-2 * eps + nf];
                        const double cbar = xfx[-4 * eps + nf];
                        const double F2 = N * dbar + Z * ubar + A * cbar;

                        sf[0] = 2 * (F1 + F2);
                        sf[1] = 2 * (F1 - F2);
                        sf[2] = 0.;
                } else if (process == ENT_PROCESS_DIS_CC_TOP) {
                        /* Charged current DIS process (only bottom). */
                        const double b = xfx[5 * eps + nf];
                        const double F1 = A * b;

                        sf[0] = 2 * F1;
                        sf[1] = 2 * F1;
                        sf[2] = 0.;
                } else {
                        /* Neutral current DIS process. */
                        const double d = xfx[1 * eps + nf];
                        const double u = xfx[2 * eps + nf];
                        const double s = xfx[3 * eps + nf];
                        const double c = xfx[4 * eps + nf];
                        const double b = xfx[5 * eps + nf];

                        const double s23 = 2. * ENT_PHYS_SIN_THETA_W_2 / 3.;
                        const double gp2 = 1. + 4. * s23 * (s23 - 1.);
                        const double gm2 = 4. * s23 * s23;
                        const double gpp2 = 1. + s23 * (s23 - 2.);
                        const double gpm2 = s23 * s23;
                        const double F1 = Z * u + N * d + A * c;
                        const double F2 = Z * d + N * u + A * (s + b);

                        const double dbar = xfx[-1 * eps + nf];
                        const double ubar = xfx[-2 * eps + nf];
                        const double sbar = xfx[-3 * eps + nf];
                        const double cbar = xfx[-4 * eps + nf];
                        const double bbar = xfx[-5 * eps + nf];
                        const double F3 = Z * ubar + N * dbar + A * cbar;
                        const double F4 = Z * dbar + N * ubar +
                                          A * (sbar + bbar);

                        const double Fp1 = 0.5 * (gp2 * F1 + gpp2 * F2 +
                            gm2 * F3 + gpm2 * F4);
                        const double Fp2 = (gm2 * F1 + gpm2 * F2 + gp2 * F3 +
                            gpp2 * F4);

                        sf[0] = Fp1 + Fp2;
                        sf[1] = Fp1 - Fp2;
                        sf[2] = 0.;
                }
        } else {
                /* Combine the tabulated SFs */
                int itab0, itab1;
                if (process == ENT_PROCESS_DIS_CC_TOP) {
                        itab0 = 8 + (projectile < 0);
                        itab1 = itab0 + 2;
                } else {
                        itab0 = (process == ENT_PROCESS_DIS_NC) +
                            2 * (projectile < 0);
                        itab1 = itab0 + 4;
                }
                const double N = A - Z; /* Number of neutrons. */

                memset(sf, 0x0, 3 * sizeof(*sf));
                if (N > 0.) {
                        const struct grid * const g = physics->sf[itab0];
                        float tmp[3];
                        grid_compute(g, x, Q2, tmp);
                        int i;
                        for (i = 0; i < 3; i++) sf[i] += N * tmp[i];
                }
                if (Z > 0.) {
                        const struct grid * const g = physics->sf[itab1];
                        float tmp[3];
                        grid_compute(g, x, Q2, tmp);
                        int i;
                        for (i = 0; i < 3; i++) sf[i] += Z * tmp[i];
                }
        }
}

/* Particle rest mass for a given PID. */
static double rest_mass(enum ent_pid pid)
{
        const double mass[25] = {
            0., 4.67E-03, 2.16E-03, 93E-03, 1.27, 4.18, 172.76, 0., 0., 0., 0.,
            ENT_MASS_ELECTRON,
            0.,
            ENT_MASS_MUON,
            0.,
            ENT_MASS_TAU,
            0.,
            0., 0., 0., 0., 0., 0.,
            ENT_MASS_Z,
            ENT_MASS_W
        };

        int aid = abs(pid);
        if (aid < sizeof(mass) / sizeof(mass[0])) {
                return mass[aid];
        } else if ((pid == ENT_PID_NEUTRON) || (pid == ENT_PID_PROTON)) {
                return ENT_MASS_NUCLEON;
        } else {
                return 0.;
        }
}

/* DCS for Deep Inelastic Scattering (DIS) on nucleons. */
static double dcs_dis(struct ent_physics * physics, enum ent_pid projectile,
    double energy, double Z, double A, enum ent_process process, double x,
    double y)
{
        if (process == ENT_PROCESS_DIS_CC_TOP) {
                /* Check kinematic threshold. */
                if (y <= 0.) return 0.;
                const double mt = physics->mteff[projectile < 0];
                const double xmax =
                    1. - mt * mt / (2 * ENT_MASS_NUCLEON * energy * y);
                if ((xmax <= 0.) || (x >= xmax)) return 0.;
        }

        /* Compute the SFs. */
        const double Q2 = 2. * x * y * ENT_MASS_NUCLEON * energy;
        double sf[3];
        dis_compute_sf(physics, projectile, Z, A, process, x, Q2, sf);

        /* Compute the DCS factor. */
        const double MX2 = (process == ENT_PROCESS_DIS_NC) ?
            ENT_MASS_Z * ENT_MASS_Z : ENT_MASS_W * ENT_MASS_W;
        const double r = MX2 / (MX2 + Q2);
        const double y1 = 1. - y;
        const double c0 = 1. + y1 * y1;
        const double c1 = 1. - y1 * y1;
        const double c2 = -y * y;
        const double factor = r * r * (c0 * sf[0] + c1 * sf[1] + c2 * sf[2]);
        if (factor <= 0.) return 0.;

        /* Return the DCS. */
        return energy * factor * (ENT_PHYS_GF * ENT_PHYS_HBC) *
            (ENT_PHYS_GF * ENT_PHYS_HBC) * ENT_MASS_NUCLEON / (2 * M_PI);
}

/* DCS for elastic scattering on electrons. */
static double dcs_elastic(
    enum ent_pid projectile, double energy, double Z, double y)
{
        const double R = 2. * ENT_PHYS_SIN_THETA_W_2;
        const double L = 2. * ENT_PHYS_SIN_THETA_W_2 - 1.;
        const double MZ2 = ENT_MASS_Z * ENT_MASS_Z;
        const double rZ = MZ2 / (MZ2 + 2. * ENT_MASS_ELECTRON * energy * y);

        double factor;
        if (projectile == ENT_PID_NU_BAR_E) {
                const double tmp1 = R * rZ;
                const double a =
                    ENT_MASS_W * ENT_MASS_W - 2. * ENT_MASS_ELECTRON * energy;
                const double b = ENT_WIDTH_W * ENT_MASS_W;
                const double c = 2. * ENT_MASS_W * ENT_MASS_W / (a * a + b * b);
                const double d = rZ + c * a;
                const double e = -c * b;
                factor = tmp1 * tmp1 + (d * d + e * e) * (1. - y) * (1. - y);
        } else if (projectile == ENT_PID_NU_E) {
                const double tmp1 = R * (1. - y) * rZ;
                const double MW2 = ENT_MASS_W * ENT_MASS_W;
                const double rW =
                    MW2 / (MW2 + 2. * ENT_MASS_ELECTRON * energy * (1. - y));
                const double tmp2 = L * rZ + 2. * rW;
                factor = tmp1 * tmp1 + tmp2 * tmp2;
        } else {
                const double tmp = (projectile > 0) ?
                    (R * R * (1. - y) * (1. - y) + L * L) :
                    (L * L * (1. - y) * (1. - y) + R * R);
                factor = tmp * rZ * rZ;
        }

        return Z * energy * ENT_MASS_ELECTRON * (ENT_PHYS_GF * ENT_PHYS_HBC) *
            (ENT_PHYS_GF * ENT_PHYS_HBC) * factor / (2 * M_PI);
}

/* DCS for inelastic scattering on electrons. */
static double dcs_inverse(enum ent_pid projectile, double energy, double Z,
    enum ent_process process, double y)
{
        if ((projectile != ENT_PID_NU_BAR_E) &&
            ((projectile != ENT_PID_NU_MU) ||
                (process != ENT_PROCESS_INVERSE_MUON)) &&
            ((projectile != ENT_PID_NU_TAU) ||
                (process != ENT_PROCESS_INVERSE_TAU)))
                return 0.;

        const double ml = (process == ENT_PROCESS_INVERSE_MUON) ?
            ENT_MASS_MUON :
            ENT_MASS_TAU;
        const double MW2 = ENT_MASS_W * ENT_MASS_W;

        double factor;
        if (projectile == ENT_PID_NU_BAR_E) {
                const double a = 1. - 2. * ENT_MASS_ELECTRON * energy / MW2;
                const double b2 = ENT_WIDTH_W * ENT_WIDTH_W / MW2;
                factor = (1. - y) * (1. - y) / (a * a + b2);
        } else {
                const double tmp =
                    MW2 / (MW2 + 2. * ENT_MASS_ELECTRON * energy * (1. - y));
                factor = tmp * tmp;
        }

        const double r = 1. -
            (ml * ml - ENT_MASS_ELECTRON * ENT_MASS_ELECTRON) /
                (2. * ENT_MASS_ELECTRON * energy);
        return Z * energy * ENT_MASS_ELECTRON * (ENT_PHYS_GF * ENT_PHYS_HBC) *
            (ENT_PHYS_GF * ENT_PHYS_HBC) * 2. * r * r * factor / M_PI;
}

/* DCS for hadron(s) production by Glashow's resonance. */
static double dcs_glashow(
    enum ent_pid projectile, double energy, double Z, double y)
{
        if (projectile != ENT_PID_NU_BAR_E) return 0.;

        return dcs_inverse(projectile, energy, ENT_PROCESS_INVERSE_MUON, Z, y) *
            (1. - ENT_BR_W_TO_ELECTRON - ENT_BR_W_TO_MUON - ENT_BR_W_TO_TAU) /
            ENT_BR_W_TO_MUON;
}

/* Compute the DCS for a given process, projectile and target. */
static double dcs_compute(struct ent_physics * physics, enum ent_pid projectile,
    double energy, double Z, double A, enum ent_process process, double x,
    double y)
{
        if ((process == ENT_PROCESS_DIS_CC_OTHER) || 
            (process == ENT_PROCESS_DIS_CC_TOP) ||
            (process == ENT_PROCESS_DIS_NC)) {
                return dcs_dis(
                    physics, projectile, energy, Z, A, process, x, y);
        } else if (process == ENT_PROCESS_DIS_CC) {
                const double d0 = dcs_dis(physics, projectile, energy, Z, A,
                    ENT_PROCESS_DIS_CC_OTHER, x, y);
                const double d1 = dcs_dis(physics, projectile, energy, Z, A,
                    ENT_PROCESS_DIS_CC_TOP, x, y);
                return d0 + d1;
        } else if (process == ENT_PROCESS_ELASTIC) {
                return dcs_elastic(projectile, energy, Z, y);
        } else if ((process == ENT_PROCESS_INVERSE_MUON) ||
            (process == ENT_PROCESS_INVERSE_TAU)) {
                return dcs_inverse(projectile, energy, Z, process, y);
        } else if (process == ENT_PROCESS_GLASHOW_HADRON) {
                return dcs_glashow(projectile, energy, Z, y);
        } else {
                return -DBL_MAX;
        }
}

/* Kinematic limit for top prodcution in DIS. */
static void dis_top_clip_x(struct ent_physics * physics,
    enum proget_index proget, double energy, double y, double * xmax)
{
        const double mt = physics->mteff[proget - 8];
        double xk =
            1. - mt * mt / (2 * ENT_MASS_NUCLEON * energy * y);
        if (xk < 0.) xk = 0.;
        if (*xmax > xk) *xmax = xk;
}

/* Map DIS DDCS support as x, for fix y. */
static void dis_compute_support_x(struct ent_physics * physics, 
    enum ent_pid projectile, double energy, double Z, double A,
    enum ent_process process, double y, double * xlim)
{
        xlim[2] = 0.;

        double xmin = physics->dis_Q2min /
            (2 * ENT_MASS_NUCLEON * energy * y);
        double xmax = physics->dis_Q2max /
            (2 * ENT_MASS_NUCLEON * energy * y);
        if (xmax > 1.) xmax = 1.;

        if (process == ENT_PROCESS_DIS_CC_TOP) {
                dis_top_clip_x(physics, 8 + (projectile < 0), energy, y, &xmax);
        }

        if (xmin >= xmax) {
                xlim[0] = xmax;
                xlim[1] = xmax;
                return;
        }

        /* Find the lower bound of the support. First let us find an initial
         * bracketing using an exponential search.
         */
        double r, x0, x1;

        r = 2.;
        x0 = x1 = xmin;
        for (;;) {
                const double d = dcs_compute(physics, projectile, energy, Z,
                    A, process, x1, y);
                if (d > 0.) break;
                x0 = x1;
                x1 *= r;
                if (x1 >= xmax) {
                        /* Restart with a lower increment. */
                        x0 = x1 = xmin;
                        r = 0.5 * (1. + r);

                        if (r - 1. < 1E-02) {
                                xlim[0] = xmax;
                                xlim[1] = xmax;
                                return;
                        }
                }
        }
        if (x1 > x0) {
                /* Refine the bracketing with a binary search. */
                double tol = 1E-03 * xmin;
                if (tol > 1E-05) tol = 1E-05;
                while (x1 - x0 > tol) {
                        const double x2 = 0.5 * (x0 + x1);
                        const double d = dcs_compute(physics, projectile,
                            energy, Z, A, process, x2, y);
                        if (d > 0) {
                                x1 = x2;
                        } else {
                                x0 = x2;
                        }
                }
                xmin = x1;
        }

        /* Find the upper bound of the support. */
        r = 2.;
        x0 = x1 = xmax;
        for (;;) {
                const double d = dcs_compute(physics, projectile, energy, Z,
                    A, process, x0, y);
                if (d > 0.) break;
                x1 = x0;
                x0 /= r;
                if (x0 <= xmin) {
                        /* Restart with a lower increment. */
                        x0 = x1 = xmax;
                        r = 0.5 * (1. + r);
                }
        }
        if (x1 > x0) {
                /* Refine the bracketing with a binary search. */
                double tol = 1E-03 * xmax;
                if (tol > 1E-05) tol = 1E-05;
                while (x1 - x0 > tol) {
                        const double x2 = 0.5 * (x0 + x1);
                        const double d = dcs_compute(physics, projectile,
                            energy, Z, A, process, x2, y);
                        if (d > 0) {
                                x0 = x2;
                        } else {
                                x1 = x2;
                        }
                }
                xmax = x0;
        }

        xlim[0] = xmin;
        xlim[1] = xmax;
}

/* Interpolate x support for DIS DDCS. */
static void dis_get_support_x(struct ent_physics * physics,
    enum proget_index proget, double energy, double y, double * xmin,
    double * xmax, double * pdf_ratio)
{
        const double ymin = physics->dis_Q2min /
            (2 * ENT_MASS_NUCLEON * energy);
        if (ymin >= 1.) {
                *xmin = *xmax = 1.;
                if (pdf_ratio != NULL) *pdf_ratio = 0.;
                return;
        }

        const double dle = log(ENERGY_MAX / ENERGY_MIN) / (ENERGY_N - 1);
        const double dly = -log(ymin) / (DIS_Y_N - 1);

        double hy = log(y / ymin) / dly;
        int iy = (int)hy;
        if (iy >= DIS_Y_N - 1) {
                iy = DIS_Y_N - 2;
                hy = 1.;
        } else if (iy < 0) {
                iy = 0;
                hy = 0.;
        } else {
                hy -= iy;
        }

        double he = log(energy / ENERGY_MIN) / dle;
        int ie = (int)he;
        if (ie >= ENERGY_N - 1) {
                /* Use an extrapolation, for backward mode. */
                *xmin = physics->dis_Q2min /
                    (2 * ENT_MASS_NUCLEON * energy * y);
                *xmax = physics->dis_Q2max /
                    (2 * ENT_MASS_NUCLEON * energy * y);
                if (*xmax > 1.) *xmax = 1.;

                if (proget >= 8) {
                        /* Kinematic threshold. */
                        dis_top_clip_x(physics, proget, energy, y, xmax);
                }
                if (*xmin >= *xmax) *xmin = *xmax;

                if (pdf_ratio) {
                        /* However, freeze PDF ratio. */
                        const double * const xlim = physics->dis_xlim +
                            3 * ((proget * ENERGY_N + ENERGY_N - 1) * DIS_Y_N);

                        *pdf_ratio = xlim[3 * iy + 2] * (1. - hy) +
                            xlim[3 * (iy + 1) + 2] * hy;
                }
                return;
        } else if (ie < 0) {
                ie = 0;
                he = 0.;
        } else {
                he -= ie;
        }

        const double * const x00 = physics->dis_xlim +
            3 * ((proget * ENERGY_N + ie) * DIS_Y_N + iy);
        const double * const x01 = physics->dis_xlim +
            3 * ((proget * ENERGY_N + ie) * DIS_Y_N + iy + 1);
        const double * const x10 = physics->dis_xlim +
            3 * ((proget * ENERGY_N + ie + 1) * DIS_Y_N + iy);
        const double * const x11 = physics->dis_xlim +
            3 * ((proget * ENERGY_N + ie + 1) * DIS_Y_N + iy + 1);

        *xmin = x00[0] * (1. - he) * (1. - hy) + x01[0] * (1. - he) * hy +
                x10[0] * he * (1. - hy) + x11[0] * he * hy;
        *xmax = x00[1] * (1. - he) * (1. - hy) + x01[1] * (1. - he) * hy +
                x10[1] * he * (1. - hy) + x11[1] * he * hy;

        if (proget >= 8) {
                /* Kinematic threshold. */
                dis_top_clip_x(physics, proget, energy, y, xmax);
        }
        if (*xmin >= *xmax) *xmin = *xmax;

        if (pdf_ratio != NULL) {
                *pdf_ratio = x00[2] * (1. - he) * (1. - hy) +
                    x01[2] * (1. - he) * hy + x10[2] * he * (1. - hy) +
                    x11[2] * he * hy;
        }
}

/* Kinematic limit for top production. */
static void dis_top_clip_y(const struct ent_physics * physics,
    enum proget_index proget, double energy, double * ymin)
{
        const double mt = physics->mteff[proget - 8];
        const double yk = mt * mt / (2 * ENT_MASS_NUCLEON * energy);
        if (*ymin < yk) {
                *ymin = (yk < 1.) ? yk : 1.;
        }
}

/* Refine DIS DCS support, in forward case. */
static void dis_compute_support_y_forward(struct ent_physics * physics, 
    enum ent_pid projectile, double energy, double Z, double A,
    enum ent_process process, int imin, int imax, double * ylim)
{
        if ((imin >= imax) || (imax <= 0)) {
                ylim[0] = ylim[1] = 1.;
                return;
        }

        const double ymin = physics->dis_Q2min /
             (2 * ENT_MASS_NUCLEON * energy);
        const double dly = -log(ymin) / (DIS_Y_N - 1);

        /* Find the lower bound of the support. */
        double y0, y1, tol;
        enum proget_index proget;
        const enum ent_pid target = (Z > 0.) ? ENT_PID_PROTON : ENT_PID_NEUTRON;
        proget_compute(projectile, target, process, &proget);

        if (imin == 0) {
                y0 = ymin;
                y1 = ymin * exp(dly);
        } else {
                y0 = ymin * exp(imin * dly);
                y1 = ymin * exp((imin + 1) * dly);
        }

        /* Refine the bracketing with a binary search. */
        tol = 1E-04 * y0;
        if (tol > 1E-06) tol = 1E-06;
        while (y1 - y0 > tol) {
                const double y2 = 0.5 * (y0 + y1);
                double xmin, xmax;
                if (y2 < 1.) {
                        dis_get_support_x(physics, proget, energy, y2, &xmin,
                            &xmax, NULL);
                } else {
                        xmin = xmax = 1.;
                }
                if (xmin < xmax) {
                        y1 = y2;
                } else {
                        y0 = y2;
                }
        }
        if (process == ENT_PROCESS_DIS_CC_TOP) {
                dis_top_clip_y(physics, proget, energy, &y1);
        }
        ylim[0] = y1;

        /* Find the upper bound of the support. */
        if (imax >= DIS_Y_N - 1) {
                y0 = ymin * exp((DIS_Y_N - 2) * dly);
                y1 = 1.;
        } else {
                y0 = ymin * exp(imax * dly);
                y1 = ymin * exp((imax + 1) * dly);
        }

        /* Refine the bracketing with a binary search. */
        tol = 1E-04 * y1;
        if (tol > 1E-06) tol = 1E-06;
        while (y1 - y0 > tol) {
                const double y2 = 0.5 * (y0 + y1);
                double xmin, xmax;
                if (y2 < 1.) {
                        dis_get_support_x(physics, proget, energy, y2, &xmin,
                            &xmax, NULL);
                } else {
                        xmin = xmax = 1.;
                }
                if (xmin < xmax) {
                        y0 = y2;
                } else {
                        y1 = y2;
                }
        }
        ylim[1] = y0;
}

/* Build the interpolation or extrapolation factors. */
static int cross_section_prepare(double * cs, double energy, double ** cs0,
    double ** cs1, double * p1, double * p2)
{
        int mode;
        if (energy < ENERGY_MIN) {
                /* Log extrapolation model below Emin. */
                mode = 1;
                *cs0 = cs;
                *cs1 = cs + PROGET_N - 1;
                *p1 = (ENERGY_N - 1) / log(ENERGY_MAX / ENERGY_MIN);
                *p2 = energy / ENERGY_MIN;
        } else if (energy >= ENERGY_MAX) {
                /* Log extrapolation model above Emax. */
                mode = 2;
                *cs0 = cs + (ENERGY_N - 2) * (PROGET_N - 1);
                *cs1 = cs + (ENERGY_N - 1) * (PROGET_N - 1);
                *p1 = (ENERGY_N - 1) / log(ENERGY_MAX / ENERGY_MIN);
                *p2 = energy / ENERGY_MAX;
        } else {
                /* Interpolation model. */
                mode = 0;
                const double dle =
                    log(ENERGY_MAX / ENERGY_MIN) / (ENERGY_N - 1);
                *p1 = log(energy / ENERGY_MIN) / dle;
                const int i0 = (int)(*p1);
                *p1 -= i0;
                *cs0 = cs + i0 * (PROGET_N - 1);
                *cs1 = cs + (i0 + 1) * (PROGET_N - 1);
                *p2 = 0.;
        }
        return mode;
}

/* Low level routine for computing a specific cross-section by interpolation
 * or by extrapolation.
 */
static double cross_section_compute(
    int mode, int proget, double * cs0, double * cs1, double p1, double p2)
{
        if (mode == 0) {
                /* interpolation case. */
                if ((cs0[proget] > 0.) && (cs1[proget] > 0.)) {
                        return exp(log(cs0[proget]) * (1. - p1) +
                                   log(cs1[proget]) * p1);
                } else {
                        return cs0[proget] * (1. - p1) + cs1[proget] * p1;
                }
        } else if (mode == 1) {
                /* Extrapolation case (below). */
                if (cs0[proget] * cs1[proget] <= 0.) {
                        return 0.;
                } else {
                        const double a = log(cs1[proget] / cs0[proget]) * p1;
                        return cs0[proget] * pow(p2, a);
                }
        } else {
                /* Extrapolation case (above). */
                const double a = log(cs1[proget] / cs0[proget]) * p1;
                return cs1[proget] * pow(p2, a);
        }
}

/* Map DIS DCS support as Ef, in backward case. */
static void dis_compute_support_y_backward(struct ent_physics * physics,
    enum ent_pid projectile, double energy, double Z, double A,
    enum ent_process process, int mode, double * ylim)
{
#define HADRON_MAX_RATIO 1E+04

        enum proget_index proget;
        const enum ent_pid target = (Z > 0.) ? ENT_PID_PROTON : ENT_PID_NEUTRON;
        proget_compute(projectile, target, process, &proget);

        double Q2min;
        if (process == ENT_PROCESS_DIS_CC_TOP) {
                const double mt = physics->mteff[proget - 8];
                Q2min = mt * mt;
        } else {
                Q2min = physics->dis_Q2min;
        }

        double ymin;
        if (mode) {
                ymin = Q2min / (Q2min + 2 * ENT_MASS_NUCLEON * energy);
        } else {
                if (energy <= Q2min / (2 * ENT_MASS_NUCLEON)) {
                        ylim[0] = 1.;
                        ylim[1] = 1.;
                        return;
                } else {
                        ymin = Q2min /
                            (2 * ENT_MASS_NUCLEON * energy * HADRON_MAX_RATIO);
                        if (ymin > 1.) ymin = 1.;
                }
        }
        double ymax = 1.;

        /* Find the lower bound of the support. First let us find an initial
         * bracketing using an exponential search.
         */
        double r, y0, y1;

        r = 2.;
        y0 = y1 = ymin;
        for (;;) {
                if (y1 < 1.) {
                        const double Ei = mode ?
                            energy / (1. - y1) : energy / y1;
                        double *csl, *csh;
                        double pl, ph;
                        const int m = cross_section_prepare(
                            physics->cs_k, Ei, &csl, &csh, &pl, &ph);
                        const double cs = cross_section_compute(
                            m, proget, csl, csh, pl, ph);
                        if (cs > 0.) break;
                }
                y0 = y1;
                y1 *= r;
                if (y1 >= ymax) {
                        /* Restart with a lower increment. */
                        y0 = y1 = ymin;
                        r = 0.5 * (1. + r);

                        if (r - 1. < 1E-04) {
                                ylim[0] = ymax;
                                ylim[1] = ymax;
                                return;
                        }
                }
        }
        if (y1 > y0) {
                /* Refine the bracketing with a binary search. */
                double tol = 1E-04 * ymin;
                if (tol > 1E-06) tol = 1E-06;
                while (y1 - y0 > tol) {
                        const double y2 = 0.5 * (y0 + y1);
                        double cs;
                        if (y2 < 1.) {
                                const double Ei = mode ?
                                    energy / (1. - y2) : energy / y2;
                                double *csl, *csh;
                                double pl, ph;
                                const int m = cross_section_prepare(
                                    physics->cs_k, Ei, &csl, &csh, &pl, &ph);
                                cs = cross_section_compute(
                                    m, proget, csl, csh, pl, ph);
                        } else {
                                cs = 0.;
                        }
                        if (cs > 0.) {
                                y1 = y2;
                        } else {
                                y0 = y2;
                        }
                }
                ymin = y1;
        }

        /* Find the upper bound of the support. */
        r = 2.;
        y0 = y1 = ymax;
        for (;;) {
                if (y0 < 1.) {
                        const double Ei = mode ?
                            energy / (1. - y0) : energy / y0;
                        double *csl, *csh;
                        double pl, ph;
                        const int m = cross_section_prepare(
                            physics->cs_k, Ei, &csl, &csh, &pl, &ph);
                        const double cs = cross_section_compute(
                            m, proget, csl, csh, pl, ph);
                        if (cs > 0.) break;
                }
                y1 = y0;
                y0 /= r;
                if (y0 <= ymin) {
                        /* Restart with a lower increment. */
                        y0 = y1 = ymax;
                        r = 0.5 * (1. + r);
                }
        }
        if (y1 > y0) {
                /* Refine the bracketing with a binary search. */
                double tol = 1E-04 * ymax;
                if (tol > 1E-06) tol = 1E-06;
                while (y1 - y0 > tol) {
                        const double y2 = 0.5 * (y0 + y1);
                        double cs;
                        if (y2 < 1.) {
                                const double Ei = mode ?
                                    energy / (1. - y2) : energy / y2;
                                double *csl, *csh;
                                double pl, ph;
                                const int m = cross_section_prepare(
                                    physics->cs_k, Ei, &csl, &csh, &pl, &ph);
                                cs = cross_section_compute(
                                    m, proget, csl, csh, pl, ph);
                        } else {
                                cs = 0.;
                        }
                        if (cs > 0.) {
                                y0 = y2;
                        } else {
                                y1 = y2;
                        }
                }
                ymax = y0;
        }

        ylim[0] = ymin;
        ylim[1] = ymax;
}

/* Interpolate y support for DIS DCS. */
static void dis_get_support_y(struct ent_physics * physics,
    enum proget_index proget, double energy, int mode, double * ymin,
    double * ymax)
{
        if (mode == -1) {
                /* Kinematic limit for top production. */
                const double mt = physics->mteff[proget - 8];
                if (energy <= mt * mt / (2 * ENT_MASS_NUCLEON)) {
                        *ymin = 1.;
                        *ymax = 1.;
                        return;
                }
        }

        const double dle = log(ENERGY_MAX / ENERGY_MIN) / (ENERGY_N - 1);

        double he = log(energy / ENERGY_MIN) / dle;
        int ie = (int)he;
        if (ie >= ENERGY_N - 1) {
                if (mode == 1) {
                        *ymin = physics->dis_Q2min /
                            (2 * ENT_MASS_NUCLEON * energy );
                        if (proget >= 8)
                                dis_top_clip_y(physics, proget, energy, ymin);
                } else if (mode == 0) {
                        double Q2min;
                        if (proget >= 8) {
                                const double mt = physics->mteff[proget - 8];
                                Q2min = mt * mt;
                        } else {
                                Q2min = physics->dis_Q2min;
                        }
                        *ymin = Q2min / (Q2min + 2 * ENT_MASS_NUCLEON * energy);
                } else {
                        *ymin = physics->dis_Q2min /
                            (2 * ENT_MASS_NUCLEON * energy * HADRON_MAX_RATIO);
                }
                *ymax = 1.;
                return;
        } else if (ie < 0) {
                ie = 0;
                he = 0.;
        } else {
                he -= ie;
        }

        int offset = 2 * proget * ENERGY_N;
        double * dis_ylim;
        if (mode == 1) dis_ylim = physics->dis_ylim_f;
        else if (mode == 0) dis_ylim = physics->dis_ylim_b;
        else dis_ylim = physics->dis_ylim_h;
        double * y0 = dis_ylim + offset + 2 * ie;
        double * y1 = dis_ylim + offset + 2 * (ie + 1);

        *ymin = (y0[0] == y1[0]) ? y0[0] : y0[0] * (1. - he) + y1[0] * he;
        *ymax = (y0[1] == y1[1]) ? y0[1] : y0[1] * (1. - he) + y1[1] * he;

        if (mode == 1) {
                if (proget < 8) {
                        const double yk = physics->dis_Q2min /
                            (2 * ENT_MASS_NUCLEON * energy);
                        if (*ymin < yk) *ymin = yk;
                } else {
                        dis_top_clip_y(physics, proget, energy, ymin);
                }
        } else if (mode == 0) {
                double Q2min;
                if (proget < 8) {
                        Q2min = physics->dis_Q2min;
                } else {
                        const double mt = physics->mteff[proget - 8];
                        Q2min = mt * mt;
                }
                const double yk =
                    Q2min / (Q2min + 2 * ENT_MASS_NUCLEON * energy);
                if (*ymin < yk) *ymin = yk;
        }

#undef HADRON_MAX_RATIO
}

/* Compute the total cross-section for a process using a Gaussian
 * quadrature.
 */
static double dcs_integrate(struct ent_physics * physics,
    enum ent_pid projectile, double energy, double Z, double A,
    enum ent_process process, double * pdf, double * cdf, double * xlim)
{
/*
 * Coefficients for the Gaussian quadrature from:
 * https://pomax.github.io/bezierinfo/legendre-gauss.html.
 */
#define N_GQ 9
        const double xGQ[N_GQ] = { 0.0000000000000000, -0.8360311073266358,
                0.8360311073266358, -0.9681602395076261, 0.9681602395076261,
                -0.3242534234038089, 0.3242534234038089, -0.6133714327005904,
                0.6133714327005904 };
        const double wGQ[N_GQ] = { 0.3302393550012598, 0.1806481606948574,
                0.1806481606948574, 0.0812743883615744, 0.0812743883615744,
                0.3123470770400029, 0.3123470770400029, 0.2606106964029354,
                0.2606106964029354 };

        if ((process == ENT_PROCESS_DIS_CC_OTHER) ||
            (process == ENT_PROCESS_DIS_CC_TOP) ||
            (process == ENT_PROCESS_DIS_NC)) {
                /* Deep Inelastic Scattering requires a double integral. */
                const enum ent_pid target =
                    (Z > 0.) ? ENT_PID_PROTON : ENT_PID_NEUTRON;
                enum proget_index proget;
                proget_compute(projectile, target, process, &proget);

                double ymin, ymax;
                dis_get_support_y(physics, proget, energy, 1, &ymin, &ymax);
                if (ymin / ymax - 1. >= -FLT_EPSILON) {
                        memset(pdf, 0x0, DIS_Y_N * sizeof(*pdf));
                        memset(cdf, 0x0, DIS_Y_N * sizeof(*pdf));
                        return 0.;
                }

                const double dly = log(ymax / ymin) / (DIS_Y_N - 1);
                double ratio[DIS_Y_N];

                int i;
                for (i = 0; i < DIS_Y_N; i++) {
                        const double y = ymin * exp(i * dly);

                        double xmin, xmax;
                        dis_get_support_x(
                            physics, proget, energy, y, &xmin, &xmax, NULL);
                        if (xmin >= xmax) {
                                pdf[i] = 0.;
                                continue;
                        }

                        const double ln10i = 0.4343; /* 1 / ln(10) */
                        double dlx = log(xmax / xmin);
                        const int n_pts_per_decade = 20;
                        int nx = (int)(dlx * n_pts_per_decade * ln10i / N_GQ);
                        if (nx < 4) nx = 4;
                        dlx /= nx;

                        const double MX = (process == ENT_PROCESS_DIS_NC) ?
                            ENT_MASS_Z : ENT_MASS_W;
                        const double x0 = 0.5 * MX * MX /
                            (energy * y * ENT_MASS_NUCLEON);
                        const double beta = DIS_X_EXPONENT;

                        double J = 0., rmax = 0.;
                        int j;
                        for (j = 0; j < nx * N_GQ; j++) {
                                const double x = xmin *
                                    exp((0.5 + 0.5 * xGQ[j % N_GQ] + j / N_GQ) *
                                        dlx);

                                const double d = dcs_compute(physics,
                                    projectile, energy, Z, A, process, x, y);
                                const double r = d / pow(1. + x / x0, -beta);
                                if (r > rmax) rmax = r;

                                J += wGQ[j % N_GQ] * d * x;
                        }
                        pdf[i] = 0.5 * J * dlx;

                        if (rmax > 0.) {
                                const double bmin =
                                    pow(1. + xmin / x0, 1. - beta);
                                const double bmax =
                                    pow(1. + xmax / x0, 1. - beta);

                                ratio[i] = rmax  * x0 * (bmin - bmax) /
                                    ((beta - 1) * pdf[i]);
                        } else {
                                ratio[i] = 0.;
                        }
                }

                double cs = 0., y0 = ymin;
                cdf[0] = 0.;
                for (i = 0; i < DIS_Y_N - 1; i++) {
                        const double y1 = ymin * exp((i + 1) * dly);
                        cs += 0.5 * (pdf[i] * y0 + pdf[i + 1] * y1) * dly;
                        cdf[i + 1] = cs;
                        y0 = y1;
                }

                if (cs > 0.) {
                        const double csi = 1. / cs;
                        for (i = 1; i < DIS_Y_N; i++) {
                                pdf[i] *= csi;
                                cdf[i] *= csi;
                        }
                }

                /* Interpolate ratio values over the large y-grid. */
                const double ymin0 = physics->dis_Q2min /
                    (2 * ENT_MASS_NUCLEON * energy);
                const double dly0 = -log(ymin0) / (DIS_Y_N - 1);
                for (i = 0; i < DIS_Y_N; i++) {
                        const double y0 = ymin0 * exp(i * dly0);
                        double h1 = log(y0 / ymin) / dly;
                        int i1 = (int)h1;
                        if (i1 > DIS_Y_N - 1) {
                                xlim[3 * i + 2] = 0.;
                                continue;
                        } else if (i1 == DIS_Y_N) {
                                h1 = 1.;
                                i1 = DIS_Y_N - 2;
                        } else if (i1 < 0) {
                                xlim[3 * i + 2] = 0.;
                                continue;
                        } else {
                                h1 -= i1;
                        }

                        xlim[3 * i + 2] = ratio[i1] * (1. - h1) +
                            ratio[i1 + 1] * h1;
                }

                return cs;
        } else {
                /* Scattering off an atomic electron. */
                double mu = ENT_MASS_ELECTRON;
                if (process == ENT_PROCESS_INVERSE_MUON)
                        mu = ENT_MASS_MUON;
                else if (process == ENT_PROCESS_INVERSE_TAU)
                        mu = ENT_MASS_TAU;
                if (energy <= mu) return 0.;
                const double ymin = mu / energy;
                const int ny = 6;
                const double dly = log(0.5 / ymin) / ny;
                double I = 0.;
                int i;
                for (i = 0; i < ny * N_GQ; i++) {
                        const double y = ymin *
                            exp((0.5 + 0.5 * xGQ[i % N_GQ] + i / N_GQ) * dly);
                        I += wGQ[i % N_GQ] * y *
                            dcs_compute(physics, projectile, energy, Z, A,
                                process, 0., y);
                }

                const double y1min = mu / energy;
                const double dly1 = log(0.5 / y1min) / ny;
                double I1 = 0.;
                for (i = 0; i < ny * N_GQ; i++) {
                        const double y1 = y1min *
                            exp((0.5 + 0.5 * xGQ[i % N_GQ] + i / N_GQ) * dly1);
                        I1 += wGQ[i % N_GQ] * y1 *
                            dcs_compute(physics, projectile, energy, Z, A,
                                process, 0., 1. - y1);
                }

                return 0.5 * (I * dly + I1 * dly1);
        }
#undef N_GQ
}

/* Tabulate cross-sections etc. */
static void physics_tabulate(struct ent_physics * physics)
{
        const enum ent_process process[PROGET_N - 1] = {
                ENT_PROCESS_DIS_CC_OTHER, ENT_PROCESS_DIS_NC,
                ENT_PROCESS_DIS_CC_OTHER, ENT_PROCESS_DIS_NC,
                ENT_PROCESS_DIS_CC_OTHER, ENT_PROCESS_DIS_NC,
                ENT_PROCESS_DIS_CC_OTHER, ENT_PROCESS_DIS_NC,
                ENT_PROCESS_DIS_CC_TOP, ENT_PROCESS_DIS_CC_TOP,
                ENT_PROCESS_DIS_CC_TOP, ENT_PROCESS_DIS_CC_TOP,
                ENT_PROCESS_ELASTIC, ENT_PROCESS_ELASTIC,
                ENT_PROCESS_ELASTIC, ENT_PROCESS_ELASTIC, ENT_PROCESS_ELASTIC,
                ENT_PROCESS_ELASTIC, ENT_PROCESS_INVERSE_MUON,
                ENT_PROCESS_INVERSE_TAU, ENT_PROCESS_INVERSE_MUON,
                ENT_PROCESS_INVERSE_TAU };
        const enum ent_pid projectile[PROGET_N - 1] = {
                ENT_PID_NU_E, ENT_PID_NU_E, ENT_PID_NU_BAR_E, ENT_PID_NU_BAR_E,
                ENT_PID_NU_E, ENT_PID_NU_E, ENT_PID_NU_BAR_E, ENT_PID_NU_BAR_E,
                ENT_PID_NU_E, ENT_PID_NU_BAR_E, ENT_PID_NU_E, ENT_PID_NU_BAR_E,
                ENT_PID_NU_E, ENT_PID_NU_BAR_E, ENT_PID_NU_MU,
                ENT_PID_NU_BAR_MU, ENT_PID_NU_TAU, ENT_PID_NU_BAR_TAU,
                ENT_PID_NU_MU, ENT_PID_NU_TAU, ENT_PID_NU_BAR_E,
                ENT_PID_NU_BAR_E };
        const double Z[PROGET_N - 1] = {
                0., 0., 0., 0., 1., 1., 1., 1., 0., 0., 1., 1.,
                1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
        };

        /* Get the PDF / SFs bounds. */
        double Q2min = DBL_MAX, Q2max = -DBL_MAX;
        if (physics->pdf != NULL) {
                struct grid * g;
                for (g = physics->pdf; g != NULL; g = g->next) {
                        if (g->Q2[0] < Q2min) Q2min = g->Q2[0];
                        const double tmp = g->Q2[g->nQ2 - 1];
                        if (tmp > Q2max) Q2max = tmp;
                }
        }

        int i;
        for (i = 0; i < PROGET_N_DIS; i++) if (physics->sf[i] != NULL) {
                struct grid * g = physics->sf[i];
                const int iq = g->nQ2 - 1;
                if (g->Q2[0] < Q2min) Q2min = g->Q2[0];
                if (g->Q2[iq] > Q2max) Q2max = g->Q2[iq];
        }

        physics->dis_Q2min = Q2min;
        physics->dis_Q2max = Q2max;

        /* Tabulate DIS supports in forward mode. */
        const double dlE = log(ENERGY_MAX / ENERGY_MIN) / (ENERGY_N - 1);
        for (i = 0; i < PROGET_N_DIS; i++) {
                int kmin[ENERGY_N], kmax[ENERGY_N];

                int j;
                for (j = 0; j < ENERGY_N; j++) {
                        const double energy = ENERGY_MIN * exp(j * dlE);
                        const double ymin = physics->dis_Q2min /
                             (2 * ENT_MASS_NUCLEON * energy);

                        const double dly = -log(ymin) / (DIS_Y_N - 1);
                        double *  xlim = physics->dis_xlim +
                            3 * (i * ENERGY_N + j) * DIS_Y_N;

                        kmin[j] = kmax[j] = 0;
                        int k, inside = 0;
                        for (k = 0; k < DIS_Y_N; k++, xlim += 3) {
                                const double y = ymin * exp(k * dly);
                                dis_compute_support_x(physics,
                                    projectile[i], energy, Z[i], 1.,
                                    process[i], y, xlim);

                                if (inside) {
                                        if (xlim[0] < xlim[1]) {
                                                kmax[j] = k;
                                        } else {
                                                inside = 0;
                                        }
                                } else {
                                        if (xlim[0] < xlim[1]) {
                                                kmin[j] = k;
                                                inside = 1;
                                        }
                                }
                        }
                }

                for (j = 0; j < ENERGY_N; j++) {
                        const double energy = ENERGY_MIN * exp(j * dlE);
                        double * ylim = physics->dis_ylim_f +
                            2 * (i * ENERGY_N + j);

                        dis_compute_support_y_forward(physics, projectile[i],
                            energy, Z[i], 1., process[i], kmin[j], kmax[j],
                            ylim);
                }
        }

        /* Tabulate cross-sections. */
        double * cs_k;
        const int np = PROGET_N - 1;
        for (i = 0, cs_k = physics->cs_k; i < ENERGY_N; i++, cs_k += np) {
                const double energy = ENERGY_MIN * exp(i * dlE);

                int j;
                for (j = 0; j < np; j++) {
                        double * pdf, * cdf, * xlim;
                        if (j < PROGET_N_DIS) {
                                pdf = physics->dis_pdf +
                                    (j * ENERGY_N + i) * DIS_Y_N;
                                cdf = physics->dis_cdf +
                                    (j * ENERGY_N + i) * DIS_Y_N;
                                xlim = physics->dis_xlim +
                                    3 * (j * ENERGY_N + i) * DIS_Y_N;
                        } else {
                                pdf = cdf = xlim = NULL;
                        }

                        /* Integrate. */
                        cs_k[j] = dcs_integrate(physics, projectile[j],
                            energy, Z[j], 1., process[j], pdf, cdf, xlim);
                }
        }

        /* Tabulate DIS DCS supports in backward mode. */
        for (i = 0; i < PROGET_N_DIS; i++) {
                double * ylim_b = physics->dis_ylim_b + 2 * i * ENERGY_N;
                double * ylim_h = physics->dis_ylim_h + 2 * i * ENERGY_N;

                int j;
                for (j = 0; j < ENERGY_N; j++, ylim_b += 2, ylim_h += 2) {
                        const double energy = ENERGY_MIN * exp(j * dlE);
                        dis_compute_support_y_backward(physics, projectile[i],
                            energy, Z[i], 1., process[i], 1, ylim_b);
                        dis_compute_support_y_backward(physics, projectile[i],
                            energy, Z[i], 1., process[i], 0, ylim_h);
                }
        }
}

/* Load cross-sections from a text file. */
static enum ent_return cs_load(FILE * stream, struct ent_physics * physics)
{
#define BUFFER_SIZE 2048

        char buffer[BUFFER_SIZE];
        int i = 0;
        double * table = physics->cs_t;
        double * cs_k = physics->cs_k;
        const double re = log(ENERGY_MAX / ENERGY_MIN) / (ENERGY_N - 1);
        while (!feof(stream)) {
                if (fgets(buffer, BUFFER_SIZE - 1, stream) == NULL)
                        return ENT_RETURN_FORMAT_ERROR;
                if (buffer[0] == '#') continue;

                double data[9];
                const int nread = sscanf(buffer,
                    "%lf %lf %lf %lf %lf %lf %lf %lf %lf", data, data + 1,
                    data + 2, data + 3, data + 4, data + 5, data + 6, data + 7,
                    data + 8);
                if (nread != 9) return ENT_RETURN_FORMAT_ERROR;

                const double energy = ENERGY_MIN * exp(i * re);
                if (fabs(data[0] / energy - 1.) > 1E-03)
                        return ENT_RETURN_FORMAT_ERROR;

                const double r2 = data[2] /
                    (cs_k[PROGET_CC_OTHER_NU_NEUTRON] +
                     cs_k[PROGET_CC_TOP_NU_NEUTRON]);
                table[PROGET_CC_OTHER_NU_NEUTRON] =
                    cs_k[PROGET_CC_OTHER_NU_NEUTRON] * r2;
                table[PROGET_CC_TOP_NU_NEUTRON] =
                    cs_k[PROGET_CC_TOP_NU_NEUTRON] * r2;

                table[PROGET_NC_NU_NEUTRON] = data[4];

                const double r6 = data[6] /
                    (cs_k[PROGET_CC_OTHER_NU_BAR_NEUTRON] +
                     cs_k[PROGET_CC_TOP_NU_BAR_NEUTRON]);
                table[PROGET_CC_OTHER_NU_BAR_NEUTRON] =
                    cs_k[PROGET_CC_OTHER_NU_BAR_NEUTRON] * r6;
                table[PROGET_CC_TOP_NU_BAR_NEUTRON] =
                    cs_k[PROGET_CC_TOP_NU_BAR_NEUTRON] * r6;

                table[PROGET_NC_NU_BAR_NEUTRON] = data[8];

                const double r1 = data[1] /
                    (cs_k[PROGET_CC_OTHER_NU_PROTON] +
                     cs_k[PROGET_CC_TOP_NU_PROTON]);
                table[PROGET_CC_OTHER_NU_PROTON] =
                    cs_k[PROGET_CC_OTHER_NU_PROTON] * r1;
                table[PROGET_CC_TOP_NU_PROTON] =
                    cs_k[PROGET_CC_TOP_NU_PROTON] * r1;

                table[PROGET_NC_NU_PROTON] = data[3];

                const double r5 = data[5] /
                    (cs_k[PROGET_CC_OTHER_NU_BAR_PROTON] +
                     cs_k[PROGET_CC_TOP_NU_BAR_PROTON]);
                table[PROGET_CC_OTHER_NU_BAR_PROTON] =
                    cs_k[PROGET_CC_OTHER_NU_BAR_PROTON] * r5;
                table[PROGET_CC_TOP_NU_BAR_PROTON] =
                    cs_k[PROGET_CC_TOP_NU_BAR_PROTON] * r5;

                table[PROGET_NC_NU_BAR_PROTON] = data[7];

                /* Use kinematic cross-sections for other processes. */
                int j;
                for (j = PROGET_N_DIS; j < PROGET_N - 1; j++)
                        table[j] = cs_k[j];

                table += PROGET_N - 1;
                cs_k += PROGET_N - 1;
                i++;
                if (i == ENERGY_N) break;
        }

        return ENT_RETURN_SUCCESS;

#undef BUFFER_SIZE
}

/* API constructor for a Physics object. */
enum ent_return ent_physics_create(struct ent_physics ** physics,
    const char * pdf_or_sf_file, const char * cs_file)
{
        ENT_ACKNOWLEDGE(ent_physics_create);

        enum ent_return rc;
        *physics = NULL;

        /* Load DIS SFs or PDFs. */
        FILE * stream;
        rc = ENT_RETURN_PATH_ERROR;
        if ((stream = fopen(pdf_or_sf_file, "rb")) == NULL) goto exit;

        char tag[sizeof(ENT_FORMAT_TAG)] = {0x0};
        fread(tag, sizeof(tag) - 1, 1, stream);
        if (strcmp(tag, ENT_FORMAT_TAG) == 0) {
                /* This is an ENT file. */
                if ((rc = physics_load(physics, stream)) != ENT_RETURN_SUCCESS)
                        goto exit;
        } else {
                /* This must be a PDF file in LHA format. */
                fseek(stream, 0, SEEK_SET);
                if ((rc = lha_load(physics, stream)) != ENT_RETURN_SUCCESS)
                        goto exit;
        }
        fclose(stream);
        stream = NULL;

        /* Build physics tabulations. */
        rc = ENT_RETURN_SUCCESS;
        physics_tabulate(*physics);

        if (cs_file != NULL) {
                /* Use provided cross-sections for the transport. */
                rc = ENT_RETURN_PATH_ERROR;
                if ((stream = fopen(cs_file, "r")) == NULL) goto exit;

                if ((rc = cs_load(stream, *physics)) != ENT_RETURN_SUCCESS)
                    goto exit;
                fclose(stream);
                stream = NULL;
        } else {
                /* Use computed cross-sections for the transport. */
                memcpy((*physics)->cs_t, (*physics)->cs_k,
                    ENERGY_N * (PROGET_N - 1) * sizeof(*(*physics)->cs_t));
        }

exit:
        if (stream != NULL) fclose(stream);
        ENT_RETURN(rc);
}

/* API destructor for a Physics object. */
void ent_physics_destroy(struct ent_physics ** physics)
{
        if ((physics == NULL) || (*physics == NULL)) return;

        struct grid * g = (*physics)->pdf;
        while (g != NULL) {
                struct grid * next = g->next;
                free(g);
                g = next;
        }

        int i;
        for (i = 0; i < PROGET_N_DIS; i++) {
                free((*physics)->sf[i]);
        }

        free(*physics);
        *physics = NULL;
}

static enum ent_return transport_cross_section(struct ent_physics * physics,
    enum ent_pid projectile, double energy, double Z, double A,
    enum ent_process process, double * cs);

/* Generic API function for accessing a specific total cross-section. */
enum ent_return ent_physics_cross_section(struct ent_physics * physics,
    enum ent_pid projectile, double energy, double Z, double A,
    enum ent_process process, double * cross_section)
{
        ENT_ACKNOWLEDGE(ent_physics_cross_section);
        *cross_section = 0.;

        enum ent_return rc;
        if (process == ENT_PROCESS_NONE) {
                /* Let us compute the total cross-section. */
                double cs[PROGET_N];
                if ((rc = transport_cross_section(physics, projectile, energy,
                         Z, A, ENT_PROCESS_NONE, cs)) != ENT_RETURN_SUCCESS)
                        ENT_RETURN(rc);
                *cross_section = cs[PROGET_N - 1];
                return ENT_RETURN_SUCCESS;
        }

        /* Build the interpolation or extrapolation factors for specific
         * cross-section.
         */
        double *cs0, *cs1;
        double p1, p2;
        int mode = cross_section_prepare(
            physics->cs_t, energy, &cs0, &cs1, &p1, &p2);

        /* Compute the relevant process-target indices. */
        if ((process == ENT_PROCESS_DIS_CC) ||
            (process == ENT_PROCESS_DIS_CC_OTHER) ||
            (process == ENT_PROCESS_DIS_CC_TOP)) {
                int proget;
                if (Z > 0.) {
                        if ((process == ENT_PROCESS_DIS_CC) ||
                            (process == ENT_PROCESS_DIS_CC_OTHER)) {
                                if ((rc = proget_compute(projectile,
                                    ENT_PID_PROTON, ENT_PROCESS_DIS_CC_OTHER,
                                    &proget)) != ENT_RETURN_SUCCESS)
                                        ENT_RETURN(rc);
                                *cross_section += Z *
                                    cross_section_compute(
                                        mode, proget, cs0, cs1, p1, p2);
                        }
                        if ((process == ENT_PROCESS_DIS_CC) ||
                            (process == ENT_PROCESS_DIS_CC_TOP)) {
                                if ((rc = proget_compute(projectile,
                                    ENT_PID_PROTON, ENT_PROCESS_DIS_CC_TOP,
                                    &proget)) != ENT_RETURN_SUCCESS)
                                        ENT_RETURN(rc);
                                *cross_section += Z *
                                    cross_section_compute(
                                        mode, proget, cs0, cs1, p1, p2);
                        }
                }
                const double N = A - Z;
                if (N > 0.) {
                        if ((process == ENT_PROCESS_DIS_CC) ||
                            (process == ENT_PROCESS_DIS_CC_OTHER)) {
                                if ((rc = proget_compute(projectile,
                                    ENT_PID_NEUTRON, ENT_PROCESS_DIS_CC_OTHER,
                                    &proget)) != ENT_RETURN_SUCCESS)
                                        ENT_RETURN(rc);
                                *cross_section += N *
                                    cross_section_compute(
                                        mode, proget, cs0, cs1, p1, p2);
                        }
                        if ((process == ENT_PROCESS_DIS_CC) ||
                            (process == ENT_PROCESS_DIS_CC_TOP)) {
                                if ((rc = proget_compute(projectile,
                                    ENT_PID_NEUTRON, ENT_PROCESS_DIS_CC_TOP,
                                    &proget)) != ENT_RETURN_SUCCESS)
                                        ENT_RETURN(rc);
                                *cross_section += N *
                                    cross_section_compute(
                                        mode, proget, cs0, cs1, p1, p2);
                        }
                }
        } else if (process == ENT_PROCESS_DIS_NC) {
                int proget;
                if (Z > 0.) {
                        if ((rc = proget_compute(projectile, ENT_PID_PROTON,
                                 process, &proget)) != ENT_RETURN_SUCCESS)
                                ENT_RETURN(rc);
                        *cross_section += Z *
                            cross_section_compute(
                                mode, proget, cs0, cs1, p1, p2);
                }
                const double N = A - Z;
                if (N > 0.) {
                        if ((rc = proget_compute(projectile, ENT_PID_NEUTRON,
                                 process, &proget)) != ENT_RETURN_SUCCESS)
                                ENT_RETURN(rc);
                        *cross_section += N *
                            cross_section_compute(
                                mode, proget, cs0, cs1, p1, p2);
                }
        } else {
                int proget;
                if ((rc = proget_compute(projectile, ENT_PID_ELECTRON, process,
                         &proget)) != ENT_RETURN_SUCCESS)
                        ENT_RETURN(rc);
                if (proget < PROGET_N - 1) {
                        *cross_section = Z *
                            cross_section_compute(
                                mode, proget, cs0, cs1, p1, p2);
                } else {
                        *cross_section = Z *
                            cross_section_compute(
                                mode, PROGET_N - 2, cs0, cs1, p1, p2) *
                            (1. - ENT_BR_W_TO_ELECTRON - ENT_BR_W_TO_MUON -
                             ENT_BR_W_TO_TAU) / ENT_BR_W_TO_MUON;
                }
        }

        return ENT_RETURN_SUCCESS;
}

/* Generic API function for computing DCSs. */
enum ent_return ent_physics_dcs(struct ent_physics * physics,
    enum ent_pid projectile, double energy, double Z, double A,
    enum ent_process process, double x, double y, double * dcs)
{
        ENT_ACKNOWLEDGE(ent_physics_dcs);

        /* Check the inputs. */
        *dcs = 0.;
        if ((x > 1.) || (x < 0.) || (y > 1.) || (y < 0.))
                ENT_RETURN(ENT_RETURN_DOMAIN_ERROR);

        /* Compute the corresponding DCS. */
        const double d =
            dcs_compute(physics, projectile, energy, Z, A, process, x, y);
        if (d == -DBL_MAX) ENT_RETURN(ENT_RETURN_DOMAIN_ERROR);
        *dcs = d;

        return ENT_RETURN_SUCCESS;
}

/* API interface to PDFs. */
enum ent_return ent_physics_pdf(struct ent_physics * physics,
    enum ent_parton parton, double x, double Q2, double * value)
{
        ENT_ACKNOWLEDGE(ent_physics_pdf);

        /* Check the inputs. */
        if ((x > 1.) || (x <= 0.) || (Q2 < 0.))
                ENT_RETURN(ENT_RETURN_DOMAIN_ERROR);

        if (physics->pdf == NULL) {
                /* No PDF data. Return 0. */
                const int nf = (ENT_N_PARTONS - 1) / 2;
                if ((parton >= -nf) && (parton <= nf)) {
                        *value = 0.;
                } else if (parton == ENT_N_PARTONS) {
                        int i;
                        for (i = 0; i < ENT_N_PARTONS; i++) {
                                value[i] = 0.;
                        }
                } else {
                        ENT_RETURN(ENT_RETURN_DOMAIN_ERROR);
                }
        } else {
                /* Compute the PDF and return. */
                const struct grid * const pdf = physics->pdf;
                float xfx[pdf->nf];
                grid_compute(pdf, x, Q2, xfx);

                const int nf = (ENT_N_PARTONS - 1) / 2;
                const int pdf_nf = (pdf->nf - 1) / 2;
                if ((parton >= -nf) && (parton <= nf)) {
                        /* Return only the requested PDF */
                        if (abs(parton) <= pdf_nf) {
                                *value = xfx[pdf_nf + parton] / x;
                        } else {
                                *value = 0.;
                        }
                } else if (parton == ENT_N_PARTONS) {
                        /* Return all PDFs */
                        int i;
                        for (i = 0; i < ENT_N_PARTONS; i++) {
                                const int ii = pdf_nf + i - nf;
                                value[i] =
                                    ((ii >= 0) && (ii < 2 * pdf_nf + 1)) ?
                                    xfx[ii] : 0.;
                        }
                } else {
                        ENT_RETURN(ENT_RETURN_DOMAIN_ERROR);
                }
        }

        return ENT_RETURN_SUCCESS;
}

/* API interface to DIS SFs. */
enum ent_return ent_physics_sf(struct ent_physics * physics,
    enum ent_pid projectile, enum ent_pid target, enum ent_process process,
    double x, double Q2, double * F2, double * F3, double * FL)
{
        ENT_ACKNOWLEDGE(ent_physics_sf);

        /* Check the inputs. */
        if ((x > 1.) || (x <= 0.) || (Q2 < 0.))
                ENT_RETURN(ENT_RETURN_DOMAIN_ERROR);

        double Z;
        if (target == ENT_PID_PROTON) {
                Z = 1.;
        } else if (target == ENT_PID_NEUTRON) {
                Z = 0.;
        } else {
                ENT_RETURN(ENT_RETURN_DOMAIN_ERROR);
        }

        const int aid = abs(projectile);
        if ((aid != 12) && (aid != 14) && (aid != 16)) {
                ENT_RETURN(ENT_RETURN_DOMAIN_ERROR);
        }

        /* Compute the SFs and return. */
        const double A = 1.;
        double sf[3];
        if ((process == ENT_PROCESS_DIS_NC) ||
            (process == ENT_PROCESS_DIS_CC_OTHER) ||
            (process == ENT_PROCESS_DIS_CC_TOP)) {
                dis_compute_sf(physics, projectile, Z, A, process, x, Q2, sf);
        } else if (process == ENT_PROCESS_DIS_CC) {
                double tmp[3];
                int i;

                dis_compute_sf(physics, projectile, Z, A,
                    ENT_PROCESS_DIS_CC_OTHER, x, Q2, sf);
                dis_compute_sf(physics, projectile, Z, A,
                    ENT_PROCESS_DIS_CC_TOP, x, Q2, tmp);

                for (i = 0; i < 3; i++) {
                        sf[i] += tmp[i];
                }
        } else {
                ENT_RETURN(ENT_RETURN_DOMAIN_ERROR);
        }

        if (F2 != NULL) *F2 = sf[0];
        if (F3 != NULL) *F3 = sf[1];
        if (FL != NULL) *FL = sf[2];

        return ENT_RETURN_SUCCESS;
}

/* Do a transport step and/or update/initialise the local properties of the
 * propagation medium.
 */
static enum ent_return transport_step(struct ent_context * context,
    struct ent_state * state, struct ent_medium ** medium, double * step,
    double * density, double grammage_max, enum ent_event * event)
{
/* Stepping resolution, in m, for locating the interface between two media. */
#define STEP_MIN 1E-06

        *event = ENT_EVENT_NONE;
        const double sgn = (context->ancestor == NULL) ? 1. : -1.;
        if (*step > 0.) {
                /* Apply any distance limit. */
                if ((context->distance_max > 0.) &&
                    (state->distance + *step > context->distance_max)) {
                        *step = context->distance_max - state->distance;
                        if (*step < 0.) *step = 0.;
                        *event = ENT_EVENT_LIMIT_DISTANCE;
                }

                /* Start with a tentative step, if not an initialisation. */
                if (*step) {
                        state->position[0] +=
                            sgn * state->direction[0] * (*step);
                        state->position[1] +=
                            sgn * state->direction[1] * (*step);
                        state->position[2] +=
                            sgn * state->direction[2] * (*step);
                }
        }

        /* Check the new medium. */
        struct ent_medium * medium1;
        double step1 = context->medium(context, state, &medium1);
        if (step1 > 0.) step1 += 0.5 * STEP_MIN;
        double ds_boundary = 0.;
        if (*medium == NULL) {
                /* initialisation step. */
                *medium = medium1;
                if (medium1 == NULL) {
                        *event = ENT_EVENT_EXIT;
                        return ENT_RETURN_SUCCESS;
                }
        } else if ((*medium != NULL) && (medium1 != *medium)) {
                /* A change of medium occured. First let us check for an exact
                 * boundary.
                 */
                double pi[3];
                memcpy(pi, state->position, sizeof(pi));
                double s1 = 0., s2 = -STEP_MIN;
                state->position[0] = pi[0] + s2 * sgn * state->direction[0];
                state->position[1] = pi[1] + s2 * sgn * state->direction[1];
                state->position[2] = pi[2] + s2 * sgn * state->direction[2];
                struct ent_medium * tmp_medium;
                const double tmp_step =
                    context->medium(context, state, &tmp_medium);
                if (tmp_medium != *medium) {
                        /* Let us locate the interface by dichotomy. */
                        step1 = tmp_step;
                        if (step1 > 0.) step1 += 0.5 * STEP_MIN;
                        if (tmp_medium != medium1) medium1 = tmp_medium;
                        s1 = s2;
                        s2 = -(*step);
                        while (fabs(s1 - s2) > STEP_MIN) {
                                double s3 = 0.5 * (s1 + s2);
                                state->position[0] =
                                    pi[0] + s3 * sgn * state->direction[0];
                                state->position[1] =
                                    pi[1] + s3 * sgn * state->direction[1];
                                state->position[2] =
                                    pi[2] + s3 * sgn * state->direction[2];
                                const double tmp_step = context->medium(
                                    context, state, &tmp_medium);
                                if (tmp_medium == *medium) {
                                        s2 = s3;
                                } else {
                                        s1 = s3;
                                        step1 = tmp_step;
                                        if (step1 > 0.) step1 += 0.5 * STEP_MIN;
                                        if (tmp_medium != medium1)
                                                medium1 = tmp_medium;
                                }
                        }
                        state->position[0] =
                            pi[0] + s2 * sgn * state->direction[0];
                        state->position[1] =
                            pi[1] + s2 * sgn * state->direction[1];
                        state->position[2] =
                            pi[2] + s2 * sgn * state->direction[2];
                        *step += s1;
                }
                ds_boundary = s1 - s2;
        }

        /* Get the end step local properties. */
        double density1, step2;
        step2 = (*medium)->density(*medium, state, &density1);
        if ((density1 <= 0.) || ((*medium)->A <= 0.))
                return ENT_RETURN_DOMAIN_ERROR;

        /* Offset the end step position if a boundary crossing occured. */
        if (ds_boundary != 0.) {
                state->position[0] += ds_boundary * sgn * state->direction[0];
                state->position[1] += ds_boundary * sgn * state->direction[1];
                state->position[2] += ds_boundary * sgn * state->direction[2];

                /* Reset any distance limit if no more valid. */
                if ((context->distance_max > 0.) &&
                    (state->distance + *step < context->distance_max))
                        *event = ENT_EVENT_NONE;
        }

        /* Check for a grammage limit. */
        if (*step > 0.) {
                double dX = 0.5 * (*density + density1) * (*step);
                if ((grammage_max > 0.) &&
                    (state->grammage + dX > grammage_max)) {
                        /* Roll back the position and update. */
                        const double dX = grammage_max - state->grammage;
                        double ds;
                        const double drho = density1 - *density;
                        if (fabs(drho) < 1E-03) {
                                ds = 2. * dX / (density1 + *density) - *step;
                        } else {
                                ds = *step *
                                    ((sqrt(density1 * density1 +
                                          2. * dX * drho / *step) -
                                         density1) /
                                            drho -
                                        1.);
                        }
                        state->grammage = grammage_max;
                        state->position[0] += ds * sgn * state->direction[0];
                        state->position[1] += ds * sgn * state->direction[1];
                        state->position[2] += ds * sgn * state->direction[2];
                        *step += ds;
                        *event = ENT_EVENT_LIMIT_GRAMMAGE;
                } else {
                        state->grammage += dX;
                }
                state->distance += *step;
        }

        /* Update the end step density if a valid change of medium occured. */
        if ((ds_boundary != 0.) && (*event == ENT_EVENT_NONE)) {
                *medium = medium1;
                if (medium1 != NULL) {
                        step2 = medium1->density(medium1, state, &density1);
                        if ((density1 <= 0.) || (medium1->A <= 0.))
                                return ENT_RETURN_DOMAIN_ERROR;
                }
        }

        /* Update and return. */
        if ((step2 > 0.) && ((step1 <= 0.) || (step2 < step1))) step1 = step2;
        *step = (step1 < 0) ? 0. : step1;
        *density = density1;
        if (*medium == NULL)
                *event = ENT_EVENT_EXIT;
        return ENT_RETURN_SUCCESS;

#undef STEP_MIN
}

/* Do a straight transport step in a uniform medium of infinite extension. */
static enum ent_event transport_straight(struct ent_context * context,
    struct ent_state * state, double density, double grammage_max)
{
        enum ent_event event;
        double ds = (grammage_max - state->grammage) / density;
        if ((context->distance_max > 0.) &&
            (state->distance + ds >= context->distance_max)) {
                ds = context->distance_max - state->distance;
                state->distance = context->distance_max;
                state->grammage += ds * density;
                event = ENT_EVENT_LIMIT_DISTANCE;
        } else {
                state->distance += ds;
                state->grammage = grammage_max;
                event = ENT_EVENT_LIMIT_GRAMMAGE;
        }

        /* Update the position. */
        const double sgn = (context->ancestor == NULL) ? 1. : -1.;
        state->position[0] += sgn * state->direction[0] * ds;
        state->position[1] += sgn * state->direction[1] * ds;
        state->position[2] += sgn * state->direction[2] * ds;

        return event;
}

/* Sample the y and Q2 parameters in a DIS event. */
static enum ent_return transport_sample_yQ2(struct ent_physics * physics,
    struct ent_context * context, const struct ent_state * neutrino, int proget,
    double * y_, double * Q2)
{
        /* Get back the projectile, the target and the process. */
        const double A = 1.;
        double Z;
        enum ent_process process;
        if (proget < 8) {
                Z = (proget / 4);
                process = (proget % 2) ?
                    ENT_PROCESS_DIS_NC : ENT_PROCESS_DIS_CC_OTHER;
        } else {
                Z = (proget - 8) / 2;
                process = ENT_PROCESS_DIS_CC_TOP;
        }

        /* Check the total cross-section and validate DCS support. */
        double cs1, ymin, ymax, dly;
        int valid;
        for (valid = 0; valid < 1; valid++) {
                double *csl, *csh;
                double pl, ph;
                const int mode = cross_section_prepare(
                     physics->cs_k, neutrino->energy, &csl, &csh, &pl, &ph);
                cs1 = cross_section_compute(mode, proget, csl, csh, pl, ph);
                if (cs1 <= 0.) break;

                dis_get_support_y(
                    physics, proget, neutrino->energy, 1, &ymin, &ymax);
                dly = log(ymax / ymin) / (DIS_Y_N - 1);

                if (ymin >= ymax) break;

                const double ym = 0.5 * (ymin + ymax);
                double xmin, xmax;
                dis_get_support_x(physics, proget, neutrino->energy, ym,
                    &xmin, &xmax, NULL);
                if (xmin >= xmax) break;

                const double xm = 0.5 * (xmin + xmax);
                if (dcs_dis(physics, neutrino->pid, neutrino->energy, Z, A,
                    process, xm, ym) <= 0.) break;
        }
        if (!valid) {
                *y_ = *Q2 = 0.;
                return ENT_RETURN_SUCCESS;
        }

        /* Sample y using inverse CDF. First, let us find a bracketing interval
         * for y.
         */
        const double dle = log(ENERGY_MAX / ENERGY_MIN) / (ENERGY_N - 1);
        double he = log(neutrino->energy / ENERGY_MIN) / dle;
        int ie = (int)he;
        if (ie < 0) {
                ie = 0;
                he = 0.;
        } else if (ie >= ENERGY_N - 1) {
                ie = ENERGY_N - 2;
                he = 1.;
        } else {
                he -= ie;
        }

        double * cdf0 = physics->dis_cdf + (proget * ENERGY_N + ie) * DIS_Y_N;
        double * cdf1 =
            physics->dis_cdf + (proget * ENERGY_N + ie + 1) * DIS_Y_N;
        double * pdf0 = physics->dis_pdf + (proget * ENERGY_N + ie) * DIS_Y_N;
        double * pdf1 =
            physics->dis_pdf + (proget * ENERGY_N + ie + 1) * DIS_Y_N;

        double y, xmin, xmax, ratio;
        for (;;) {
                const double ry = context->random(context);
                int i0 = 0, i1 = DIS_Y_N - 1;
                while (i1 - i0 > 1) {
                        const int i2 = (i0 + i1) / 2;
                        const double cdf2 =
                            cdf0[i2] * (1. - he) + cdf1[i2] * he;
                        if (ry > cdf2) {
                                i0 = i2;
                        } else {
                                i1 = i2;
                        }
                }

                const double y0 = ymin * exp(i0 * dly);
                const double y1 = ymin * exp(i1 * dly);

                /* Secondly, let us interpolate the CDF over the bracketing
                 * interval.
                 */
                const double pdfy0 = pdf0[i0] * (1. - he) + pdf1[i0] * he;
                const double f0 = pdfy0 * y0;
                const double pdfy1 = pdf0[i1] * (1. - he) + pdf1[i1] * he;
                const double f1 = pdfy1 * y1;

                const double r0 = cdf0[i0] * (1. - he) + cdf1[i0] * he;
                const double r1 = cdf0[i1] * (1. - he) + cdf1[i1] * he;
                if (r1 <= r0) continue;
                const double dI = (ry - r0) / (r1 - r0) * 0.5 * (f0 + f1);

                const double ay = 0.5 * (f1 - f0);
                const double by = f0;
                const double cy = -dI;
                double dy = by * by - 4 * ay * cy;
                dy = (dy < 0.) ? 0. : sqrt(dy);
                double hy = (dy - by) / (2 * ay);
                if (hy < 0.) hy = 0.;
                else if (hy > 1.) hy = 1.;

                y = ymin * exp(dly * (i0 + hy));

                /* Check that the result is consistent with x support. */
                dis_get_support_x(
                    physics, proget, neutrino->energy, y, &xmin, &xmax, &ratio);
                if (xmin >= xmax) continue;

                if (dcs_dis(physics, neutrino->pid, neutrino->energy, Z, A,
                    process, 0.5 * (xmin + xmax), y) <= 0.) continue;

                /* Interpolate the PDF as y. */
                const double pdfy = pdfy0 * (1. - hy) + pdfy1 * hy;
                if (pdfy <= 0.) continue;

                /* Forward the selected y value. */
                *y_ = y;

                /* Sample x using an envelope. */
                const double MX = (process == ENT_PROCESS_DIS_NC) ?
                    ENT_MASS_Z : ENT_MASS_W;
                const double x0 = 0.5 * MX * MX /
                    (neutrino->energy * y * ENT_MASS_NUCLEON);

                const double beta = DIS_X_EXPONENT;
                const double bmin = pow(1. + xmin / x0, 1. - beta);
                const double bmax = pow(1. + xmax / x0, 1. - beta);

                /* Safeguard in case of bad estimate of [xmin, xmax] range. */
                const int mintrials = 100;
                const int maxtrials = 1000;
                int trials = (int)(100 * ratio);
                if (trials < mintrials) trials = mintrials;
                else if (trials > maxtrials) trials = maxtrials;

                /* Rejection sampling of x*/
                for (; trials > 0; trials--) {
                        double rx, x;
                        for (;;) {
                                rx = (bmin - bmax) *
                                    context->random(context) + bmax;
                                x = x0 * (pow(rx, 1. / (1. - beta)) - 1.);
                                if ((x > xmin) || (x < xmax)) break;
                        }
                        const double pdf0 =
                            (beta - 1.) * rx / ((bmin - bmax) * (x0 + x));

                        /* Compute the true PDF. */
                        const double dcs1 = dcs_dis(physics, neutrino->pid,
                            neutrino->energy, Z, A, process, x, y);

                        const double pdf1 = dcs1 / (cs1 * pdfy);

                        if (context->random(context) * ratio <= pdf1 / pdf0) {
                                *Q2 = 2 * ENT_MASS_NUCLEON * neutrino->energy *
                                    x * y;
                                break;
                        }
                }
                if (trials == 0) continue;

                break;
        }

        return ENT_RETURN_SUCCESS;
}

/* Sample the E and Q2 parameters in a backward DIS event. */
static enum ent_return backward_sample_EQ2(struct ent_physics * physics,
    struct ent_context * context, struct ent_state * state, int proget,
    enum ent_pid mother, double * E, double * Q2)
{
        /* Get back the projectile, the target and the process. */
        const double A = 1.;
        double Z;
        enum ent_process process;
        if (proget < 8) {
                Z = (proget / 4);
                process = (proget % 2) ?
                    ENT_PROCESS_DIS_NC : ENT_PROCESS_DIS_CC_OTHER;
        } else {
                Z = (proget - 8) / 2;
                process = ENT_PROCESS_DIS_CC_TOP;
        }
        int pid, mode_y;
        if (mother == ENT_PID_NONE) {
                mode_y = 0;
                if (process == ENT_PROCESS_DIS_NC)
                        pid = state->pid;
                else
                        pid = (state->pid > 0) ?
                            state->pid + 1 : state->pid - 1;
        } else {
                mode_y = -1;
                pid = mother;
        }

        /* Sample y using a bias PDF as 1 / y^alpha over [ymin, ymax], where
         * the support is determined from tabulated value as function of the
         * final energy.
         */
        double ymin, ymax;
        dis_get_support_y(physics, proget, state->energy, mode_y, &ymin, &ymax);
        if (ymin >= ymax) {
                /* This is not expected to occur. It likely implies that
                 * the support tabulation routine failed.
                 */
                return ENT_RETURN_DOMAIN_ERROR;
        }

        const double alpha = 0.5;
        double ry, y;
        const double bmin = pow(ymin, 1. - alpha);
        const double bmax = pow(ymax, 1. - alpha);
        double xmin, xmax, cs1 = 0.;
        int trials;
        const int max_trials = 20;
        for (trials = 0; trials < max_trials;) {
                if (mode_y == -1) trials++; /* Close to the kinematic threshold
                                             * the support resolution is
                                             * numerically delicate. This is
                                             * a safeguard in case former check
                                             * failled.
                                             */

                ry = (bmax - bmin) * context->random(context) + bmin;
                y = pow(ry, 1. / (1. - alpha));
                if ((y <= ymin) || (y >= ymax)) continue;

                /* Check cross-section consistency. */
                *E = (mother == ENT_PID_NONE) ?
                    state->energy / (1. - y) : state->energy / y;
                double *csl, *csh;
                double pl, ph;
                const int mode = cross_section_prepare(
                     physics->cs_k, *E, &csl, &csh, &pl, &ph);
                cs1 = cross_section_compute(mode, proget, csl, csh, pl, ph);
                if (cs1 <= 0.) continue;

                /* Check consistency of x support. */
                dis_get_support_x(physics, proget, *E, y, &xmin, &xmax, NULL);
                if ((xmax <= 0.) || (xmin / xmax - 1. >= -FLT_EPSILON)) {
                        continue;
                }

                const double xm = 0.5 * (xmin + xmax);
                if (dcs_dis(physics, pid, *E, Z, A, process, xm, y) <= 0.)
                        continue;

                break;
        }
        if (trials == max_trials) {
                *E = *Q2 = 0.;
                state->weight = 0.;

                return ENT_RETURN_SUCCESS;
        }

        double pdf0 = (1. - alpha) * ry / (y * (bmax - bmin));

        /* Sample x. */
        const double MX =
            (process == ENT_PROCESS_DIS_NC) ? ENT_MASS_Z : ENT_MASS_W;
        const double x0 = 0.5 * MX * MX / (*E * y * ENT_MASS_NUCLEON);

        double pdf1;
        {
                /* Sample x assuming an asymptotic small x PDF. */
                const double beta = DIS_X_EXPONENT;
                const double bmin = pow(1. + xmin / x0, 1. - beta);
                const double bmax = pow(1. + xmax / x0, 1. - beta);
                double rx, x;
                for (;;) {
                        rx = (bmin - bmax) * context->random(context) + bmax;
                        x = x0 * (pow(rx, 1. / (1. - beta)) - 1.);
                        if ((x > xmin) || (x < xmax)) break;
                }
                pdf0 *= (beta - 1.) * rx / ((bmin - bmax) * (x0 + x));
                *Q2 = 2. * ENT_MASS_NUCLEON * (*E) * x * y;

                /* Compute the true PDF. */
                const double dcs1 =
                    dcs_dis(physics, pid, *E, Z, A, process, x, y);
                pdf1 = dcs1 / cs1;
        }

        /* Check and update the BMC weight. */
        double w = pdf1 / pdf0;
        if (w <= 0.) {
                *E = *Q2 = 0.;
                state->weight = 0.;
        } else {
                if (mother == ENT_PID_NONE) {
                        w /= 1. - y;
                } else {
                        w /= y;
                }

                state->weight *= w;
        }

        return ENT_RETURN_SUCCESS;
}

/* Sample the inelasticity, _y_, for an interaction with an electron. The
 * sampling is done with a combination of inverse and rejection sampling.
 */
static double transport_sample_y(struct ent_context * context, double energy,
    int proget, enum ent_pid * ejectile, enum ent_pid * recoil, double * mu)
{
        const double yZ =
            0.5 * ENT_MASS_Z * ENT_MASS_Z / (ENT_MASS_ELECTRON * energy);
        const double yW =
            0.5 * ENT_MASS_W * ENT_MASS_W / (ENT_MASS_ELECTRON * energy);
        const double Re2 = 4. * ENT_PHYS_SIN_THETA_W_2 * ENT_PHYS_SIN_THETA_W_2;
        const double Le = 2. * ENT_PHYS_SIN_THETA_W_2 - 1.;
        const double Le2 = Le * Le;

        *mu = ENT_MASS_ELECTRON;
        *recoil = ENT_PID_ELECTRON;
        double y;
        if (proget == PROGET_ELASTIC_NU_E) {
                const double Cz = yZ / (1. + yZ);
                const double Cw = yW / (1. + yW);
                const double r1 = Cz * Re2;
                const double r2 = r1 + Cz * Le2;
                const double r3 = r2 + 4. * Cw;
                for (;;) {
                        const double r = context->random(context) * r3;
                        const double u = context->random(context);
                        const double t = context->random(context);
                        if (r <= r2)
                                y = yZ * u / (1. - u + yZ);
                        else
                                y = (1. + yW) * u / (yW + u);
                        if (r <= r1) {
                                if (t > (1. - y) * (1. - y)) continue;
                        } else {
                                const double dZ = Le * yZ / (y + yZ);
                                const double dW = 2. * yW / (1. - y + yW);
                                if (t * (dZ * dZ + dW * dW) >
                                    (dZ + dW) * (dZ + dW))
                                        continue;
                        }
                        break;
                }
        } else if (proget == PROGET_ELASTIC_NU_BAR_E) {
                const double rW = ENT_WIDTH_W / ENT_MASS_W;
                const double di =
                    1. / ((yW - 1.) * (yW - 1.) + rW * rW * yW * yW);
                const double a = (yW + yW) * (yW - 1.) * di;
                const double b = -2. * yW * yW * rW * di;
                const double Cz = yZ / (1. + yZ);
                const double r1 = Cz * (Re2 + Le2);
                const double r2 = r1 + (a * a + b * b) / 3.;
                for (;;) {
                        const double r = context->random(context) * r2;
                        const double u = context->random(context);
                        const double t = context->random(context);
                        if (r <= r1)
                                y = yZ * u / (1. - u + yZ);
                        else
                                y = 1. - pow(u, 1. / 3);
                        const double dZ = yZ / (y + yZ);
                        const double dy2 = (1. - y) * (1. - y);
                        const double t0 =
                            dZ * dZ * (Re2 + Le2) + (a * a + b * b) * dy2;
                        const double t1 = dZ * dZ * Re2 +
                            ((Le * dZ + a) * (Le * dZ + a) + b * b) * dy2;
                        if (t * t0 <= t1) break;
                }

        } else if (proget <= PROGET_ELASTIC_NU_BAR_TAU) {
                const double r1 = ((proget == PROGET_ELASTIC_NU_MU) ||
                                      (proget == PROGET_ELASTIC_NU_TAU)) ?
                    Re2 / (Re2 + Le2) :
                    Le2 / (Re2 + Le2);
                for (;;) {
                        const double u = context->random(context);
                        const double r = context->random(context);
                        y = yZ * u / (1. - u + yZ);
                        if (r <= r1) {
                                if (context->random(context) >
                                    (1. - y) * (1. - y))
                                        continue;
                        }
                        break;
                }
        } else if (proget <= PROGET_INVERSE_NU_TAU_TAU) {
                if (ejectile != NULL) *ejectile = ENT_PID_NU_E;
                if (proget == PROGET_INVERSE_NU_MU_MU) {
                        *recoil = ENT_PID_MUON;
                        *mu = ENT_MASS_MUON;
                } else {
                        *recoil = ENT_PID_TAU;
                        *mu = ENT_MASS_TAU;
                }
                const double x = energy - *mu;
                const double w =
                    0.5 * ENT_MASS_W * ENT_MASS_W / ENT_MASS_ELECTRON;
                const double z = context->random(context);
                const double e = w * z * x / (w + x * (1. - z));
                y = 1. - e / energy;
        } else {
                const double u = context->random(context);
                y = 1. - pow(u, 1. / 3.);
                if (proget == PROGET_INVERSE_NU_BAR_E_MU) {
                        if (ejectile != NULL) *ejectile = ENT_PID_NU_BAR_MU;
                        *recoil = ENT_PID_MUON;
                        *mu = ENT_MASS_MUON;
                } else {
                        if (ejectile != NULL) *ejectile = ENT_PID_NU_BAR_TAU;
                        *recoil = ENT_PID_TAU;
                        *mu = ENT_MASS_TAU;
                }
        }

        return y;
}

/* Sample the initial energy in a backward interaction with an electron. An
 * adjoint like method is used where the inelasticity is sampled at E=E_f
 * instead of E_i.
 */
static enum ent_return backward_sample_E(struct ent_physics * physics,
    struct ent_context * context, struct ent_state * state, int proget,
    enum ent_pid * pid0, enum ent_pid * pid1, double * E0, double * E1,
    double * m1)
{
        /* Sample the inelasticity with the forward model but setting
         * E = E_f.
         */
        double y, mu;
        enum ent_pid recoil;
        for (;;) {
                y = transport_sample_y(
                    context, state->energy, proget, NULL, &recoil, &mu);
                const double Ei = (state->pid == recoil) ?
                    state->energy / y :
                    (state->energy - ENT_MASS_ELECTRON) / (1. - y);
                if ((y >= mu / Ei) && (y <= 1.)) break;
        }

        /* Configure the initial state and the missing product */
        enum ent_process process;
        if (state->pid == recoil) {
                /* The final state is the recoiling charged lepton. */
                *E0 = state->energy / y;
                *E1 = state->energy * (1. / y - 1.) + ENT_MASS_ELECTRON;
                if (proget == PROGET_ELASTIC_NU_E) {
                        process = ENT_PROCESS_ELASTIC;
                        *pid1 = *pid0 = ENT_PID_NU_E;
                } else if (proget == PROGET_ELASTIC_NU_BAR_E) {
                        process = ENT_PROCESS_ELASTIC;
                        *pid1 = *pid0 = ENT_PID_NU_BAR_E;
                } else if (proget == PROGET_ELASTIC_NU_MU) {
                        process = ENT_PROCESS_ELASTIC;
                        *pid1 = *pid0 = ENT_PID_NU_MU;
                } else if (proget == PROGET_ELASTIC_NU_BAR_MU) {
                        process = ENT_PROCESS_ELASTIC;
                        *pid1 = *pid0 = ENT_PID_NU_BAR_MU;
                } else if (proget == PROGET_ELASTIC_NU_TAU) {
                        process = ENT_PROCESS_ELASTIC;
                        *pid1 = *pid0 = ENT_PID_NU_TAU;
                } else if (proget == PROGET_INVERSE_NU_MU_MU) {
                        process = ENT_PROCESS_INVERSE_MUON;
                        *pid0 = ENT_PID_NU_MU;
                        *pid1 = ENT_PID_NU_E;
                } else if (proget == PROGET_INVERSE_NU_TAU_TAU) {
                        process = ENT_PROCESS_INVERSE_TAU;
                        *pid0 = ENT_PID_NU_TAU;
                        *pid1 = ENT_PID_NU_E;
                } else if (proget == PROGET_INVERSE_NU_BAR_E_MU) {
                        process = ENT_PROCESS_INVERSE_MUON;
                        *pid0 = ENT_PID_NU_BAR_E;
                        *pid1 = ENT_PID_NU_BAR_MU;
                } else if (proget == PROGET_INVERSE_NU_BAR_E_TAU) {
                        process = ENT_PROCESS_INVERSE_TAU;
                        *pid0 = ENT_PID_NU_BAR_E;
                        *pid1 = ENT_PID_NU_BAR_TAU;
                } else
                        return ENT_RETURN_DOMAIN_ERROR;
                *m1 = 0.;
        } else {
                /* The final state is the neutrino ejectile. */
                *E0 = (state->energy - ENT_MASS_ELECTRON) / (1. - y);
                *E1 = *E0 * y;
                if ((proget >= PROGET_ELASTIC_NU_E) &&
                    (proget <= PROGET_ELASTIC_NU_BAR_TAU)) {
                        process = ENT_PROCESS_ELASTIC;
                        *pid0 = state->pid;
                } else if (proget == PROGET_INVERSE_NU_MU_MU) {
                        process = ENT_PROCESS_INVERSE_MUON;
                        *pid0 = ENT_PID_NU_MU;
                } else if (proget == PROGET_INVERSE_NU_TAU_TAU) {
                        process = ENT_PROCESS_INVERSE_TAU;
                        *pid0 = ENT_PID_NU_TAU;
                } else if (proget == PROGET_INVERSE_NU_BAR_E_MU) {
                        process = ENT_PROCESS_INVERSE_MUON;
                        *pid0 = ENT_PID_NU_BAR_E;
                } else if (proget == PROGET_INVERSE_NU_BAR_E_TAU) {
                        process = ENT_PROCESS_INVERSE_TAU;
                        *pid0 = ENT_PID_NU_BAR_E;
                } else
                        return ENT_RETURN_DOMAIN_ERROR;
                *pid1 = recoil;
                *m1 = mu;
        }

        /* Compute the PDF at E=E_f. */
        const double dcs0 =
            dcs_compute(physics, *pid0, state->energy, 1., 1., process, 0., y);
        double *csl, *csh;
        double pl, ph;
        int mode = cross_section_prepare(
            physics->cs_k, state->energy, &csl, &csh, &pl, &ph);
        const double cs0 =
            cross_section_compute(mode, proget, csl, csh, pl, ph);
        const double pdf0 = dcs0 / cs0;

        /* Compute the true PDF at E=E_i. */
        const double dcs1 =
            dcs_compute(physics, *pid0, *E0, 1., 1., process, 0., y);
        mode = cross_section_prepare(physics->cs_k, *E0, &csl, &csh, &pl, &ph);
        const double cs1 =
            cross_section_compute(mode, proget, csl, csh, pl, ph);
        const double pdf1 = dcs1 / cs1;

        /* Reweight. */
        state->weight *= pdf1 * *E0 / (pdf0 * state->energy);

        return ENT_RETURN_SUCCESS;
}

/* Rotate the state direction. */
static void transport_rotate(
    double direction[3], double cos_theta, double cos_phi, double sin_phi)
{
        /* Check the numerical sine. */
        const double stsq = 1. - cos_theta * cos_theta;
        if (stsq <= 0.) return;
        const double st = sqrt(stsq);

        /* select the co-vectors for the local basis. */
        double u0x = 0., u0y = 0., u0z = 0.;
        const double a0 = fabs(direction[0]);
        const double a1 = fabs(direction[1]);
        const double a2 = fabs(direction[2]);
        if (a0 > a1) {
                if (a0 > a2) {
                        const double nrm = 1. /
                            sqrt(direction[0] * direction[0] +
                                direction[2] * direction[2]);
                        u0x = -direction[2] * nrm, u0z = direction[0] * nrm;
                } else {
                        const double nrm = 1. /
                            sqrt(direction[1] * direction[1] +
                                direction[2] * direction[2]);
                        u0y = direction[2] * nrm, u0z = -direction[1] * nrm;
                }
        } else {
                if (a1 > a2) {
                        const double nrm = 1. /
                            sqrt(direction[0] * direction[0] +
                                direction[1] * direction[1]);
                        u0x = direction[1] * nrm, u0y = -direction[0] * nrm;
                } else {
                        const double nrm = 1. /
                            sqrt(direction[1] * direction[1] +
                                direction[2] * direction[2]);
                        u0y = direction[2] * nrm, u0z = -direction[1] * nrm;
                }
        }
        const double u1x = u0y * direction[2] - u0z * direction[1];
        const double u1y = u0z * direction[0] - u0x * direction[2];
        const double u1z = u0x * direction[1] - u0y * direction[0];

        /* Apply the rotation. */
        direction[0] =
            cos_theta * direction[0] + st * (cos_phi * u0x + sin_phi * u1x);
        direction[1] =
            cos_theta * direction[1] + st * (cos_phi * u0y + sin_phi * u1y);
        direction[2] =
            cos_theta * direction[2] + st * (cos_phi * u0z + sin_phi * u1z);
}

/* Compute the polar angles for an interaction with an atomic electron. */
static void polar_electron(
    double Ep, double Er, double Ee, double mu, double * ce, double * cr)
{
        const double tmp = mu / Er;
        const double pr = Er * sqrt(1. - tmp * tmp);
        *ce = 1. - 0.5 * (pr + Ep - Ee) * (pr + Ee - Ep) / (Ep * Ee);
        if (*ce < -1.)
                *ce = -1.;
        else if (*ce > 1.)
                *ce = 1.;
        *cr = 1. - 0.5 * (Ee + Ep - pr) * (Ee - pr + Ep) / (Ep * pr);
        if (*cr < -1.)
                *cr = -1.;
        else if (*cr > 1.)
                *cr = 1.;
}

/* Two body decay (in CM frame). */
static void decay_two_body_cm(struct ent_context * context,
    enum ent_pid pid0, enum ent_pid pid1, enum ent_pid pid2,
    const double direction[3], double * M_, double * E1s_,
    double p1s[3])
{
        const double M = rest_mass(pid0);
        const double m1 = rest_mass(pid1);
        const double m2 = rest_mass(pid2);
        const double M2 = M * M;
        const double m12 = m1 * m1;
        const double m22 = m2 * m2;

        const double p = 0.5 * sqrt(M2 * M2 + m12 * m12 + m22 * m22
            -2 * ((m12 + m22) * M2 + m12 * m22)) / M;
        const double p2 = p * p;
        const double E1s = sqrt(p2 + m12);

        double ct = 2. * context->random(context) - 1.;
        if (ct > 1.) {
                ct = 1.;
        } else if (ct < -1.) {
                ct = -1.;
        }
        const double phi = 2 * M_PI * context->random(context);

        memcpy(p1s, direction, 3 * sizeof(*p1s));
        transport_rotate(p1s, ct, cos(phi), sin(phi));
        int i;
        for (i = 0; i < 3; i++) p1s[i] *= p;

        *M_ = M;
        *E1s_ = E1s;
}

/* Two body decay. */
static void decay_two_body(struct ent_context * context,
    const struct ent_state * mother, struct ent_state * daughter1,
    struct ent_state * daughter2)
{
        /* CM decay. */
        double M, E1s, p1s[3];
        decay_two_body_cm(context, mother->pid, daughter1->pid, daughter2->pid,
            mother->direction, &M, &E1s, p1s);

        /* Boost to the Lab frame. */
        const double E0 = mother->energy;
        const double p0 = sqrt(E0 * E0 - M * M);
        const double beta = p0 / E0;
        const double gamma = E0 / M;
        const double b[3] = {beta * mother->direction[0],
            beta * mother->direction[1], beta * mother->direction[2]};
        const double bp1s = b[0] * p1s[0] + b[1] * p1s[1] + b[2] * p1s[2];
        const double gr = gamma * gamma / (gamma + 1.);

        daughter1->energy = gamma * (E1s + bp1s);
        daughter2->energy = E0 - daughter1->energy;

        int i;
        for (i = 0; i < 3; i++) {
                daughter1->direction[i] =
                    p1s[i] + (gr * bp1s + gamma * E1s) * b[i];
                daughter2->direction[i] = p0 * mother->direction[i] -
                    daughter1->direction[i];
        }

        const double nrm1 = 1. / sqrt(
            daughter1->direction[0] * daughter1->direction[0] +
            daughter1->direction[1] * daughter1->direction[1] +
            daughter1->direction[2] * daughter1->direction[2]);
        for (i = 0; i < 3; i++) {
                daughter1->direction[i] *= nrm1;
        }
        memcpy(daughter1->position, mother->position,
            sizeof(daughter1->position));
        daughter1->weight = mother->weight;

        const double nrm2 = 1. / sqrt(
            daughter2->direction[0] * daughter2->direction[0] +
            daughter2->direction[1] * daughter2->direction[1] +
            daughter2->direction[2] * daughter2->direction[2]);
        for (i = 0; i < 3; i++) {
                daughter2->direction[i] *= nrm2;
        }
        memcpy(daughter2->position, mother->position,
            sizeof(daughter2->position));
        daughter2->weight = mother->weight;
}

/* Two body backward decay. */
static void undecay_two_body(struct ent_context * context,
    const struct ent_state * daughter1, struct ent_state * mother,
    struct ent_state * daughter2)
{
        /* CM decay. */
        double M, E1s, p1s[3];
        decay_two_body_cm(context, mother->pid, daughter1->pid, daughter2->pid,
            daughter1->direction, &M, &E1s, p1s);

        /* Compute Lorentz transform parameters. */
        const double m1 = rest_mass(daughter1->pid);
        const double p1 = sqrt(daughter1->energy * daughter1->energy - m1 * m1);
        const double dp[3] = {
                p1 * daughter1->direction[0] - p1s[0],
                p1 * daughter1->direction[1] - p1s[1],
                p1 * daughter1->direction[2] - p1s[2]
        };
        const double dp2 = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];
        const double gamma = 1. + dp2 /
            (daughter1->energy * E1s + p1 * (daughter1->direction[0] * p1s[0] +
             daughter1->direction[1] * p1s[1] +
             daughter1->direction[2] * p1s[2]) + m1 * m1);
        const double b1 = (gamma + 1.) / (gamma * (daughter1->energy + E1s));
        const double b[3] = {b1 * dp[0], b1 * dp[1], b1 * dp[2]};

        /* Apply transform. */
        mother->energy = gamma * M;
        daughter2->energy = mother->energy - daughter1->energy;

        int i;
        const double nrm0 = 1. / sqrt(dp2);
        for (i = 0; i < 3; i++) {
                mother->direction[i] = dp[i] * nrm0;
                daughter2->direction[i] =
                    mother->energy * b[i] - p1 * daughter1->direction[i];
        }
        memcpy(mother->position, daughter1->position,
            sizeof(mother->position));
        mother->weight = daughter1->weight;

        const double nrm2 = 1. / sqrt(
            daughter2->direction[0] * daughter2->direction[0] +
            daughter2->direction[1] * daughter2->direction[1] +
            daughter2->direction[2] * daughter2->direction[2]);
        for (i = 0; i < 3; i++) {
                daughter2->direction[i] *= nrm2;
        }
        memcpy(daughter2->position, mother->position,
            sizeof(daughter2->position));
        daughter2->weight = mother->weight;

        /* Apply BMC weight. */
        const double s0 = mother->energy + M;
        const double s1 = daughter1->energy + E1s;
        const double p0 = sqrt(mother->energy * mother->energy - M * M);
        const double w = s0 * s0 * p1 / (s1 * s1 * p0);

        mother->weight *= w;
        daughter2->weight *= w;
}

/* Decay DIS products, whenever a top is produced. */
static void process_dis_products(const struct ent_physics * physics,
    struct ent_context * context, enum ent_pid neutrino,
    struct ent_state * struck_quark, struct ent_state * products)
{
        /* Check kinematic limit for top production. */
        const double mq = (abs(struck_quark->pid) == ENT_PID_TOP) ?
            physics->mteff[neutrino < 0] : rest_mass(struck_quark->pid);
        if (struck_quark->energy <= mq * mq / (2 * ENT_MASS_NUCLEON)) {
                struck_quark->pid = ENT_PID_NONE;
        }

        /* Check top case with a W leptonic decay. */
        const double rd = (abs(struck_quark->pid) == ENT_PID_TOP) ?
            context->random(context) : 1.;
        if (rd < ENT_BR_W_TO_ELECTRON + ENT_BR_W_TO_MUON +
            ENT_BR_W_TO_TAU) {
                /* Top production with leptonic decay of W.
                 * First let us decay the top quark. */
                struct ent_state w_boson, b_quark;
                if (struck_quark->pid > 0) {
                        w_boson.pid = ENT_PID_W_PLUS;
                        b_quark.pid = ENT_PID_BOTTOM;
                } else {
                        w_boson.pid = ENT_PID_W_MINUS;
                        b_quark.pid = ENT_PID_BOTTOM_BAR;
                }
                decay_two_body(context, struck_quark, &w_boson, &b_quark);

                /* Then, decay the W. */
                struct ent_state l_companion, n_companion;
                if (rd < ENT_BR_W_TO_ELECTRON) {
                        l_companion.pid = ENT_PID_POSITRON;
                        n_companion.pid = ENT_PID_NU_E;
                } else if (rd < ENT_BR_W_TO_ELECTRON +
                    ENT_BR_W_TO_MUON) {
                        l_companion.pid = ENT_PID_MUON_BAR;
                        n_companion.pid = ENT_PID_NU_MU;
                } else {
                        l_companion.pid = ENT_PID_TAU_BAR;
                        n_companion.pid = ENT_PID_NU_TAU;
                }
                if (w_boson.pid < 0) {
                        l_companion.pid = -l_companion.pid;
                        n_companion.pid = -n_companion.pid;
                }
                decay_two_body(context, &w_boson, &l_companion, &n_companion);

                memcpy(products, &n_companion, sizeof n_companion);
                memcpy(products + 1, &l_companion, sizeof l_companion);

                b_quark.pid = ENT_PID_HADRON; /* Hadronize. */
                memcpy(products + 2, &b_quark, sizeof b_quark);
        } else {
                struck_quark->pid = ENT_PID_HADRON; /* Hadronize. */
                memcpy(products, struck_quark, sizeof(*struck_quark));
        }
}

/* Process a forward interaction vertex. */
static enum ent_return transport_vertex_forward(struct ent_physics * physics,
    struct ent_context * context, int proget, struct ent_state * state,
    struct ent_state * products)
{
        /* Process the corresponding vertex. */
        if (proget < PROGET_N_DIS) {
                /* Backup the initial neutrino state. */
                struct ent_state neutrino;
                memcpy(&neutrino, state, sizeof(neutrino));

                double y, Q2;
                enum ent_return rc = transport_sample_yQ2(
                    physics, context, state, proget, &y, &Q2);
                if (rc != ENT_RETURN_SUCCESS) return rc;

                if ((y <= 0.) && (Q2 <= 0.)) {
                        /* Collision did not occur, e.g. because the projectile
                         * is below the kinematic threshold.
                         */
                        if (products != NULL) {
                                memset(products, 0x0, sizeof(*products));
                        }
                        return ENT_RETURN_SUCCESS;
                }

                if ((proget >= 8) || ((proget % 2) == 0)) {
                        /* Charged current event: update the lepton PID. */
                        state->pid += (state->pid > 0) ? -1 : 1;
                }

                const double mu = rest_mass(state->pid);
                if (y < 1.) {
                        /* Compute the product lepton's energy and its
                         * direction.
                         */
                        const double Emu = neutrino.energy * (1. - y);
                        state->energy = (Emu > mu) ? Emu : mu;

                        /* Always approximate out-going lepton as mass-less
                         * when computing the scattering angle.
                         */
                        double
                            ct = 1. - 0.5 * Q2 / (neutrino.energy * Emu);
                        if (ct < 1.) {
                                if (ct < -1.) ct = -1.;
                                const double phi = 2. * M_PI *
                                     context->random(context);
                                const double cp = cos(phi);
                                const double sp = sin(phi);
                                transport_rotate(state->direction, ct, cp, sp);
                        }
                } else {
                        /* This is a total conversion. Let's update the
                         * energy.
                         */
                        state->energy = mu;
                }

                /* Copy back the product data if requested. */
                if (products != NULL) {
                        struct ent_state struck_quark = { 0. };
                        struck_quark.energy = neutrino.energy - state->energy;
                        struck_quark.weight = state->weight;
                        memcpy(struck_quark.position, state->position,
                            sizeof struck_quark.position);
                        if (proget >= 8) {
                                struck_quark.pid = (state->pid > 0) ?
                                    ENT_PID_TOP : ENT_PID_TOP_BAR;
                        } else {
                                struck_quark.pid = ENT_PID_NONE;
                        }

                        /* Compute the struck quark final direction
                         * from momentum conservation. */
                        const double pmu = (state->energy > mu) ?
                            sqrt(state->energy * state->energy - mu * mu) :
                            0.;
                        int i;
                        double d = 0.;
                        for (i = 0; i < 3; i++) {
                                const double tmp =
                                    neutrino.energy * neutrino.direction[i] -
                                    pmu * state->direction[i];
                                struck_quark.direction[i] = tmp;
                                d += tmp * tmp;
                        }
                        if (d > 0.) {
                                d = 1. / sqrt(d);
                                for (i = 0; i < 3; i++) {
                                        struck_quark.direction[i] *= d;
                                }
                        }

                        /* Update energies. */
                        process_dis_products(physics, context, state->pid,
                            &struck_quark, products);
                }
        } else if (proget < PROGET_GLASHOW_HADRONS) {
                /* This is an interaction with an atomic electron and a
                 * neutrino in the final states.
                 */
                enum ent_pid ejectile = state->pid;
                enum ent_pid recoil = ENT_PID_NONE;
                double y, mu, ce, cr;
                for (;;) {
                        /* Let's first sample the inelasticity _y_. */
                        y = transport_sample_y(context, state->energy,
                            proget, &ejectile, &recoil, &mu);
                        if ((y >= mu / state->energy) && (y <= 1.)) break;
                }

                /* Then, compute the cosines of the polar angles. */
                const double Ep = state->energy;
                const double Er = state->energy * y;
                const double Ee = Ep + ENT_MASS_ELECTRON - Er;
                polar_electron(Ep, Er, Ee, mu, &ce, &cr);

                /* Update the particles states. */
                const double phi = 2. * M_PI * context->random(context);
                const double cp = cos(phi);
                const double sp = sin(phi);
                if ((products != NULL) && (recoil != ENT_PID_NONE)) {
                        memcpy(products, state, sizeof(*products));
                        products->pid = recoil;
                        products->energy = Er;
                        transport_rotate(products->direction, cr, -cp, -sp);
                }
                state->pid = ejectile;
                state->energy = Ee;
                transport_rotate(state->direction, ce, cp, sp);
        } else if (proget == PROGET_GLASHOW_HADRONS) {
                /* This is a total conversion of a anti nu_e neutrino on an
                 * atomic electron. */
                state->pid = ENT_PID_HADRON;
        } else {
                return ENT_RETURN_DOMAIN_ERROR;
        }

        return ENT_RETURN_SUCCESS;
}

/* Process a BMC vertex. */
static enum ent_return transport_vertex_backward(struct ent_physics * physics,
    struct ent_context * context, int proget, struct ent_state * state,
    struct ent_state * products)

{
        /* Process the corresponding vertex. */
        if (proget < 0) {
                /* This is a backward decay from a muon or from a tau. It must
                 * be randomised with an external package. Let us flag this
                 * case by modifying the PID.
                 */
                return ENT_RETURN_SUCCESS;
        }
        if (proget < PROGET_N_DIS) {
                /* Backup the initial state. */
                struct ent_state lepton;
                if (products != NULL) {
                        memcpy(&lepton, state, sizeof(lepton));
                }

                /* Sample the energy loss. */
                enum ent_return rc;
                double Enu, Q2;
                const int ntrials = 20;
                int i;
                for (i = 0; i < ntrials; i++) {
                        if ((rc = backward_sample_EQ2(physics, context, state,
                            proget, ENT_PID_NONE, &Enu, &Q2)) ==
                            ENT_RETURN_SUCCESS)
                                break;
                        /* This should not occur, except due to numeric rounding
                         * errors. Thus, whenever it happens, let us generate
                         * a new event.
                         */
                }
                if (rc != ENT_RETURN_SUCCESS) return rc;
                if (state->weight <= 0.) {
                        return ENT_RETURN_SUCCESS;
                }

                /* Update the MC state. */
                if ((proget >= 8) || ((proget % 2) == 0)) {
                        /* Charged current event. */
                        state->pid += (state->pid > 0) ? 1 : -1;
                }
                const double Emu = state->energy;
                state->energy = Enu;
                double ct = 1. - 0.5 * Q2 / (Enu * Emu);
                if (ct < 1.) {
                        if (ct < -1.) ct = -1.;
                        const double phi = 2. * M_PI * context->random(context);
                        const double cp = cos(phi);
                        const double sp = sin(phi);
                        transport_rotate(state->direction, ct, cp, sp);
                }

                /* Process collision side products. */
                if (products != NULL) {
                        struct ent_state struck_quark = { 0. };
                        struck_quark.energy = state->energy - lepton.energy;
                        struck_quark.weight = state->weight;
                        memcpy(struck_quark.position, state->position,
                            sizeof struck_quark.position);
                        if (proget >= 8) {
                                struck_quark.pid = (state->pid > 0) ?
                                    ENT_PID_TOP : ENT_PID_TOP_BAR;
                        } else {
                                struck_quark.pid = ENT_PID_NONE;
                        }

                        /* Compute the struck quark final direction from
                         * momentum conservation.
                         */
                        const double mu = rest_mass(state->pid);
                        const double pmu = (Emu > mu) ?
                            sqrt(Emu * Emu - mu * mu) : 0.;
                        int i;
                        double d = 0.;
                        for (i = 0; i < 3; i++) {
                                const double tmp =
                                    state->direction[i] * state->energy -
                                    pmu * lepton.direction[i];
                                struck_quark.direction[i] = tmp;
                                d += tmp * tmp;
                        }
                        if (d > 0.) {
                                d = 1. / sqrt(d);
                                for (i = 0; i < 3; i++) {
                                        struck_quark.direction[i] *= d;
                                }
                        }

                        /* Process decay products. */
                        process_dis_products(physics, context, state->pid,
                            &struck_quark, products);
                }
        } else if (proget < PROGET_GLASHOW_HADRONS) {
                /* This is an interaction with an atomic electron and a
                 * neutrino in the final states.  Let's first sample the
                 * initial energy E.
                 */
                enum ent_pid pid0, pid1;
                double E0, E1, m1;
                enum ent_return rc;
                if ((rc = backward_sample_E(physics, context, state, proget,
                         &pid0, &pid1, &E0, &E1, &m1)) != ENT_RETURN_SUCCESS)
                        return rc;

                /* Then, compute the cosines of the polar angles. */
                double Er, Ee, mu, c0, c1, *ce, *cr;
                if (m1 > 0.) {
                        /* The final state is the neutrino ejectile. */
                        Ee = state->energy;
                        ce = &c0;
                        Er = E1;
                        mu = m1;
                        cr = &c1;
                } else {
                        /* The final state is the recoiling charged lepton. */
                        Er = state->energy;
                        if (state->pid == ENT_PID_ELECTRON)
                                mu = ENT_MASS_ELECTRON;
                        else if (state->pid == ENT_PID_MUON)
                                mu = ENT_MASS_MUON;
                        else if (state->pid == ENT_PID_TAU)
                                mu = ENT_MASS_TAU;
                        else
                                return ENT_RETURN_DOMAIN_ERROR;
                        cr = &c0;
                        Ee = E1;
                        ce = &c1;
                }
                polar_electron(E0, Er, Ee, mu, ce, cr);

                /* Update the particles states. */
                const double phi = 2. * M_PI * context->random(context);
                const double cp = cos(phi);
                const double sp = sin(phi);
                state->pid = pid0;
                state->energy = E0;
                transport_rotate(state->direction, c0, cp, sp);
                if (products != NULL) {
                        memcpy(products, state, sizeof(*products));
                        products->pid = pid1;
                        products->energy = E1;
                        transport_rotate(products->direction, c1, -cp, -sp);
                }
        } else
                return ENT_RETURN_DOMAIN_ERROR;

        return ENT_RETURN_SUCCESS;
}

/* Process a BMC vertex from a top decay. */
static enum ent_return transport_vertex_backward_top_decay(
    struct ent_physics * physics, struct ent_context * context, int proget,
    struct ent_state * state, struct ent_state * products)
{
        /* Backup the lepton companion state. */
        struct ent_state l_companion;
        memcpy(&l_companion, state, sizeof l_companion);

        /* Undecay the daughter particle. First, undecay to a W boson. */
        struct ent_state w_boson = { 0. }, n_companion = { 0. };
        w_boson.pid = ENT_PID_W_PLUS;
        const int aid = abs(state->pid);
        if (aid % 2) {
                /* Charged lepton case. */
                n_companion.pid = aid + 1;
                if (state->pid > 0) {
                        w_boson.pid = -w_boson.pid;
                        n_companion.pid = -n_companion.pid;
                }
        } else {
                /* Neutrino case. */
                n_companion.pid = -aid + 1;
                if (state->pid < 0) {
                        w_boson.pid = -w_boson.pid;
                        n_companion.pid = -n_companion.pid;
                }
        }
        undecay_two_body(context, &l_companion, &w_boson, &n_companion);

        /* Secondly, undecay the W boson to a top quark. */
        struct ent_state t_quark = { 0. }, b_quark = { 0. };
        if (w_boson.pid > 0) {
                t_quark.pid = ENT_PID_TOP;
                b_quark.pid = ENT_PID_BOTTOM;
        } else {
                t_quark.pid = ENT_PID_TOP_BAR;
                b_quark.pid = ENT_PID_BOTTOM_BAR;
        }
        undecay_two_body(context, &w_boson, &t_quark, &b_quark);

        /* Check the kinematic limit. */
        const double mt = physics->mteff[proget - 8];
        if (t_quark.energy <= mt * mt / (2 * ENT_MASS_NUCLEON)) {
                state->weight = 0.;
                return ENT_RETURN_SUCCESS;
        }

        /* Sample the primary neutrino flavour from the hadron product.
         * XXX any weight?
         */
        double p[3];
        enum ent_pid pid[3];
        int i;
        for (i = 0; i < 3; i++) {
                if (state->pid % 2) {
                        /* Charged lepton product is opposite CP. */
                        pid[i] = (state->pid > 0) ? -12 - 2 * i : 12 + 2 * i;
                } else {
                        /* Neutrino product is same CP. */
                        pid[i] = (state->pid > 0) ? 12 + 2 * i : -12 - 2 * i;
                }
                const double d = context->ancestor(context, pid[i], state);
                p[i] = ((i > 0) ? p[i - 1] : 0.) + ((d > 0.) ? d : 0.);
        }
        const double r = context->random(context) * p[2];

        enum ent_pid ancestor;
        if (r < p[0]) ancestor = pid[0];
        else if (r < p[1]) ancestor = pid[1];
        else ancestor = pid[2];

        /* Backward sample the energy loss from the t-quark. */
        enum ent_return rc;
        double Enu, Q2;
        const int ntrials = 20;
        for (i = 0; i < ntrials; i++) {
                if ((rc = backward_sample_EQ2(physics, context, &t_quark,
                         proget, ancestor, &Enu, &Q2)) == ENT_RETURN_SUCCESS)
                        break;
                /* This should not occur, except due to numeric rounding
                 * errors. Thus, whenever it happens, let us generate
                 * a new event.
                 */
        }
        if (rc != ENT_RETURN_SUCCESS) return rc;
        if (t_quark.weight <= 0.) {
                state->weight = 0.;
                return ENT_RETURN_SUCCESS;
        }

        /* Set the neutrino state. */
        state->pid = ancestor;
        state->energy = Enu;
        state->weight = t_quark.weight;
        memcpy(state->direction, &t_quark.direction, sizeof state->direction);

        /* Randomise the neutrino direction. */
        const double Emu = Enu - t_quark.energy;
        const double ph = sqrt((Enu - Emu) * (Enu - Emu) + Q2);
        double ct = 1. - 0.5 * Q2 / (Enu * Emu);
        ct = (Enu - Emu * ct) / ph; /* lepton to hadron angle. */
        if (ct < 1.) {
                if (ct < -1.) ct = -1.;
                const double phi = 2. * M_PI * context->random(context);
                const double cp = cos(phi);
                const double sp = sin(phi);
                transport_rotate(state->direction, ct, cp, sp);
        }

        /* Copy back the products. */
        if (products != NULL) {
                memset(products, 0x0, sizeof *products);
                products->pid = (ancestor > 0) ? ancestor - 1 : ancestor + 1;
                const double mu = rest_mass(products->pid);
                products->energy = (Emu > mu) ? Emu : mu;
                products->weight = state->weight;
                memcpy(products->position, state->position,
                    sizeof(products->position));

                /* Set the primary lepton direction from momentum conservation.
                 */
                for (i = 0; i < 3; i++) {
                        products->direction[i] = Enu * state->direction[i] -
                            t_quark.energy * t_quark.direction[i];
                }
                const double nrm = 1. / sqrt(
                    products->direction[0] * products->direction[0] +
                    products->direction[1] * products->direction[1] +
                    products->direction[2] * products->direction[2]);
                for (i = 0; i < 3; i++) {
                        products->direction[i] *= nrm;
                }

                n_companion.weight = state->weight;
                memcpy(products + 1, &n_companion, sizeof n_companion);

                b_quark.pid = ENT_PID_HADRON; /* Hadronization ... */
                b_quark.weight = state->weight;
                memcpy(products + 2, &b_quark, sizeof b_quark);
        }

        return ENT_RETURN_SUCCESS;
}

/* Process an interaction vertex. */
static enum ent_return transport_vertex(struct ent_physics * physics,
    struct ent_context * context, int proget, struct ent_state * state,
    struct ent_state * products)
{
        if (context->ancestor == NULL)
                return transport_vertex_forward(
                    physics, context, proget, state, products);
        else
                return transport_vertex_backward(
                    physics, context, proget, state, products);
}

/* Compute the tranport cross-sections for a given projectile and
 * medium. */
static enum ent_return transport_cross_section(struct ent_physics * physics,
    enum ent_pid projectile, double energy, double Z, double A,
    enum ent_process process, double * cs)
{
        /* Build the interpolation or extrapolation factors. */
        double *cs0, *cs1;
        double p1, p2;
        const int mode = cross_section_prepare(
            physics->cs_t, energy, &cs0, &cs1, &p1, &p2);

        /* Build the table of cumulative cross-section values. */
        double N0, N1, Z0, Z1;
        if (projectile > 0) {
                Z0 = Z;
                Z1 = 0.;
                N0 = A - Z;
                N1 = 0.;
        } else {
                Z0 = 0.;
                Z1 = Z;
                N0 = 0.;
                N1 = A - Z;
        }

        memset(cs, 0x0, PROGET_N * sizeof(*cs));

        if ((process == ENT_PROCESS_NONE) ||
            (process == ENT_PROCESS_DIS_CC) ||
            (process == ENT_PROCESS_DIS_CC_OTHER)) {
                /* Charged current DIS of a neutrino on a neutron. */
                cs[0] = ((N0 > 0.) ?
                    N0 * cross_section_compute(mode, 0, cs0, cs1, p1, p2) : 0.);

                /* Charged current DIS of an anti-neutrino on a neutron. */
                cs[2] = ((N1 > 0.) ?
                    N1 * cross_section_compute(mode, 2, cs0, cs1, p1, p2) : 0.);

                /* Charged current DIS of a neutrino on a proton. */
                cs[4] = ((Z0 > 0.) ?
                    Z0 * cross_section_compute(mode, 4, cs0, cs1, p1, p2) : 0.);

                /* Charged current DIS of an anti-neutrino on a proton. */
                cs[6] = ((Z1 > 0.) ?
                    Z1 * cross_section_compute(mode, 6, cs0, cs1, p1, p2) : 0.);
        }


        if ((process == ENT_PROCESS_NONE) ||
            (process == ENT_PROCESS_DIS_NC)) {
                /* Neutral current DIS of a neutrino on a neutron. */
                cs[1] = ((N0 > 0.) ?
                    N0 * cross_section_compute(mode, 1, cs0, cs1, p1, p2) : 0.);

                /* Neutral current DIS of an anti-neutrino on a neutron. */
                cs[3] = ((N1 > 0.) ?
                    N1 * cross_section_compute(mode, 3, cs0, cs1, p1, p2) : 0.);

                /* Neutral current DIS of a neutrino on a proton. */
                cs[5] = ((Z0 > 0.) ?
                    Z0 * cross_section_compute(mode, 5, cs0, cs1, p1, p2) : 0.);

                /* Neutral current DIS of an anti-neutrino on a proton. */
                cs[7] = ((Z1 > 0.) ?
                    Z1 * cross_section_compute(mode, 7, cs0, cs1, p1, p2) : 0.);
        }


        if ((process == ENT_PROCESS_NONE) ||
            (process == ENT_PROCESS_DIS_CC) ||
            (process == ENT_PROCESS_DIS_CC_TOP)) {
                /* Top production for a CC neutrino on a neutron. */
                cs[8] = ((N0 > 0.) ?
                    N0 * cross_section_compute(mode, 8, cs0, cs1, p1, p2) :
                    0.);

                /* Top production for a CC anti-neutrino on a neutron. */
                cs[9] = ((N1 > 0.) ?
                    N1 * cross_section_compute(mode, 9, cs0, cs1, p1, p2) :
                    0.);

                /* Top production for a CC neutrino on a proton. */
                cs[10] = ((Z0 > 0.) ?
                    Z0 * cross_section_compute(mode, 10, cs0, cs1, p1, p2) :
                    0.);

                /* Top production for a CC anti-neutrino on a proton. */
                cs[11] = ((Z1 > 0.) ?
                    Z1 * cross_section_compute(mode, 11, cs0, cs1, p1, p2) :
                    0.);
        }

        if ((process == ENT_PROCESS_NONE) ||
            (process == ENT_PROCESS_ELASTIC)) {
                if (projectile == ENT_PID_NU_E) {
                        /* Elastic scattering of a nu_e on a an electron. */
                        cs[12] = ((Z > 0.) ?
                            Z * cross_section_compute(
                                mode, 12, cs0, cs1, p1, p2) : 0.);
                } else if (projectile == ENT_PID_NU_BAR_E) {
                        /* Elastic scattering of an anti nu_e on a an
                         * electron.
                         */
                        cs[13] = ((Z > 0.) ?
                            Z * cross_section_compute(
                                mode, 13, cs0, cs1, p1, p2) : 0.);
                } else if (projectile == ENT_PID_NU_MU) {
                        /* Elastic scattering of a nu_mu on a an electron. */
                        cs[14] = ((Z > 0.) ?
                            Z * cross_section_compute(
                                mode, 14, cs0, cs1, p1, p2) : 0.);
                } else if (projectile == ENT_PID_NU_BAR_MU) {
                        /* Elastic scattering of an anti nu_mu on a an
                         * electron.
                         */
                        cs[15] = ((Z > 0.) ?
                            Z * cross_section_compute(
                                mode, 15, cs0, cs1, p1, p2) : 0.);
                } else if (projectile == ENT_PID_NU_TAU) {
                        /* Elastic scattering of a nu_tau on a an electron. */
                        cs[16] = ((Z > 0.) ?
                            Z * cross_section_compute(
                                mode, 16, cs0, cs1, p1, p2) : 0.);
                } else if (projectile == ENT_PID_NU_BAR_TAU) {
                        /* Elastic scattering of an anti nu_tau on a an
                         * electron.
                         */
                        cs[17] = ((Z > 0.) ?
                            Z * cross_section_compute(
                                mode, 17, cs0, cs1, p1, p2) : 0.);
                }
        }

        if ((process == ENT_PROCESS_NONE) ||
            (process == ENT_PROCESS_INVERSE_MUON)) {
                /* Inverse muon decay with a nu_mu projectile. */
                if (projectile == ENT_PID_NU_MU) {
                        cs[18] = ((Z > 0.) ?
                            Z * cross_section_compute(
                                mode, 18, cs0, cs1, p1, p2) : 0.);
                } else if (projectile == ENT_PID_NU_BAR_E) {
                        cs[20] = ((Z > 0.) ?
                            Z * cross_section_compute(
                                mode, 20, cs0, cs1, p1, p2) : 0.);
                }
        }

        if ((process == ENT_PROCESS_NONE) ||
            (process == ENT_PROCESS_INVERSE_TAU)) {
                /* Inverse tau decay with a nu_tau projectile. */
                if (projectile == ENT_PID_NU_TAU) {
                        cs[19] = ((Z > 0.) ?
                            Z * cross_section_compute(
                                mode, 19, cs0, cs1, p1, p2) : 0.);
                } else if (projectile == ENT_PID_NU_BAR_E) {
                        cs[21] = ((Z > 0.) ?
                            Z * cross_section_compute(
                                mode, 21, cs0, cs1, p1, p2) : 0.);
                }
        }

        if ((process == ENT_PROCESS_NONE) ||
            (process == ENT_PROCESS_GLASHOW_HADRON)) {
                if (projectile == ENT_PID_NU_BAR_E) {
                        /* Anti nu_e projectile on an electron with hadrons
                         * production.
                         */
                        const double d = ((Z > 0.) ?
                            Z * cross_section_compute(
                                mode, 20, cs0, cs1, p1, p2) : 0.);

                        cs[22] = d * (1. - ENT_BR_W_TO_ELECTRON -
                             ENT_BR_W_TO_MUON - ENT_BR_W_TO_TAU) /
                             ENT_BR_W_TO_MUON;
                }
        }

        /* Build the cumulative sum. */
        int i;
        for (i = 1; i < PROGET_N; i++) {
                cs[i] += cs[i - 1];
        }

        /* Check the result and return. */
        if (cs[PROGET_N - 1] <= 0.) return ENT_RETURN_DOMAIN_ERROR;
        return ENT_RETURN_SUCCESS;
}

/* Generic builder for the ancestor likeliness. */
static void ancestor_likeliness_fill(struct ent_context * context,
    struct ent_state * daughter, double Z, int * np, int * proget_v, double * p,
    enum ent_pid * ancestor_v, int mode, double * cs0, double * cs1, double p1,
    double p2, const int * pid, const int * proget, int n)
{
        int i;
        for (i = 0; i < n; i++) {
                const double rho0 =
                    context->ancestor(context, pid[i], daughter);
                if (rho0 > 0.) {
                        proget_v[*np] = proget[i];
                        const double p0 = (*np > 0) ? p[*np - 1] : 0.;
                        p[*np] = p0 +
                            rho0 * Z *
                                cross_section_compute(
                                    mode, proget[i], cs0, cs1, p1, p2);
                        ancestor_v[(*np)++] = pid[i];
                }
        }
}

/* Build the ancestor likeliness for an elastic processes with an electron
 * final state.
 */
static void ancestor_electron_elastic(struct ent_context * context,
    struct ent_state * daughter, double Z, int * np, int * proget_v, double * p,
    enum ent_pid * ancestor_v, int mode, double * cs0, double * cs1, double p1,
    double p2)
{
        const int pid[6] = { ENT_PID_NU_E, ENT_PID_NU_BAR_E, ENT_PID_NU_MU,
                ENT_PID_NU_BAR_MU, ENT_PID_NU_TAU, ENT_PID_NU_BAR_TAU };
        const int proget[6] = { PROGET_ELASTIC_NU_E, PROGET_ELASTIC_NU_BAR_E,
                PROGET_ELASTIC_NU_MU, PROGET_ELASTIC_NU_BAR_MU,
                PROGET_ELASTIC_NU_TAU, PROGET_ELASTIC_NU_BAR_TAU };

        ancestor_likeliness_fill(context, daughter, Z, np, proget_v, p,
            ancestor_v, mode, cs0, cs1, p1, p2, pid, proget, 6);
}

/* Build the ancestor likeliness for an inverse muon decay process with a
 * muon final state.
 */
static void ancestor_muon_inverse(struct ent_context * context,
    struct ent_state * daughter, double Z, int * np, int * proget_v, double * p,
    enum ent_pid * ancestor_v, int mode, double * cs0, double * cs1, double p1,
    double p2)
{
        const int pid[2] = { ENT_PID_NU_MU, ENT_PID_NU_BAR_E };
        const int proget[2] = { PROGET_INVERSE_NU_MU_MU,
                PROGET_INVERSE_NU_BAR_E_MU };

        ancestor_likeliness_fill(context, daughter, Z, np, proget_v, p,
            ancestor_v, mode, cs0, cs1, p1, p2, pid, proget, 2);
}

/* Build the ancestor likeliness for an inverse tau decay process with a
 * tau final state.
 */
static void ancestor_tau_inverse(struct ent_context * context,
    struct ent_state * daughter, double Z, int * np, int * proget_v, double * p,
    enum ent_pid * ancestor_v, int mode, double * cs0, double * cs1, double p1,
    double p2)
{
        const int pid[2] = { ENT_PID_NU_TAU, ENT_PID_NU_BAR_E };
        const int proget[2] = { PROGET_INVERSE_NU_TAU_TAU,
                PROGET_INVERSE_NU_BAR_E_TAU };

        ancestor_likeliness_fill(context, daughter, Z, np, proget_v, p,
            ancestor_v, mode, cs0, cs1, p1, p2, pid, proget, 2);
}

static void ancestor_decay(struct ent_context * context,
    struct ent_state * daughter, enum ent_pid mother, double A, double density,
    int * np, int * proget_v, double * p, enum ent_pid * ancestor_v)
{
        const double rho0 = context->ancestor(context, mother, daughter);
        if (rho0 > 0.) {
                double c;
                if (abs(mother) == ENT_PID_MUON) {
                        const double g = daughter->energy / ENT_MASS_MUON;
                        c = ENT_CTAU_MUON * sqrt(g * (2. + g));
                        proget_v[*np] = PROGET_BACKWARD_DECAY_MUON;
                } else {
                        const double g = daughter->energy / ENT_MASS_TAU;
                        c = ENT_CTAU_TAU * sqrt(g * (2. + g));
                        proget_v[*np] = PROGET_BACKWARD_DECAY_TAU;
                }
                const double p0 = (*np > 0) ? p[*np - 1] : 0.;
                p[*np] = p0 + rho0 * A * 1E-03 / (ENT_PHYS_NA * density * c);
                ancestor_v[(*np)++] = mother;
        }
}

/* Randomise the ancestor and the interaction process. */
static enum ent_return ancestor_draw(struct ent_context * context,
    struct ent_state * daughter, int np, int * proget_v, double * p,
    enum ent_pid * ancestor_v, enum ent_pid * ancestor, int * proget)
{
        if (np == 0) {
                /* No valid process. */
                *ancestor = ENT_PID_NONE;
                *proget = PROGET_N;
        } else if (np == 1) {
                *ancestor = ancestor_v[0];
                *proget = proget_v[0];
        } else {
                const double r = context->random(context) * p[np - 1];
                int i;
                for (i = 0; i < np; i++) {
                        if (r <= p[i]) break;
                }
                if (i > np - 1) i = np - 1;

                *ancestor = ancestor_v[i];
                *proget = proget_v[i];
                const double dp = (i == 0) ? p[0] : p[i] - p[i - 1];
                daughter->weight *= p[np - 1] / dp;
        }

        return ENT_RETURN_SUCCESS;
}

/* Randomise the ancestor at a backward vertex. */
static enum ent_return transport_ancestor_draw(struct ent_physics * physics,
    struct ent_context * context, struct ent_state * daughter,
    struct ent_medium * medium, double density, enum ent_process process,
    enum ent_pid * ancestor, int * proget, double * br_ptr)
{
        /* Build the interpolation or extrapolation factors for cross-sections.
         */
        double *cs0, *cs1;
        double p1, p2;
        const double mt = physics->mteff[4];
        const double Emin = 3 * mt * mt / (4 * ENT_MASS_NUCLEON);
        /* The top production cross-section is null below kinematic cut, in
         * which case the implemented BMC procedure would fail. Let us freeze
         * the final energy whenever this happens.
         */
        const double Ef =
            (daughter->energy > Emin) ? daughter->energy: Emin;
        const int mode = cross_section_prepare(
            physics->cs_t, Ef, &cs0, &cs1, &p1, &p2);

        /* Check the valid backward processes and compute their a priori
         * probabilities of occurence.
         */
        int proget_v[13] =
            { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
        int ancestor_v[13] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double p[13] =
            { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
        int np = 0;
        *br_ptr = 1.;

        const double Z = medium->Z;
        const double N = medium->A - Z;
        int apid = abs(daughter->pid);
        if (apid == ENT_PID_NU_E || apid == ENT_PID_NU_MU ||
            apid == ENT_PID_NU_TAU) {
                double rho0 = -1.;
                if ((process == ENT_PROCESS_NONE) ||
                    (process == ENT_PROCESS_DIS_NC)) {
                        /* Neutral current events. */
                        if (rho0 < 0.) {
                                rho0 = context->ancestor(
                                    context, daughter->pid, daughter);
                        }
                        if (rho0 > 0.) {
                                if (N > 0.) {
                                        /* Neutron target. */
                                        proget_v[np] = (daughter->pid > 0) ?
                                            PROGET_NC_NU_NEUTRON :
                                            PROGET_NC_NU_BAR_NEUTRON;
                                        ancestor_v[np] = daughter->pid;
                                        p[np] = (np ? p[np - 1] : 0.) +
                                            rho0 * N * cross_section_compute(
                                                mode, proget_v[np], cs0, cs1,
                                                p1, p2);
                                        if (p[np] > 0.) np++;
                                }

                                /* Proton target. */
                                if (Z > 0.) {
                                        proget_v[np] = (daughter->pid > 0) ?
                                            PROGET_NC_NU_PROTON :
                                            PROGET_NC_NU_BAR_PROTON;
                                        ancestor_v[np] = daughter->pid;
                                        p[np] = (np ? p[np - 1] : 0.) +
                                            rho0 * Z * cross_section_compute(
                                                mode, proget_v[np], cs0, cs1,
                                                p1, p2);
                                        if (p[np] > 0.) np++;
                                }
                        }
                }

                if (Z > 0.) {
                        /* Atomic electron(s) target. */
                        if ((process == ENT_PROCESS_NONE) ||
                            (process == ENT_PROCESS_ELASTIC)) {
                                /* Elastic event on an atomic electron. */
                                if (rho0 < 0.) {
                                        rho0 = context->ancestor(
                                            context, daughter->pid, daughter);
                                }
                                if (rho0 > 0.) {
                                        if (apid == ENT_PID_NU_E)
                                                proget_v[np] =
                                                    PROGET_ELASTIC_NU_E;
                                        else if (apid == ENT_PID_NU_MU)
                                                proget_v[np] =
                                                    PROGET_ELASTIC_NU_MU;
                                        else
                                                proget_v[np] =
                                                    PROGET_ELASTIC_NU_TAU;
                                        if (daughter->pid < 0) proget_v[np]++;
                                        ancestor_v[np] = daughter->pid;
                                        p[np] = (np ? p[np - 1] : 0.) +
                                            rho0 * Z * cross_section_compute(
                                                mode, proget_v[np], cs0, cs1,
                                                p1, p2);
                                        if (p[np] > 0.) np++;
                                }
                        }

                        /* Inverse decay processes. */
                        if (daughter->pid == ENT_PID_NU_E) {
                                double rho1 = 0.;
                                if ((process == ENT_PROCESS_NONE) ||
                                    (process == ENT_PROCESS_INVERSE_MUON)) {
                                        rho1 = context->ancestor(
                                            context, ENT_PID_NU_MU, daughter);
                                }
                                if (rho1 > 0.) {
                                        proget_v[np] = PROGET_INVERSE_NU_MU_MU;
                                        ancestor_v[np] = ENT_PID_NU_MU;
                                        p[np] = (np ? p[np - 1] : 0.) +
                                            rho1 * Z * cross_section_compute(
                                                mode, proget_v[np], cs0, cs1,
                                                p1, p2);
                                        if (p[np] > 0.) np++;
                                }
                                double rho2 = 0.;
                                if ((process == ENT_PROCESS_NONE) ||
                                    (process == ENT_PROCESS_INVERSE_TAU)) {
                                        rho2 = context->ancestor(
                                            context, ENT_PID_NU_TAU, daughter);
                                }
                                if (rho2 > 0.) {
                                        proget_v[np] =
                                            PROGET_INVERSE_NU_TAU_TAU;
                                        ancestor_v[np] = ENT_PID_NU_TAU;
                                        p[np] = (np ? p[np - 1] : 0.) +
                                            rho2 * Z * cross_section_compute(
                                                mode, proget_v[np], cs0, cs1,
                                                p1, p2);
                                        if (p[np] > 0.) np++;
                                }
                        } else if (daughter->pid == ENT_PID_NU_BAR_MU) {
                                double rho1 = 0.;
                                if ((process == ENT_PROCESS_NONE) ||
                                    (process == ENT_PROCESS_INVERSE_MUON)) {
                                        rho1 = context->ancestor(context,
                                            ENT_PID_NU_BAR_E, daughter);
                                }
                                if (rho1 > 0.) {
                                        proget_v[np] =
                                            PROGET_INVERSE_NU_BAR_E_MU;
                                        ancestor_v[np] = ENT_PID_NU_BAR_E;
                                        p[np] = (np ? p[np - 1] : 0.) +
                                            rho0 * Z * cross_section_compute(
                                                mode, proget_v[np], cs0, cs1,
                                                p1, p2);
                                        if (p[np] > 0.) np++;
                                }
                        } else if (daughter->pid == ENT_PID_NU_BAR_TAU) {
                                double rho1 = 0.;
                                if ((process == ENT_PROCESS_NONE) ||
                                    (process == ENT_PROCESS_INVERSE_TAU)) {
                                        rho1 = context->ancestor(context,
                                            ENT_PID_NU_BAR_E, daughter);
                                }
                                if (rho1 > 0.) {
                                        proget_v[np] =
                                            PROGET_INVERSE_NU_BAR_E_TAU;
                                        ancestor_v[np] = ENT_PID_NU_BAR_E;
                                        p[np] = (np ? p[np - 1] : 0.) +
                                            rho0 * Z * cross_section_compute(
                                                mode, proget_v[np], cs0, cs1,
                                                p1, p2);
                                        if (p[np] > 0.) np++;
                                }
                        }
                }

                /* Neutrino from top decay in a CC DIS. */
                if ((process == ENT_PROCESS_NONE) ||
                    (process == ENT_PROCESS_DIS_CC) ||
                    (process == ENT_PROCESS_DIS_CC_TOP)) {
                        double rho1 = 0.;
                        int j;
                        for (j = 0; j < 3; j++) {
                                const int pid =
                                    (daughter->pid > 0) ? 12 + 2 * j :
                                        -12 - 2 * j;
                                const double d = context->ancestor(
                                    context, pid, daughter);
                                if (d > 0.) rho1 += d;
                        }
                        if (rho1 > 0.) {
                                /* CC processes with top production and
                                 * leptonic W decay.
                                 */
                                double br;
                                if (apid == ENT_PID_NU_E) {
                                        br = ENT_BR_W_TO_ELECTRON;
                                } else if (apid == ENT_PID_NU_MU) {
                                        br = ENT_BR_W_TO_MUON;
                                } else {
                                        br = ENT_BR_W_TO_TAU;
                                }
                                *br_ptr = br;

                                /* Neutron target. */
                                if (N > 0.) {
                                        if (daughter->pid > 0) {
                                                proget_v[np] =
                                                    PROGET_CC_TOP_NU_NEUTRON;
                                                ancestor_v[np] =
                                                    ENT_PID_TOP;
                                        } else {
                                                proget_v[np] =
                                                    PROGET_CC_TOP_NU_BAR_NEUTRON;
                                                ancestor_v[np] =
                                                    ENT_PID_TOP_BAR;
                                        }
                                        p[np] = (np > 0 ? p[np - 1] : 0.) +
                                            rho1 * N * br *
                                            cross_section_compute(
                                                mode, proget_v[np], cs0, cs1,
                                                p1, p2);
                                        if (p[np] > 0.) np++;
                                }

                                /* Proton target. */
                                if (Z > 0.) {
                                        if (daughter->pid > 0) {
                                                proget_v[np] =
                                                    PROGET_CC_TOP_NU_PROTON;
                                                ancestor_v[np] =
                                                    ENT_PID_TOP;
                                        } else {
                                                proget_v[np] =
                                                    PROGET_CC_TOP_NU_BAR_PROTON;
                                                ancestor_v[np] =
                                                    ENT_PID_TOP_BAR;
                                        }
                                        p[np] = (np > 0 ? p[np - 1] : 0.) +
                                            rho1 * Z * br *
                                            cross_section_compute(
                                                mode, proget_v[np], cs0, cs1,
                                                p1, p2);
                                        if (p[np] > 0.) np++;
                                }
                        }
                }

                /* True decay processes. */
                if (medium->density != NULL) {
                        /* XXX Check ancestor density? */
                        if (density <= 0.) {
                                medium->density(medium, daughter, &density);
                                if (density <= 0.)
                                        return ENT_RETURN_DOMAIN_ERROR;
                        }
                        if (apid != ENT_PID_NU_TAU) {
                                int mother;
                                if (apid == ENT_PID_MUON) {
                                        mother = (daughter->pid > 0) ?
                                            ENT_PID_MUON :
                                            ENT_PID_MUON_BAR;
                                } else {
                                        mother = (daughter->pid > 0) ?
                                            ENT_PID_MUON_BAR :
                                            ENT_PID_MUON;
                                }
                                ancestor_decay(context, daughter, mother,
                                    medium->A, density, &np, proget_v, p,
                                    ancestor_v);
                        }

                        /* True decay process from a tau. */
                        int mother;
                        if (apid == ENT_PID_NU_TAU) {
                                mother = (daughter->pid > 0) ? ENT_PID_TAU :
                                                               ENT_PID_TAU_BAR;
                        } else {
                                mother = (daughter->pid > 0) ? ENT_PID_TAU_BAR :
                                                               ENT_PID_TAU;
                        }
                        ancestor_decay(context, daughter, mother, medium->A,
                            density, &np, proget_v, p, ancestor_v);
                }
        } else if ((apid == ENT_PID_ELECTRON) || (apid == ENT_PID_MUON) ||
            (apid == ENT_PID_TAU)) {
                const int npid = (daughter->pid > 0) ? apid + 1 : -apid - 1;
                double rho0 = -1.;
                if ((process == ENT_PROCESS_NONE) ||
                    (process == ENT_PROCESS_DIS_CC) ||
                    (process == ENT_PROCESS_DIS_CC_OTHER)) {
                        /* Charged current processes (no top). */
                        if (rho0 < 0) {
                                rho0 = context->ancestor(
                                    context, npid, daughter);
                        }
                        if (rho0 > 0.) {
                                /* Neutron target. */
                                if (N > 0.) {
                                        proget_v[np] = (daughter->pid > 0) ?
                                            PROGET_CC_OTHER_NU_NEUTRON :
                                            PROGET_CC_OTHER_NU_BAR_NEUTRON;
                                        ancestor_v[np] = npid;
                                        p[np] = (np ? p[np - 1] : 0.) +
                                            rho0 * N * cross_section_compute(
                                                mode, proget_v[np], cs0, cs1,
                                                p1, p2);
                                        if (p[np] > 0.) np++;
                                }

                                /* Proton target. */
                                if (Z > 0.) {
                                        proget_v[np] = (daughter->pid > 0) ?
                                            PROGET_CC_OTHER_NU_PROTON :
                                            PROGET_CC_OTHER_NU_BAR_PROTON;
                                        ancestor_v[np] = npid;
                                        p[np] = (np ? p[np - 1] : 0.) +
                                            rho0 * Z * cross_section_compute(
                                                mode, proget_v[np], cs0, cs1,
                                                p1, p2);
                                        if (p[np] > 0.) np++;
                                }
                        }
                }

                if ((process == ENT_PROCESS_NONE) ||
                    (process == ENT_PROCESS_DIS_CC) ||
                    (process == ENT_PROCESS_DIS_CC_TOP)) {
                        /* Charged current processes (w/ top). */
                        if (rho0 < 0) {
                                rho0 = context->ancestor(
                                    context, npid, daughter);
                        }
                        if (rho0 > 0.) {
                                /* Neutron target. */
                                if (N > 0.) {
                                        proget_v[np] = (daughter->pid > 0) ?
                                            PROGET_CC_TOP_NU_NEUTRON :
                                            PROGET_CC_TOP_NU_BAR_NEUTRON;
                                        ancestor_v[np] = npid;
                                        p[np] = (np ? p[np - 1] : 0.) +
                                            rho0 * N * cross_section_compute(
                                                mode, proget_v[np], cs0, cs1,
                                                p1, p2);
                                        if (p[np] > 0.) np++;
                                }

                                /* Proton target. */
                                if (Z > 0.) {
                                        proget_v[np] = (daughter->pid > 0) ?
                                            PROGET_CC_TOP_NU_PROTON :
                                            PROGET_CC_TOP_NU_BAR_PROTON;
                                        ancestor_v[np] = npid;
                                        p[np] = (np ? p[np - 1] : 0.) +
                                            rho0 * Z * cross_section_compute(
                                                mode, proget_v[np], cs0, cs1,
                                                p1, p2);
                                        if (p[np] > 0.) np++;
                                }
                        }

                        /* Lepton from top decay in a CC DIS event. */
                        double rho1 = 0.;
                        int j;
                        for (j = 0; j < 3; j++) {
                                const int pid =
                                    (daughter->pid > 0) ? -12 - 2 * j :
                                        12 + 2 * j;
                                const double d = context->ancestor(
                                    context, pid, daughter);
                                if (d > 0.) rho1 += d;
                        }
                        if (rho1 > 0.) {
                                /* CC processes with top production and
                                 * leptonic W decay.
                                 */
                                double br;
                                if (apid == ENT_PID_ELECTRON) {
                                        br = ENT_BR_W_TO_ELECTRON;
                                } else if (apid == ENT_PID_MUON) {
                                        br = ENT_BR_W_TO_MUON;
                                } else {
                                        br = ENT_BR_W_TO_TAU;
                                }
                                *br_ptr = br;

                                /* Neutron target. */
                                if (N > 0.) {
                                        if (daughter->pid > 0) {
                                                proget_v[np] =
                                                    PROGET_CC_TOP_NU_BAR_NEUTRON;
                                                ancestor_v[np] =
                                                    ENT_PID_TOP_BAR;
                                        } else {
                                                proget_v[np] =
                                                    PROGET_CC_TOP_NU_NEUTRON;
                                                ancestor_v[np] =
                                                    ENT_PID_TOP;
                                        }
                                        p[np] = (np > 0 ? p[np - 1] : 0.) +
                                            rho1 * N * br *
                                            cross_section_compute(
                                                mode, proget_v[np], cs0, cs1,
                                                p1, p2);
                                        if (p[np] > 0.) np++;
                                }

                                /* Proton target. */
                                if (Z > 0.) {
                                        if (daughter->pid > 0) {
                                                proget_v[np] =
                                                    PROGET_CC_TOP_NU_BAR_PROTON;
                                                ancestor_v[np] =
                                                    ENT_PID_TOP_BAR;
                                        } else {
                                                proget_v[np] =
                                                    PROGET_CC_TOP_NU_PROTON;
                                                ancestor_v[np] =
                                                    ENT_PID_TOP;
                                        }
                                        p[np] = (np > 0 ? p[np - 1] : 0.) +
                                            rho1 * Z * br *
                                            cross_section_compute(
                                                mode, proget_v[np], cs0, cs1,
                                                p1, p2);
                                        if (p[np] > 0.) np++;
                                }
                        }
                }

                if (Z > 0.) {
                        /* Collision with atomic electron(s). */
                        if (daughter->pid == ENT_PID_ELECTRON) {
                                /* Elastic processes on an atomic electron. */
                                if ((process == ENT_PROCESS_NONE) ||
                                    (process == ENT_PROCESS_ELASTIC)) {
                                        ancestor_electron_elastic(context,
                                            daughter, Z, &np, proget_v, p,
                                            ancestor_v, mode, cs0, cs1, p1, p2);
                                }
                        } else if (daughter->pid == ENT_PID_MUON) {
                                /* Inverse muon decay process. */
                                if ((process == ENT_PROCESS_NONE) ||
                                    (process == ENT_PROCESS_INVERSE_MUON)) {
                                        ancestor_muon_inverse(context, daughter,
                                            Z, &np, proget_v, p, ancestor_v,
                                            mode, cs0, cs1, p1, p2);
                                }
                        } else if (daughter->pid == ENT_PID_TAU) {
                                /* Inverse tau decay process. */
                                if ((process == ENT_PROCESS_NONE) ||
                                    (process == ENT_PROCESS_INVERSE_TAU)) {
                                        ancestor_tau_inverse(context, daughter,
                                            Z, &np, proget_v, p, ancestor_v,
                                            mode, cs0, cs1, p1, p2);
                                }
                        }
                }
        } else
                return ENT_RETURN_DOMAIN_ERROR;

        /* Randomise the ancestor and the interaction process. */
        return ancestor_draw(
            context, daughter, np, proget_v, p, ancestor_v, ancestor, proget);
}

/* Compute the p and n cross-sections for a DIS standalone vertex. */
static enum ent_return vertex_dis_compute(struct ent_physics * physics,
    enum ent_pid pid, double energy, struct ent_medium * medium,
    enum ent_process process, int * proget_p_, int * proget_n_, double * cs_p,
    double * cs_n)
{
        enum ent_return rc;
        int proget_p, proget_n;
        if ((rc = proget_compute(pid, ENT_PID_PROTON, process, &proget_p)) !=
            ENT_RETURN_SUCCESS)
                return rc;
        *proget_p_ = proget_p;
        if ((rc = proget_compute(pid, ENT_PID_NEUTRON, process, &proget_n)) !=
            ENT_RETURN_SUCCESS)
                return rc;
        *proget_n_ = proget_n;

        /* Let us build the interpolation or extrapolation
         * factors for the cross-sections.
         */
        double *cs0, *cs1;
        double p1, p2;
        int mode = cross_section_prepare(
            physics->cs_t, energy, &cs0, &cs1, &p1, &p2);

        /* Compute the relevant cross-sections and randomise the
         * target accordingly.
         */
        *cs_p = (medium->Z > 0.) ? medium->Z *
                cross_section_compute(mode, proget_p, cs0, cs1, p1, p2) :
                                   0.;
        const double N = medium->A - medium->Z;
        *cs_n = (N > 0.) ?
            N * cross_section_compute(mode, proget_n, cs0, cs1, p1, p2) :
            0.;

        if (*cs_p + *cs_n <= 0.) { /* Fallback in case of null cross-section. */
                *cs_p = medium-> Z / medium->A;
                *cs_n = 1. - *cs_p;
        }

        return ENT_RETURN_SUCCESS;
}

/* Randomise the proget index for a DIS standalone vertex.  */
static enum ent_return vertex_dis_randomise(struct ent_physics * physics,
    struct ent_context * context, struct ent_state * state,
    struct ent_medium * medium, enum ent_process process, int * proget)
{
        /* Get the mother's pid. */
        enum ent_pid pid;
        if ((process == ENT_PROCESS_DIS_NC) || (context->ancestor == NULL))
                pid = state->pid;
        else {
                pid = abs(state->pid) + 1;
                if (state->pid < 0) pid = -pid;
        }

        /* Compute the p and n cross-sections. */
        enum ent_return rc;
        int proget_p, proget_n;
        double cs_p, cs_n;
        if ((rc = vertex_dis_compute(physics, pid, state->energy, medium,
                 process, &proget_p, &proget_n, &cs_p, &cs_n)) !=
            ENT_RETURN_SUCCESS)
                return rc;

        /* Randomise the target accordingly. */
        *proget = (context->random(context) < cs_p / (cs_p + cs_n)) ? proget_p :
                                                                      proget_n;

        if (context->ancestor != NULL) {
                /* Update the BMC weight. */
                const double cs = (*proget == proget_p) ? cs_p : cs_n;
                state->weight *= (cs_p + cs_n) / cs;
        }

        return ENT_RETURN_SUCCESS;
}

/* Vertex randomisation in forward Monte-Carlo. */
static enum ent_return vertex_forward(struct ent_physics * physics,
    struct ent_context * context, struct ent_state * state,
    struct ent_medium * medium, enum ent_process process,
    struct ent_state * products)
{
        /* Get the process and target index. */
        int proget;
        if ((process == ENT_PROCESS_NONE) ||
            (process == ENT_PROCESS_DIS_CC)) {
                /* Randomise the process and the target if not specified.
                 * First let's compute the relevant cross-sections.
                 */
                enum ent_return rc;
                double cs[PROGET_N];
                if ((rc = transport_cross_section(physics, state->pid,
                         state->energy, medium->Z, medium->A, process, cs)) !=
                    ENT_RETURN_SUCCESS)
                        return rc;

                /* Then, let us randomise the interaction process and its
                 * corresponding target.
                 */
                const double r = cs[PROGET_N - 1] * context->random(context);
                if (r < 0.) return ENT_RETURN_DOMAIN_ERROR;
                for (proget = 0; proget < PROGET_N; proget++)
                        if (r <= cs[proget]) break;
        } else if ((process == ENT_PROCESS_DIS_CC_OTHER) ||
            (process == ENT_PROCESS_DIS_CC_TOP) ||
            (process == ENT_PROCESS_DIS_NC)) {
                /* This is a DIS event on a proton or neutron. Let's get the
                 * corresponding indices.
                 */
                enum ent_return rc;
                if ((rc = vertex_dis_randomise(physics, context, state, medium,
                         process, &proget)) != ENT_RETURN_SUCCESS)
                        return rc;
        } else {
                /* This is an interaction with an atomic electron. Let's get
                 * the corresponding index.
                 */
                enum ent_return rc;
                if ((rc = proget_compute(state->pid, ENT_PID_ELECTRON, process,
                         &proget)) != ENT_RETURN_SUCCESS)
                        return rc;
        }

        /* Process the vertex. */
        return transport_vertex(physics, context, proget, state, products);
}

/* Vertex randomisation in backward Monte-Carlo. */
static enum ent_return vertex_backward(struct ent_physics * physics,
    struct ent_context * context, struct ent_state * state,
    struct ent_medium * medium, double density, enum ent_process process,
    struct ent_state * products, enum proget_index * proget_)
{
        /* Get the ancestor, the process and the target. */
        enum ent_pid ancestor;
        int proget;
        double br;
        /* Randomise the ancestor, the process and the target. */
        enum ent_return rc;
        if ((rc = transport_ancestor_draw(physics, context, state,
                 medium, density, process, &ancestor, &proget, &br)) !=
            ENT_RETURN_SUCCESS)
                return rc;

        if (proget == PROGET_N) {
                /* No valid process. */
                return ENT_RETURN_SUCCESS;
        }

        /* Process the vertex. */
        if ((proget >= 8) && (proget < PROGET_N_DIS) &&
            abs(ancestor) == ENT_PID_TOP) {
                /* Backward decay from a top, from DIS CC. */
                rc = transport_vertex_backward_top_decay(
                    physics, context, proget, state, products);
        } else {
                rc = transport_vertex(physics, context, proget, state,
                    products);
        }
        if (rc != ENT_RETURN_SUCCESS) {
                return rc;
        }

        /* Apply any biasing weight for the ancestor and for the process
         * randomisation.
         *
         * N.B.: if the ancestor is a muon or a tau, then the only possible
         * source process is a decay. Thus p_true = 1. and there is no need to
         * further correct the BMC weight.
         */
        if ((proget >= 0) && (state->weight > 0.)) {
                enum ent_return rc;
                double cs[PROGET_N];
                if ((rc = transport_cross_section(physics, state->pid,
                         state->energy, medium->Z, medium->A, process, cs)) !=
                    ENT_RETURN_SUCCESS) {
                        return rc;
                }

                const double d =
                    (proget == 0) ? cs[proget] : cs[proget] - cs[proget - 1];
                state->weight *= br * d / cs[PROGET_N - 1];
        }

        if (products != NULL) {
                int i;
                for (i = 0; i < ENT_PRODUCTS_SIZE; i++) {
                        products[i].weight = (products->pid != ENT_PID_NONE) ?
                            state->weight : 0.;
                }
        }

        if (proget_ != NULL) *proget_ = proget;

        return ENT_RETURN_SUCCESS;
}

/* API interface for a Monte-Carlo collision. */
enum ent_return ent_collide(struct ent_physics * physics,
    struct ent_context * context, struct ent_state * state,
    struct ent_medium * medium, enum ent_process process,
    struct ent_state * products)
{
        ENT_ACKNOWLEDGE(ent_collide);

        if (products != NULL) {
                memset(products, 0x0, ENT_PRODUCTS_SIZE * sizeof(*products));
        }

        /* Check and format the inputs. */
        if ((physics == NULL) || (context == NULL) ||
            (context->random == NULL) || (state == NULL) || (medium == NULL))
                ENT_RETURN(ENT_RETURN_BAD_ADDRESS);
        if (medium->A < medium->Z) ENT_RETURN(ENT_RETURN_DOMAIN_ERROR);

        /* Process the vertex. */
        if (context->ancestor == NULL)
                ENT_RETURN(vertex_forward(
                    physics, context, state, medium, process, products));
        else
                ENT_RETURN(vertex_backward(physics, context, state, medium, -1.,
                    process, products, NULL));

        return ENT_RETURN_SUCCESS;
}

/* API interface for a Monte-Carlo tranport. */
enum ent_return ent_transport(struct ent_physics * physics,
    struct ent_context * context, struct ent_state * state,
    struct ent_state * products, enum ent_event * event)
{
        ENT_ACKNOWLEDGE(ent_transport);

        if (products != NULL) {
                memset(products, 0x0, ENT_PRODUCTS_SIZE * sizeof(*products));
        }
        *event = ENT_EVENT_NONE;

        /* Check and format the inputs. */
        if ((context == NULL) || (context->medium == NULL) ||
            (context->random == NULL) || (state == NULL))
                ENT_RETURN(ENT_RETURN_BAD_ADDRESS);

        /* Check for an initial limit violation. */
        enum ent_return rc = ENT_RETURN_SUCCESS;
        enum ent_event event_ = ENT_EVENT_NONE;
        if ((context->distance_max > 0.) &&
            (state->distance >= context->distance_max)) {
                event_ = ENT_EVENT_LIMIT_DISTANCE;
                goto exit;
        }
        if ((context->grammage_max > 0.) &&
            (state->grammage >= context->grammage_max)) {
                event_ = ENT_EVENT_LIMIT_GRAMMAGE;
                goto exit;
        }

        /* Initialise the transport. */
        double step = 0., density;
        struct ent_medium * medium = NULL;
        if ((rc = transport_step(context, state, &medium, &step, &density,
                 context->grammage_max, &event_)) != ENT_RETURN_SUCCESS)
                goto exit;
        if (event_ != ENT_EVENT_NONE) goto exit;

        /* Set any grammage limit. */
        enum ent_event foreseen = ENT_EVENT_NONE;
        double cs[PROGET_N], Xlim = 0., Xint = 0.;
        if (physics != NULL) {
                /* Compute the cross-sections. */
                if ((rc = transport_cross_section(physics, state->pid,
                         state->energy, medium->Z, medium->A, ENT_PROCESS_NONE,
                         cs)) !=
                    ENT_RETURN_SUCCESS)
                        goto exit;

                /* Randomise the depth of the next interaction. */
                Xint = medium->A * 1E-03 / (cs[PROGET_N - 1] * ENT_PHYS_NA);
                Xlim = state->grammage - Xint * log(context->random(context));

                if ((context->grammage_max > 0.) &&
                    (context->grammage_max < Xlim)) {
                        Xlim = context->grammage_max;
                        foreseen = ENT_EVENT_LIMIT_GRAMMAGE;
                }
        } else if (context->grammage_max > 0.) {
                Xlim = context->grammage_max;
                foreseen = ENT_EVENT_LIMIT_GRAMMAGE;
        }

        /* Call any custom stepping action. */
        if (context->stepping_action != NULL) {
                if ((rc = context->stepping_action(context, medium, state)) !=
                    ENT_RETURN_SUCCESS)
                        goto exit;
        }

        /* Do the stepping. */
        if (step)
                for (;;) {
                        /* Step until an event occurs. */
                        struct ent_medium * m = medium;
                        if ((rc = transport_step(context, state, &medium, &step,
                                 &density, Xlim, &event_)) !=
                            ENT_RETURN_SUCCESS)
                                goto exit;
                        if (event_ != ENT_EVENT_NONE) break;
                        if ((medium != NULL) && (step == 0.)) {
                                rc = ENT_RETURN_DOMAIN_ERROR;
                                goto exit;
                        }

                        /* Call any custom stepping action if a medium
                         * change occured.
                         */
                        if ((context->stepping_action != NULL) &&
                            (medium != m)) {
                                if ((rc = context->stepping_action(context,
                                         medium, state)) != ENT_RETURN_SUCCESS)
                                        goto exit;
                        }
                }
        else {
                /* This is a uniform medium of infinite extension.
                 * Let's do a single straight step.
                 */
                event_ = transport_straight(context, state, density, Xlim);
        }

        /* Call any custom stepping action. */
        if (context->stepping_action != NULL) {
                if ((rc = context->stepping_action(context, medium, state)) !=
                    ENT_RETURN_SUCCESS)
                        goto exit;
        }

        /* Process any interaction. */
        if ((event_ == ENT_EVENT_LIMIT_GRAMMAGE) &&
            (foreseen == ENT_EVENT_NONE)) {
                event_ = ENT_EVENT_INTERACTION;
                if (context->ancestor == NULL) {
                        /* This is a forward transport. Let us randomise the
                         * interaction process and the target.
                         */
                        const double r =
                            cs[PROGET_N - 1] * context->random(context);
                        if (r < 0.)
                                rc = ENT_RETURN_DOMAIN_ERROR;
                        else {
                                int proget;
                                for (proget = 0; proget < PROGET_N; proget++)
                                        if (r <= cs[proget]) break;
                                rc = transport_vertex(
                                    physics, context, proget, state, products);
                        }
                } else {
                        /* This is a backward transport. First let us randomise
                         * the ancestor and the source process.
                         */
                        enum proget_index proget;
                        rc = vertex_backward(physics, context, state, medium,
                            density, ENT_PROCESS_NONE, products, &proget);

                        /* Then let us apply the effective weight for the
                         * transport.
                         */
                        if (state->weight > 0.) {
                                double X0;
                                if (proget >= 0) {
                                        /* The ancestor is a neutrino. Let us
                                         * apply a flux like boundary condition
                                         * at the vertex.
                                         */
                                        if ((rc = transport_cross_section(
                                            physics, state->pid, state->energy,
                                            medium->Z, medium->A,
                                            ENT_PROCESS_NONE, cs)) !=
                                            ENT_RETURN_SUCCESS)
                                                goto exit;
                                        X0 = medium->A * 1E-03 /
                                            (cs[PROGET_N - 1] * ENT_PHYS_NA);
                                } else {
                                        /* The neutrino originates from a muon
                                         * or tau decay. A vertex boundary
                                         * condition is used since the initial
                                         * state is not known at that point.
                                         */
                                        event_ = (proget ==
                                            PROGET_BACKWARD_DECAY_MUON) ?
                                            ENT_EVENT_DECAY_MUON :
                                            ENT_EVENT_DECAY_TAU;
                                        X0 = density;
                                }
                                state->weight *= Xint / X0;
                        }
                }
        }

exit:
        if (event != NULL) *event = event_;
        ENT_RETURN(rc);
}

#if (_USE_GDB)
/*  For debugging with gdb, on linux. */
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <fenv.h>

/* Flag for floating point exceptions. */
static int fe_status;

/* Library initialisation. */
void __attribute__((constructor)) __init(void)
{
        /* Save the floating points exceptions status and enable them.
         */
        fe_status = fegetexcept();
        feclearexcept(FE_ALL_EXCEPT);
        feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
}

/* Library finalisation. */
void __attribute__((destructor)) __fini(void)
{
        /* Restore the floating points exceptions status. */
        feclearexcept(FE_ALL_EXCEPT);
        feenableexcept(fe_status);
}
#endif
