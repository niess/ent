/*
 *  an Engine for Neutrinos Transport (ENT)
 *  Copyright (C) 2016  Valentin Niess
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRENTY; without even the implied warranty of
 *  MERCHENTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* Standard library includes. */
#include <ctype.h>
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
/* The W boson partial width to muon + nu_mu, in GeV/c^2. */
#define ENT_WIDTH_W_TO_MUON 0.22164
/* The Fermi coupling constant GF/(hbar*c)^3, in GeV^-2. */
#define ENT_PHYS_GF 1.1663787E-05
/* The Planck constant as hbar*c, in GeV * m. */
#define ENT_PHYS_HBC 1.97326978E-16
/* The Weinberg angle at MZ, as sin(theta_W)^2. */
#define ENT_PHYS_SIN_THETA_W_2 0.231295

/* Kinetic range for tabulations. */
#define KINETIC_MIN 1E+03
#define KINETIC_MAX 1E+12
#define KINETIC_N 181

#ifndef M_PI
/* Define pi, if unknown. */
#define M_PI 3.14159265358979323846
#endif

/* Opaque Physics object. */
struct ent_physics {
        /* Index to the PDF data. */
        struct lha_pdf * pdf;
        /* Entry point for the tabulations. */
        double * table;
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
static void * physics_create(int extra_size)
{
        void * v = malloc(sizeof(struct ent_physics) + extra_size +
            35 * KINETIC_N * sizeof(double));
        if (v == NULL) return NULL;
        struct ent_physics * p = v;
        p->table = (double *)(p->data + extra_size);
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
        REGISTER_FUNCTION(ent_physics_dcs);
        REGISTER_FUNCTION(ent_physics_pdf);

        /* Other API functions. */
        REGISTER_FUNCTION(ent_physics_destroy);
        REGISTER_FUNCTION(ent_error_string);
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
                        void * tmp =
                            realloc(*buffer, sizeof(**buffer) + (*buffer)->size
                                + FILE_BLOCK_SIZE);
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

/* Container for LHAPDFs data. */
struct lha_pdf {
        /* The table size. */
        int nx, nQ2, nf;
        /* Links to the data tables. */
        float * x, * Q2, * lambda, * xfx;
        /* Placeholder for variable size data. */
        float data[];
};

/* Load PDFs from a .dat file in lhagrid1 format. */
static enum ent_return lha_load(FILE * stream, struct ent_physics ** physics)
{
#define LHAPDF_NF_MAX 13
        *physics = NULL;
        enum ent_return rc;
        struct file_buffer * buffer = NULL;

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

        /* Locate the data segment. */
        for (;;) {
                file_get_line(&buffer, 0);
                if (strlen(buffer->cursor) < 3) continue;
                if ((buffer->cursor[0] == '-') && (buffer->cursor[1] == '-') &&
                    (buffer->cursor[2] == '-')) break;
        }
        long int pos = ftell(stream);

        /* Parse the table format. */
        file_get_line(&buffer, 0);
        const int nx = file_count_words(buffer);
        file_get_line(&buffer, 0);
        const int nQ = file_count_words(buffer);
        file_get_line(&buffer, 0);
        const int nf = file_count_words(buffer);
        if ((nf != LHAPDF_NF_MAX-2) && (nf != LHAPDF_NF_MAX)) {
                rc = ENT_RETURN_FORMAT_ERROR;
                goto exit;
        }

        /* Allocate and map the memory for the tables. */
        const unsigned int size_x =
            memory_padded_size(nx * sizeof(float), sizeof(double));
        const unsigned int size_Q =
            memory_padded_size(nQ * sizeof(float), sizeof(double));
        const unsigned int size_lambda =
                memory_padded_size(nf * nQ * sizeof(float), sizeof(double));
        const unsigned int size_xfx =
            memory_padded_size(nf * nx * nQ * sizeof(float), sizeof(double));
        unsigned int size = sizeof(struct lha_pdf) + size_x + size_Q +
            size_lambda + size_xfx;
        if ((*physics = physics_create(size)) == NULL) {
                rc = ENT_RETURN_MEMORY_ERROR;
                goto exit;
        }
        (*physics)->pdf = (struct lha_pdf *)(*physics)->data;
        struct lha_pdf * const pdf = (*physics)->pdf;
        pdf->nx = nx;
        pdf->nQ2 = nQ;
        pdf->nf = (nf - 1) / 2;
        pdf->x = pdf->data;
        pdf->Q2 = (float *)((char *)pdf->data + size_x);
        pdf->lambda = (float *)((char *)pdf->Q2 + size_Q);
        pdf->xfx = (float *)((char *)pdf->lambda + size_lambda);

        /* Roll back in the file and copy the data to memory. */
        if (fseek(stream, pos, SEEK_SET) != 0) {
                rc = ENT_RETURN_IO_ERROR;
                goto exit;
        }
        file_get_line(&buffer, 0);

        float row[LHAPDF_NF_MAX];
        int index[LHAPDF_NF_MAX];
        file_get_table(&buffer, nx, pdf->x);
        file_get_table(&buffer, nQ, pdf->Q2);
        file_get_table(&buffer, nf, row);

        int i;
        for (i = 0; i < nQ; i++) pdf->Q2[i] *= pdf->Q2[i];
        for (i = 0; i < nf; i++) {
                if (row[i] == 21.) index[i] = pdf->nf;
                else index[i] = (int)row[i] + pdf->nf;
                if ((index[i] < 0) || (index[i] >= nf)) {
                        rc = ENT_RETURN_FORMAT_ERROR;
                        goto exit;
                }
        }

        float * table;
        for (i = 0, table = pdf->xfx; i < nx * nQ; i++, table += nf) {
                int j;
                file_get_table(&buffer, nf, row);
                for (j = 0; j < nf; j++) table[index[j]] = row[j];
        }

        /* Tabulate the lambda exponents for the small x extrapolation. */
        memset(pdf->lambda, 0x0, size_lambda);
        if (pdf->x[0] <= 0.) goto exit;
        float lxi = log(pdf->x[1] / pdf->x[0]);
        if (lxi <= 0.) {
                rc = ENT_RETURN_FORMAT_ERROR;
                goto exit;
        }
        lxi = 1. / lxi;
        for (i = 0, table = pdf->lambda; i < nQ; i++, table += nf) {
                const float * const xfx0 = pdf->xfx + i * nf;
                const float * const xfx1 = pdf->xfx + (nQ + i) * nf;
                int j;
                for (j = 0; j < nf; j++) {
                        const float f0 = xfx0[j];
                        const float f1 = xfx1[j];
                        if ((f0 > 0.) && (f1 > 0.))
                                table[j] = log(f0 / f1) * lxi;
                }
        }

exit:
        free(buffer);
        if (rc != ENT_RETURN_SUCCESS) {
                free(*physics);
                *physics = NULL;
        }
        return rc;
}

/* Recursive bracketing of a table value, using a dichotomy. */
static void table_bracket(
    const float * table, float value, int * p1, int * p2)
{
        int i3 = (*p1 + *p2) / 2;
        if (value >= table[i3]) *p1 = i3;
        else *p2 = i3;
        if (*p2 - *p1 >= 2) table_bracket(table, value, p1, p2);
}

static void lha_pdf_compute(
    const struct lha_pdf * pdf, float x, float Q2, float * xfx)
{
        memset(xfx, 0x0, LHAPDF_NF_MAX * sizeof(*xfx));

        /* Check the bounds. */
        if ((Q2 < pdf->Q2[0]) || (Q2 >= pdf->Q2[pdf->nQ2-1]) ||
            (x >= pdf->x[pdf->nx-1])) return;

        if (x < pdf->x[0]) {
                /* Extrapolate with a power law for small x values, i.e.
                 * x * f(x) ~ 1 / x**lambda(Q2).
                 */
                int iQ0, iQ1;
                float hQ;
                if (Q2 >= pdf->Q2[pdf->nQ2-1]) {
                        iQ0 = iQ1 = pdf->nQ2 -1;
                        hQ = 1.;
                }
                else {
                        iQ0 = 0;
                        iQ1 = pdf->nQ2 -1;
                        table_bracket(pdf->Q2, Q2, &iQ0, &iQ1);
                        hQ = (Q2 - pdf->Q2[iQ0]) / (pdf->Q2[iQ1] - pdf->Q2[iQ0]);
                }

                int i, nf = 2 * pdf->nf + 1;
                const float * const lambda0 = pdf->lambda + iQ0 * nf;
                const float * const lambda1 = pdf->lambda + iQ1 * nf;
                const float * const xfx0 = pdf->xfx + iQ0 * nf;
                const float * const xfx1 = pdf->xfx + iQ1 * nf;
                for (i = 0; i < nf; i++) {
                        const float y0 = lambda0[i];
                        const float y1 = lambda1[i];
                        if ((y0 <= 0.) || (y1 <= 0.)) xfx[i] = 0.;
                        else {
                                const float lambda = y0 * (1. - hQ) + y1 * hQ;
                                xfx[i] = (xfx0[i] * (1. - hQ) + xfx1[i] * hQ) *
                                    pow(x / pdf->x[0], -lambda);
                        }
                }
        }
        else {
                /* We are within the table bounds. Let's locate the bracketing
                 * table rows.
                 */
                int ix0 = 0, ix1 = pdf->nx - 1, iQ0 = 0, iQ1 = pdf->nQ2 -1;
                table_bracket(pdf->x, x, &ix0, &ix1);
                table_bracket(pdf->Q2, Q2, &iQ0, &iQ1);
                const float hx = (x - pdf->x[ix0]) / (pdf->x[ix1] - pdf->x[ix0]);
                const float hQ = (Q2 - pdf->Q2[iQ0]) / (pdf->Q2[iQ1] - pdf->Q2[iQ0]);

                /* Interpolate the PDFs. */
                int i, nf = 2 * pdf->nf + 1;
                const float * const xfx00 = pdf->xfx + (ix0 * pdf->nQ2 + iQ0) * nf;
                const float * const xfx01 = pdf->xfx + (ix0 * pdf->nQ2 + iQ1) * nf;
                const float * const xfx10 = pdf->xfx + (ix1 * pdf->nQ2 + iQ0) * nf;
                const float * const xfx11 = pdf->xfx + (ix1 * pdf->nQ2 + iQ1) * nf;
                const float r00 = (1. - hx) * (1. - hQ);
                const float r01 = (1. - hx) * hQ;
                const float r10 = hx * (1. - hQ);
                const float r11 = hx * hQ;
                for (i = 0; i < nf; i++) xfx[i] = r00 * xfx00[i] + r01 * xfx01[i] + r10 * xfx10[i] + r11 * xfx11[i];
        }
}

/* DCS for Deep Inelastic Scattering (DIS) on nucleons. */
static double dcs_dis(struct ent_physics * physics,
    enum ent_projectile projectile, double energy, double Z, double A,
    enum ent_process process, double x, double y)
{
        /* Get the PDFs. */
        const double Q2 = 2. * x * y * ENT_MASS_NUCLEON * energy;
        const struct lha_pdf * const pdf = physics->pdf;
        float xfx[LHAPDF_NF_MAX];
        lha_pdf_compute(pdf, x, Q2, xfx);

        /* Compute the relevant structure functions. */
        const int eps = (projectile > 0) ? 1 : -1; /* CP? */
        const int nf = pdf->nf;
        const double N = A - Z;                    /* Number of neutrons. */
        double factor;
        if (process == ENT_PROCESS_DIS_CC) {
                /* Charged current DIS process. */
                const double d =
                    (Z <= 0.) ? 0. : xfx[1 * eps + nf];
                const double u =
                    (N <= 0.) ? 0. : xfx[2 * eps + nf];
                const double s = xfx[3 * eps + nf];
                const double b = xfx[5 * eps + nf];
                const double F1 = Z * d + N * u + A * (s + b);

                const double dbar =
                    (N <= 0.) ? 0. : xfx[-1 * eps + nf];
                const double ubar =
                    (Z <= 0.) ? 0. : xfx[-2 * eps + nf];
                const double cbar = xfx[-4 * eps + nf];
                const double F2 = N * dbar + Z * ubar + A * cbar;

                const double y1 = 1. - y;
                const double F = F1 + F2 * y1 * y1;
                const double MW2 = ENT_MASS_W * ENT_MASS_W;
                const double r = MW2 / (MW2 + Q2);
                factor = 2 * F * r * r;
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
                const double y1 = 1. - y;
                const double y12 = y1 * y1;
                const double F1 = (gp2 + gm2 * y12) * (Z * u + N * d + A * c);
                const double F2 =
                    (gpp2 + gpm2 * y12) * (Z * d + N * u + A * (s + b));

                const double dbar = xfx[-1 * eps + nf];
                const double ubar = xfx[-2 * eps + nf];
                const double sbar = xfx[-3 * eps + nf];
                const double cbar = xfx[-4 * eps + nf];
                const double bbar = xfx[-5 * eps + nf];
                const double F3 =
                    (gm2 + gp2 * y12) * (Z * ubar + N * dbar + A * cbar);
                const double F4 = (gpm2 + gpp2 * y12) *
                    (Z * dbar + N * ubar + A * (sbar + bbar));

                const double F = F1 + F2 + F3 + F4;
                const double MZ2 = ENT_MASS_Z * ENT_MASS_Z;
                const double r = MZ2 / (MZ2 + Q2);
                factor = 0.5 * F * r * r;
        }

        /* Return the DCS. */
        return energy * factor * (ENT_PHYS_GF * ENT_PHYS_HBC) *
            (ENT_PHYS_GF * ENT_PHYS_HBC) * ENT_MASS_NUCLEON / M_PI;
}

/* DCS for elastic scattering on electrons. */
static double dcs_elastic(
    enum ent_projectile projectile, double energy, double Z, double y)
{
        const double R = 2. * ENT_PHYS_SIN_THETA_W_2;
        const double L = 2. * ENT_PHYS_SIN_THETA_W_2 - 1.;
        const double MZ2 = ENT_MASS_Z * ENT_MASS_Z;
        const double rZ = MZ2 / (MZ2 + 2. * ENT_MASS_ELECTRON * energy * y);

        double factor;
        if (projectile == ENT_PROJECTILE_NU_E_BAR) {
                const double tmp1 = R * rZ;
                const double a =
                    ENT_MASS_W * ENT_MASS_W - 2. * ENT_MASS_ELECTRON * energy;
                const double b = ENT_WIDTH_W * ENT_MASS_W;
                const double c = 2. * ENT_MASS_W * ENT_MASS_W / (a * a + b * b);
                const double d = rZ + c * a;
                const double e = -c * b;
                factor = tmp1 * tmp1 + (d * d + e * e) * (1. - y) * (1. - y);
        } else if (projectile == ENT_PROJECTILE_NU_E) {
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
static double dcs_inverse(enum ent_projectile projectile, double energy,
    enum ent_process process, double Z, double y)
{
        if ((projectile != ENT_PROJECTILE_NU_E_BAR) &&
            ((projectile != ENT_PROJECTILE_NU_MU) ||
                (process != ENT_PROCESS_INVERSE_MUON)) &&
            ((projectile != ENT_PROJECTILE_NU_TAU) ||
                (process != ENT_PROCESS_INVERSE_TAU)))
                return 0.;

        const double ml = (process == ENT_PROCESS_INVERSE_MUON) ?
            ENT_MASS_MUON :
            ENT_MASS_TAU;
        const double MW2 = ENT_MASS_W * ENT_MASS_W;

        double factor;
        if (projectile == ENT_PROJECTILE_NU_E_BAR) {
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
static double dcs_glashow(enum ent_projectile projectile, double energy,
    enum ent_process process, double Z, double y)
{
        if (projectile != ENT_PROJECTILE_NU_E_BAR) return 0.;

        return (ENT_WIDTH_W / ENT_WIDTH_W_TO_MUON - 2.) *
            dcs_inverse(projectile, energy, ENT_PROCESS_INVERSE_MUON, Z, y);
}

/* Compute the total cross-section for DIS processes using a Gaussian
 * quadrature.
 */
static double dcs_dis_integrate(struct ent_physics * physics,
    enum ent_projectile projectile, double energy, double Z, double A,
    enum ent_process process)
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

        const double ymin = 1E-12;
        const double xmin = 1E-12;
        const int nx = 3;
        const int ny = 3;
        const double dlx = -log(xmin) / nx;
        const double dly = -log(ymin) / ny;
        double I = 0.;
        int i, j;
        for (i = 0; i < ny * N_GQ; i++) {
                const double y = ymin * exp((0.5 + 0.5 * xGQ[i % N_GQ] + i / N_GQ) * dly);
                double J = 0.;
                for (j = 0; j < nx * N_GQ; j++) {
                        const double x = xmin * exp((0.5 + 0.5 * xGQ[j % N_GQ] + j / N_GQ) * dlx);
                        J += wGQ[j % N_GQ] * x * dcs_dis(physics, projectile, energy,
                            Z, A, process, x, y);
                }
                I += wGQ[i % N_GQ] * y * J;
        }
        return 0.25 * I * dlx * dly;
#undef  N_GQ
}

/* Tabulate the interactions lengths and processes weights. */
static void physics_tabulate(struct ent_physics * physics)
{
        const enum ent_projectile projectile = ENT_PROJECTILE_NU_TAU;
        const enum ent_process process = ENT_PROCESS_DIS_CC;
        const double energy = 1E+12;
        const double Z = 0.5;
        const double A = 1.;

        const double I = dcs_dis_integrate(physics, projectile,
            energy, Z, A, process);
        printf("I = %.5lE\n", I);
}

/* API constructor for a Physics object. */
enum ent_return ent_physics_create(
    struct ent_physics ** physics, const char * pdf_file)
{
        ENT_ACKNOWLEDGE(ent_physics_create);

        enum ent_return rc;
        *physics = NULL;

        /* Load the PDFs. */
        FILE * stream;
        rc = ENT_RETURN_PATH_ERROR;
        if ((stream = fopen(pdf_file, "r")) == NULL) goto exit;

        if ((rc = lha_load(stream, physics)) != ENT_RETURN_SUCCESS) goto exit;
        rc = ENT_RETURN_SUCCESS;

        /* Build the tabulations. */
        physics_tabulate(*physics);

exit:
        if (stream != NULL) fclose(stream);
        ENT_RETURN(rc);
}

/* API destructor for a Physics object. */
void ent_physics_destroy(struct ent_physics ** physics)
{
        if ((physics == NULL) || (*physics == NULL)) return;

        free(*physics);
        *physics = NULL;
}

/* Generic API function for computing DCSs. */
enum ent_return ent_physics_dcs(struct ent_physics * physics,
    enum ent_projectile projectile, double energy, double Z, double A,
    enum ent_process process, double x, double y, double * dcs)
{
        ENT_ACKNOWLEDGE(ent_physics_dcs);

        /* Check the inputs. */
        *dcs = 0.;
        if ((x > 1.) || (x < 0.) || (y > 1.) || (y < 0.))
                ENT_RETURN(ENT_RETURN_DOMAIN_ERROR);

        /* Compute the corresponding DCS. */
        if (process == ENT_PROCESS_ELASTIC)
                *dcs = dcs_elastic(projectile, energy, Z, y);
        else if ((process == ENT_PROCESS_DIS_CC) ||
            (process == ENT_PROCESS_DIS_NC))
                *dcs =
                    dcs_dis(physics, projectile, energy, Z, A, process, x, y);
        else if ((process == ENT_PROCESS_INVERSE_MUON) ||
            (process == ENT_PROCESS_INVERSE_TAU))
                *dcs = dcs_inverse(projectile, energy, process, Z, y);
        else if (process == ENT_PROCESS_GLASHOW_HADRON)
                *dcs = dcs_glashow(projectile, energy, process, Z, y);
        else
                ENT_RETURN(ENT_RETURN_DOMAIN_ERROR);

        return ENT_RETURN_SUCCESS;
}

/* API interface to PDFs. */
enum ent_return ent_physics_pdf(struct ent_physics * physics,
    enum ent_parton parton, double x, double Q2, double * value)
{
        ENT_ACKNOWLEDGE(ent_physics_pdf);

        /* Check the inputs. */
        *value = 0.;
        const struct lha_pdf * const pdf = physics->pdf;
        if ((x > 1.) || (x <= 0.) || (Q2 < 0.) || (abs(parton) > pdf->nf))
                ENT_RETURN(ENT_RETURN_DOMAIN_ERROR);

        /* Compute the PDF and return. */
        float xfx[LHAPDF_NF_MAX];
        lha_pdf_compute(pdf, x, Q2, xfx);
        *value = xfx[pdf->nf + parton] / x;
        return ENT_RETURN_SUCCESS;
}

#if (GDB_MODE)
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
        /* Save the floating points exceptions status and enable them. */
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
