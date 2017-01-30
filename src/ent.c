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
/* The W boson partial width to muon + nu_mu, in GeV/c^2. */
#define ENT_WIDTH_W_TO_MUON 0.22164
/* The Fermi coupling constant GF/(hbar*c)^3, in GeV^-2. */
#define ENT_PHYS_GF 1.1663787E-05
/* The Planck constant as hbar*c, in GeV * m. */
#define ENT_PHYS_HBC 1.97326978E-16
/* The Weinberg angle at MZ, as sin(theta_W)^2. */
#define ENT_PHYS_SIN_THETA_W_2 0.231295
/* Avogadro's number, in mol^-1. */
#define ENT_PHYS_NA 6.022E+23

/* Energy range for tabulations. */
#define ENERGY_MIN 1E+03
#define ENERGY_MAX 1E+12
#define ENERGY_N 181

/* Sampling for the tabulations of DIS cumulative cross-sections. */
#define DIS_Q2_N 100
#define DIS_X_N 100

#ifndef M_PI
/* Define pi, if unknown. */
#define M_PI 3.14159265358979323846
#endif

/* Indices for Physics processes with a specific projectile and target. */
enum proget_index {
        /* Charged current DIS of a neutrino on a neutron. */
        PROGET_CC_NU_NEUTRON = 0,
        /* Neutral current DIS of a neutrino on a neutron. */
        PROGET_NC_NU_NEUTRON,
        /* Charged current DIS of an anti-neutrino on a neutron. */
        PROGET_CC_NU_BAR_NEUTRON,
        /* Neutral current DIS of an anti-neutrino on a neutron. */
        PROGET_NC_NU_BAR_NEUTRON,
        /* Charged current DIS of a neutrino on a proton. */
        PROGET_CC_NU_PROTON,
        /* Neutral current DIS of a neutrino on a proton. */
        PROGET_NC_NU_PROTON,
        /* Charged current DIS of an anti-neutrino on a proton. */
        PROGET_CC_NU_BAR_PROTON,
        /* Neutral current DIS of an anti-neutrino on a proton. */
        PROGET_NC_NU_BAR_PROTON,
        /* Elastic scattering of a nu_e on a an electron. */
        PROGET_ELASTIC_NU_E,
        /* Elastic scattering of an anti nu_e on a an electron. */
        PROGET_ELASTIC_NU_BAR_E,
        /* Elastic scattering of a nu_X={mu, tau} on a an electron. */
        PROGET_ELASTIC_NU_X,
        /* Elastic scattering of an anti nu_X={mu, tau} on a an electron. */
        PROGET_ELASTIC_NU_BAR_X,
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

/* Opaque Physics object. */
struct ent_physics {
        /* Index to the PDF data. */
        struct lha_pdf * pdf;
        /* Entry point for the total cross-sections. */
        double * cs;
        /* Entry point for the probability density function for DIS. */
        double * dis_pdf;
        /* Entry point for the cumulative density function for DIS. */
        double * dis_cdf;
        /* Sampling factors for x, in DIS. */
        double dis_xmin, dis_dlx, dis_rx;
        /* Sampling factors for Q2, in DIS. */
        double dis_Q2min, dis_dlQ2, dis_rQ2;
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
            ((PROGET_N - 1) * ENERGY_N + 8 * DIS_Q2_N * DIS_X_N +
                8 * (DIS_Q2_N - 1) * (DIS_X_N - 1)) *
                sizeof(double));
        if (v == NULL) return NULL;
        struct ent_physics * p = v;
        p->cs = (double *)(p->data + extra_size);
        p->dis_pdf = p->cs + (PROGET_N - 1) * ENERGY_N;
        p->dis_cdf = p->dis_pdf + 8 * DIS_Q2_N * DIS_X_N;
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
        REGISTER_FUNCTION(ent_physics_pdf);
        REGISTER_FUNCTION(ent_vertex);
        REGISTER_FUNCTION(ent_transport);

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
                        void * tmp = realloc(*buffer, sizeof(**buffer) +
                                (*buffer)->size + FILE_BLOCK_SIZE);
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
        float *x, *Q2, *lambda, *xfx;
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
                    (buffer->cursor[2] == '-'))
                        break;
        }
        long int pos = ftell(stream);

        /* Parse the table format. */
        file_get_line(&buffer, 0);
        const int nx = file_count_words(buffer);
        file_get_line(&buffer, 0);
        const int nQ = file_count_words(buffer);
        file_get_line(&buffer, 0);
        const int nf = file_count_words(buffer);
        if ((nf != LHAPDF_NF_MAX - 2) && (nf != LHAPDF_NF_MAX)) {
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
        unsigned int size =
            sizeof(struct lha_pdf) + size_x + size_Q + size_lambda + size_xfx;
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
                if (row[i] == 21.)
                        index[i] = pdf->nf;
                else
                        index[i] = (int)row[i] + pdf->nf;
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
static void table_bracket(const float * table, float value, int * p1, int * p2)
{
        int i3 = (*p1 + *p2) / 2;
        if (value >= table[i3])
                *p1 = i3;
        else
                *p2 = i3;
        if (*p2 - *p1 >= 2) table_bracket(table, value, p1, p2);
}

static void lha_pdf_compute(
    const struct lha_pdf * pdf, float x, float Q2, float * xfx)
{
        memset(xfx, 0x0, LHAPDF_NF_MAX * sizeof(*xfx));

        /* Check the bounds. */
        if ((Q2 < pdf->Q2[0]) || (Q2 >= pdf->Q2[pdf->nQ2 - 1]) ||
            (x >= pdf->x[pdf->nx - 1]))
                return;

        if (x < pdf->x[0]) {
                /* Extrapolate with a power law for small x values, i.e.
                 * x * f(x) ~ 1 / x**lambda(Q2).
                 */
                int iQ0, iQ1;
                float hQ;
                if (Q2 >= pdf->Q2[pdf->nQ2 - 1]) {
                        iQ0 = iQ1 = pdf->nQ2 - 1;
                        hQ = 1.;
                } else {
                        iQ0 = 0;
                        iQ1 = pdf->nQ2 - 1;
                        table_bracket(pdf->Q2, Q2, &iQ0, &iQ1);
                        hQ =
                            (Q2 - pdf->Q2[iQ0]) / (pdf->Q2[iQ1] - pdf->Q2[iQ0]);
                }

                int i, nf = 2 * pdf->nf + 1;
                const float * const lambda0 = pdf->lambda + iQ0 * nf;
                const float * const lambda1 = pdf->lambda + iQ1 * nf;
                const float * const xfx0 = pdf->xfx + iQ0 * nf;
                const float * const xfx1 = pdf->xfx + iQ1 * nf;
                for (i = 0; i < nf; i++) {
                        const float y0 = lambda0[i];
                        const float y1 = lambda1[i];
                        if ((y0 <= 0.) || (y1 <= 0.))
                                xfx[i] = 0.;
                        else {
                                const float lambda = y0 * (1. - hQ) + y1 * hQ;
                                xfx[i] = (xfx0[i] * (1. - hQ) + xfx1[i] * hQ) *
                                    pow(x / pdf->x[0], -lambda);
                        }
                }
        } else {
                /* We are within the table bounds. Let's locate the bracketing
                 * table rows.
                 */
                int ix0 = 0, ix1 = pdf->nx - 1, iQ0 = 0, iQ1 = pdf->nQ2 - 1;
                table_bracket(pdf->x, x, &ix0, &ix1);
                table_bracket(pdf->Q2, Q2, &iQ0, &iQ1);
                const float hx =
                    (x - pdf->x[ix0]) / (pdf->x[ix1] - pdf->x[ix0]);
                const float hQ =
                    (Q2 - pdf->Q2[iQ0]) / (pdf->Q2[iQ1] - pdf->Q2[iQ0]);

                /* Interpolate the PDFs. */
                int i, nf = 2 * pdf->nf + 1;
                const float * const xfx00 =
                    pdf->xfx + (ix0 * pdf->nQ2 + iQ0) * nf;
                const float * const xfx01 =
                    pdf->xfx + (ix0 * pdf->nQ2 + iQ1) * nf;
                const float * const xfx10 =
                    pdf->xfx + (ix1 * pdf->nQ2 + iQ0) * nf;
                const float * const xfx11 =
                    pdf->xfx + (ix1 * pdf->nQ2 + iQ1) * nf;
                const float r00 = (1. - hx) * (1. - hQ);
                const float r01 = (1. - hx) * hQ;
                const float r10 = hx * (1. - hQ);
                const float r11 = hx * hQ;
                for (i = 0; i < nf; i++)
                        xfx[i] = r00 * xfx00[i] + r01 * xfx01[i] +
                            r10 * xfx10[i] + r11 * xfx11[i];
        }
}

/* DCS for Deep Inelastic Scattering (DIS) on nucleons. */
static double dcs_dis(struct ent_physics * physics, enum ent_pid projectile,
    double energy, double Z, double A, enum ent_process process, double x,
    double y)
{
        /* Get the PDFs. */
        const double Q2 = 2. * x * y * ENT_MASS_NUCLEON * energy;
        const struct lha_pdf * const pdf = physics->pdf;
        float xfx[LHAPDF_NF_MAX];
        lha_pdf_compute(pdf, x, Q2, xfx);

        /* Compute the relevant structure functions. */
        const int eps = (projectile > 0) ? 1 : -1; /* CP? */
        const int nf = pdf->nf;
        const double N = A - Z; /* Number of neutrons. */
        double factor;
        if (process == ENT_PROCESS_DIS_CC) {
                /* Charged current DIS process. */
                const double d = (Z <= 0.) ? 0. : xfx[1 * eps + nf];
                const double u = (N <= 0.) ? 0. : xfx[2 * eps + nf];
                const double s = xfx[3 * eps + nf];
                const double b = xfx[5 * eps + nf];
                const double F1 = Z * d + N * u + A * (s + b);

                const double dbar = (N <= 0.) ? 0. : xfx[-1 * eps + nf];
                const double ubar = (Z <= 0.) ? 0. : xfx[-2 * eps + nf];
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

        return (ENT_WIDTH_W / ENT_WIDTH_W_TO_MUON - 2.) *
            dcs_inverse(projectile, energy, ENT_PROCESS_INVERSE_MUON, Z, y);
}

/* Compute the DCS for a given process, projectile and target. */
double dcs_compute(struct ent_physics * physics, enum ent_pid projectile,
    double energy, double Z, double A, enum ent_process process, double x,
    double y)
{
        if ((process == ENT_PROCESS_DIS_CC) || (process == ENT_PROCESS_DIS_NC))
                return dcs_dis(
                    physics, projectile, energy, Z, A, process, x, y);
        else if (process == ENT_PROCESS_ELASTIC)
                return dcs_elastic(projectile, energy, Z, y);
        else if ((process == ENT_PROCESS_INVERSE_MUON) ||
            (process == ENT_PROCESS_INVERSE_TAU))
                return dcs_inverse(projectile, energy, Z, process, y);
        else if (process == ENT_PROCESS_GLASHOW_HADRON)
                return dcs_glashow(projectile, energy, Z, y);
        return -DBL_MAX;
}

/* Compute the total cross-section for a processes using a Gaussian
 * quadrature.
 */
static double dcs_integrate(struct ent_physics * physics,
    enum ent_pid projectile, double energy, double Z, double A,
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

        if ((process == ENT_PROCESS_DIS_CC) ||
            (process == ENT_PROCESS_DIS_NC)) {
                /* Deep Inelastic Scattering requires a double integral. */
                const double ymin = 1E-12;
                const double xmin = 1E-12;
                const int nx = 3;
                const int ny = 3;
                const double dlx = -log(xmin) / nx;
                const double dly = -log(ymin) / ny;
                double I = 0.;
                int i, j;
                for (i = 0; i < ny * N_GQ; i++) {
                        const double y = ymin *
                            exp((0.5 + 0.5 * xGQ[i % N_GQ] + i / N_GQ) * dly);
                        double J = 0.;
                        for (j = 0; j < nx * N_GQ; j++) {
                                const double x = xmin *
                                    exp((0.5 + 0.5 * xGQ[j % N_GQ] + j / N_GQ) *
                                                     dlx);
                                J += wGQ[j % N_GQ] * x *
                                    dcs_compute(physics, projectile, energy, Z,
                                         A, process, x, y);
                        }
                        I += wGQ[i % N_GQ] * y * J;
                }
                return 0.25 * I * dlx * dly;
        } else {
                /* Scattering on an atomic electron. */
                const double ymin = 1E-06;
                const int ny = 6;
                const double dly = -log(ymin) / ny;
                double I = 0.;
                int i;
                for (i = 0; i < ny * N_GQ; i++) {
                        const double y = ymin *
                            exp((0.5 + 0.5 * xGQ[i % N_GQ] + i / N_GQ) * dly);
                        I += wGQ[i % N_GQ] * y * dcs_compute(physics,
                                                     projectile, energy, Z, A,
                                                     process, 0., y);
                }
                return 0.5 * I * dly;
        }
#undef N_GQ
}

/* Tabulate the interactions lengths and processes weights. */
static void physics_tabulate_cs(struct ent_physics * physics)
{
        const enum ent_process process[PROGET_N - 1] = { ENT_PROCESS_DIS_CC,
                ENT_PROCESS_DIS_NC, ENT_PROCESS_DIS_CC, ENT_PROCESS_DIS_NC,
                ENT_PROCESS_DIS_CC, ENT_PROCESS_DIS_NC, ENT_PROCESS_DIS_CC,
                ENT_PROCESS_DIS_NC, ENT_PROCESS_ELASTIC, ENT_PROCESS_ELASTIC,
                ENT_PROCESS_ELASTIC, ENT_PROCESS_ELASTIC,
                ENT_PROCESS_INVERSE_MUON, ENT_PROCESS_INVERSE_TAU,
                ENT_PROCESS_INVERSE_MUON, ENT_PROCESS_INVERSE_TAU };
        const enum ent_pid projectile[PROGET_N - 1] = { ENT_PID_NU_E,
                ENT_PID_NU_E, ENT_PID_NU_BAR_E, ENT_PID_NU_BAR_E, ENT_PID_NU_E,
                ENT_PID_NU_E, ENT_PID_NU_BAR_E, ENT_PID_NU_BAR_E, ENT_PID_NU_E,
                ENT_PID_NU_BAR_E, ENT_PID_NU_MU, ENT_PID_NU_BAR_MU,
                ENT_PID_NU_MU, ENT_PID_NU_TAU, ENT_PID_NU_BAR_E,
                ENT_PID_NU_BAR_E };

        const double dlE = log(ENERGY_MAX / ENERGY_MIN) / (ENERGY_N - 1);
        double * table;
        int i;
        for (i = 0, table = physics->cs; i < ENERGY_N;
             i++, table += PROGET_N - 1) {
                const double energy = ENERGY_MIN * exp(i * dlE);
                int j;
                for (j = 0; j < PROGET_N - 1; j++) {
                        const double Z =
                            (j <= PROGET_NC_NU_BAR_NEUTRON) ? 0. : 1.;
                        table[j] = dcs_integrate(
                            physics, projectile[j], energy, Z, 1., process[j]);
                }
        }
}

/* Tabulate the PDF and CDF for DIS events. */
static void physics_tabulate_dis(struct ent_physics * physics)
{
        /* Compute the sampling factors. */
        physics->dis_xmin = 1. / ENERGY_MAX;
        if (physics->dis_xmin > physics->pdf->x[0])
                physics->dis_xmin = physics->pdf->x[0];
        physics->dis_dlx =
            log(physics->pdf->x[physics->pdf->nx - 1] / physics->dis_xmin) /
            (DIS_X_N - 1);
        physics->dis_rx = exp(physics->dis_dlx);
        physics->dis_Q2min = physics->pdf->Q2[0];
        physics->dis_dlQ2 =
            log(physics->pdf->Q2[physics->pdf->nQ2 - 1] / physics->dis_Q2min) /
            (DIS_Q2_N - 1);
        physics->dis_rQ2 = exp(physics->dis_dlQ2);

        /* Compute the asymptotic PDF, i.e. y = 1, for the differential
         * cross-section as (ln(x), ln(Q2)). The corresponding PDF depends only
         * on (x, Q2) and it provides a majoration of the true DIS PDF.
         */
        const int nf = physics->pdf->nf;
        const double c =
            (ENT_PHYS_GF * ENT_PHYS_HBC) * (ENT_PHYS_GF * ENT_PHYS_HBC) / M_PI;
        const double MZ2 = ENT_MASS_Z * ENT_MASS_Z;
        const double MW2 = ENT_MASS_W * ENT_MASS_W;
        double * pdf = physics->dis_pdf;
        int iQ2;
        for (iQ2 = 0; iQ2 < DIS_Q2_N; iQ2++) {
                const double Q2 =
                    physics->dis_Q2min * exp(iQ2 * physics->dis_dlQ2);
                double rZ = MZ2 / (MZ2 + Q2);
                rZ *= rZ;
                double rW = MW2 / (MW2 + Q2);
                rW *= rW;
                int ix;
                for (ix = 0; ix < DIS_X_N; ix++, pdf++) {
                        const double x =
                            physics->dis_xmin * exp(ix * physics->dis_dlx);
                        float xfx[LHAPDF_NF_MAX];
                        lha_pdf_compute(physics->pdf, x, Q2, xfx);
                        int k;
                        double factor;
                        for (k = 0; k < 8; k++) {
                                const int eps = (k / 2) % 2 ? -1 : 1;
                                const double Z = (k / 4);
                                const double N = 1. - Z;
                                /* Compute the relevant structure functions. */
                                if ((k % 2) == 0) {
                                        /* Charged current DIS process. */
                                        const double d =
                                            (Z <= 0.) ? 0. : xfx[1 * eps + nf];
                                        const double u =
                                            (N <= 0.) ? 0. : xfx[2 * eps + nf];
                                        const double s = xfx[3 * eps + nf];
                                        const double b = xfx[5 * eps + nf];
                                        const double F1 = Z * d + N * u + s + b;

                                        const double dbar =
                                            (N <= 0.) ? 0. : xfx[-1 * eps + nf];
                                        const double ubar =
                                            (Z <= 0.) ? 0. : xfx[-2 * eps + nf];
                                        const double cbar = xfx[-4 * eps + nf];
                                        const double F2 =
                                            N * dbar + Z * ubar + cbar;

                                        factor = (F1 + F2) * rW;
                                } else {
                                        /* Neutral current DIS process. */
                                        const double d = xfx[1 * eps + nf];
                                        const double u = xfx[2 * eps + nf];
                                        const double s = xfx[3 * eps + nf];
                                        const double c = xfx[4 * eps + nf];
                                        const double b = xfx[5 * eps + nf];

                                        const double s23 =
                                            2. * ENT_PHYS_SIN_THETA_W_2 / 3.;
                                        const double gp2 =
                                            1. + 4. * s23 * (s23 - 1.);
                                        const double gm2 = 4. * s23 * s23;
                                        const double gpp2 =
                                            1. + s23 * (s23 - 2.);
                                        const double gpm2 = s23 * s23;
                                        const double F1 =
                                            (gp2 + gm2) * (Z * u + N * d + c);
                                        const double F2 = (gpp2 + gpm2) *
                                            (Z * d + N * u + s + b);

                                        const double dbar = xfx[-1 * eps + nf];
                                        const double ubar = xfx[-2 * eps + nf];
                                        const double sbar = xfx[-3 * eps + nf];
                                        const double cbar = xfx[-4 * eps + nf];
                                        const double bbar = xfx[-5 * eps + nf];
                                        const double F3 = (gm2 + gp2) *
                                            (Z * ubar + N * dbar + cbar);
                                        const double F4 = (gpm2 + gpp2) *
                                            (Z * dbar + N * ubar + sbar + bbar);

                                        factor =
                                            0.25 * (F1 + F2 + F3 + F4) * rZ;
                                }
                                pdf[k * DIS_X_N * DIS_Q2_N] = c * factor * Q2;
                        }
                }
        }

        const double c2 = 0.25 * physics->dis_dlQ2 * physics->dis_dlx;
        int proget;
        for (proget = 0; proget < 8; proget++) {
                /* Compute the corresponding CDF using a bilinear
                 * interpolation.
                 */
                pdf = physics->dis_pdf + proget * DIS_Q2_N * DIS_X_N;
                double * cdf =
                    physics->dis_cdf + proget * (DIS_Q2_N - 1) * (DIS_X_N - 1);
                double cdf0 = 0.;
                for (iQ2 = 0; iQ2 < DIS_Q2_N - 1; iQ2++) {
                        int ix;
                        for (ix = 0; ix < DIS_X_N - 1; ix++, cdf++) {
                                const double f00 = pdf[iQ2 * DIS_X_N + ix];
                                const double f01 = pdf[iQ2 * DIS_X_N + ix + 1];
                                const double f10 =
                                    pdf[(iQ2 + 1) * DIS_X_N + ix];
                                const double f11 =
                                    pdf[(iQ2 + 1) * DIS_X_N + ix + 1];
                                cdf0 += c2 * (f00 + f01 + f10 + f11);
                                *cdf = cdf0;
                        }
                }
        }
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
        physics_tabulate_cs(*physics);
        physics_tabulate_dis(*physics);

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
                if (target == ENT_PID_NEUTRON)
                        *proget = (projectile > 0) ? PROGET_CC_NU_NEUTRON :
                                                     PROGET_CC_NU_BAR_NEUTRON;
                else if (target == ENT_PID_PROTON)
                        *proget = (projectile > 0) ? PROGET_CC_NU_PROTON :
                                                     PROGET_CC_NU_BAR_PROTON;
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
                        else if ((projectile == ENT_PID_NU_MU) ||
                            (projectile == ENT_PID_NU_TAU))
                                *proget = PROGET_ELASTIC_NU_X;
                        else
                                *proget = PROGET_ELASTIC_NU_BAR_X;
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

/* Build the interpolation or extrapolation factors. */
static int cross_section_prepare(struct ent_physics * physics, double energy,
    double ** cs0, double ** cs1, double * p1, double * p2)
{
        int mode;
        if (energy < ENERGY_MIN) {
                /* Log extrapolation model below Emin. */
                mode = 1;
                *cs0 = physics->cs;
                *cs1 = physics->cs + PROGET_N - 1;
                *p1 = (ENERGY_N - 1) / log(ENERGY_MAX / ENERGY_MIN);
                *p2 = energy / ENERGY_MIN;
        } else if (energy >= ENERGY_MAX) {
                /* Log extrapolation model above Emax. */
                mode = 1;
                *cs0 = physics->cs + (ENERGY_N - 2) * (PROGET_N - 1);
                *cs1 = physics->cs + (ENERGY_N - 1) * (PROGET_N - 1);
                *p1 = (ENERGY_N - 1) / log(ENERGY_MAX / ENERGY_MIN);
                *p2 = energy / ENERGY_MIN;
        } else {
                /* Interpolation model. */
                mode = 0;
                const double dle =
                    log(ENERGY_MAX / ENERGY_MIN) / (ENERGY_N - 1);
                *p1 = log(energy / ENERGY_MIN) / dle;
                const int i0 = (int)(*p1);
                *p1 -= i0;
                *cs0 = physics->cs + i0 * (PROGET_N - 1);
                *cs1 = physics->cs + (i0 + 1) * (PROGET_N - 1);
                *p2 = 0.;
        }
        return mode;
}

/* Low level routine for computing a specific cross-section by interpolation
 * or extrapolation.
 */
static double cross_section_compute(
    int mode, int proget, double * cs0, double * cs1, double p1, double p2)
{
        if (mode == 0) {
                /* interpolation case. */
                return cs0[proget] * (1. - p1) + cs1[proget] * p1;
        } else {
                /* Extrapolation case. */
                const double a = log(cs1[proget] / cs0[proget]) * p1;
                return cs0[proget] * pow(p2, a);
        }
}

/* Generic API function for accessing a specific total cross-section. */
enum ent_return ent_physics_cross_section(struct ent_physics * physics,
    enum ent_pid projectile, double energy, double Z, double A,
    enum ent_process process, double * cross_section)
{
        ENT_ACKNOWLEDGE(ent_physics_cross_section);
        *cross_section = 0.;

        /* Build the interpolation or extrapolation factors. */
        double *cs0, *cs1;
        double p1, p2;
        int mode = cross_section_prepare(physics, energy, &cs0, &cs1, &p1, &p2);

        /* Compute the relevant process-target indices. */
        enum ent_return rc;
        if ((process == ENT_PROCESS_DIS_CC) ||
            (process == ENT_PROCESS_DIS_NC)) {
                int proget;
                if (Z > 0.) {
                        if ((rc = proget_compute(projectile, ENT_PID_PROTON,
                                 process, &proget)) != ENT_RETURN_SUCCESS)
                                ENT_RETURN(rc);
                        *cross_section += Z * cross_section_compute(mode,
                                                  proget, cs0, cs1, p1, p2);
                }
                const double N = A - Z;
                if (N > 0.) {
                        if ((rc = proget_compute(projectile, ENT_PID_NEUTRON,
                                 process, &proget)) != ENT_RETURN_SUCCESS)
                                ENT_RETURN(rc);
                        *cross_section += N * cross_section_compute(mode,
                                                  proget, cs0, cs1, p1, p2);
                }
        } else {
                int proget;
                if ((rc = proget_compute(projectile, ENT_PID_ELECTRON, process,
                         &proget)) != ENT_RETURN_SUCCESS)
                        ENT_RETURN(rc);
                if (proget < PROGET_N - 1) {
                        *cross_section = Z * cross_section_compute(mode, proget,
                                                 cs0, cs1, p1, p2);
                } else {
                        *cross_section = Z *
                            cross_section_compute(mode, PROGET_N - 2, cs0, cs1,
                                             p1, p2) *
                            (ENT_WIDTH_W / ENT_WIDTH_W_TO_MUON - 1.);
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
        const double sgn = context->forward ? 1. : -1.;
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
        double step1 = context->medium(state, &medium1);
        double ds_boundary = 0.;
        if (*medium == NULL) {
                /* initialisation step. */
                *medium = medium1;
                if (medium1 == NULL) {
                        *event = ENT_EVENT_EXIT;
                        return ENT_RETURN_SUCCESS;
                }
        } else if ((*medium != NULL) && (medium1 != *medium)) {
                /* A change of medium occured. Let's locate the interface
                 * by dichotomy.
                 */
                double pi[3];
                memcpy(pi, state->position, sizeof(pi));
                double s1 = 0., s2 = -(*step);
                while (fabs(s1 - s2) > STEP_MIN) {
                        double s3 = 0.5 * (s1 + s2);
                        state->position[0] =
                            pi[0] + s3 * sgn * state->direction[0];
                        state->position[1] =
                            pi[1] + s3 * sgn * state->direction[1];
                        state->position[2] =
                            pi[2] + s3 * sgn * state->direction[2];
                        step1 = context->medium(state, &medium1);
                        if (medium1 == *medium)
                                s1 = s3;
                        else
                                s2 = s3;
                }
                state->position[0] = pi[0] + s2 * sgn * state->direction[0];
                state->position[1] = pi[1] + s2 * sgn * state->direction[1];
                state->position[2] = pi[2] + s2 * sgn * state->direction[2];
                *step += s1;
                ds_boundary = s1 - s2;
        }

        /* Get the end step local properties. */
        double density1;
        const double step2 = (*medium)->density(*medium, state, &density1);
        if ((density1 < 0.) || ((*medium)->A <= 0.))
                return ENT_RETURN_DOMAIN_ERROR;
        if (step2 > 0.) {
                if ((step1 <= 0.) || (step2 < step1)) step1 = step2;
        }

        /* Offset the end step position if a boundary crossing occured. */
        if (ds_boundary != 0.) {
                state->position[0] += ds_boundary * sgn * state->direction[0];
                state->position[1] += ds_boundary * sgn * state->direction[1];
                state->position[2] += ds_boundary * sgn * state->direction[2];
                *medium = medium1;
                *step += ds_boundary;

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
                                ds = *step * ((sqrt(density1 * density1 +
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

        /* Update and return. */
        *step = (step1 < 0) ? 0. : step1;
        *density = density1;
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
        const double sgn = context->forward ? 1. : -1.;
        state->position[0] += sgn * state->direction[0] * ds;
        state->position[1] += sgn * state->direction[1] * ds;
        state->position[2] += sgn * state->direction[2] * ds;

        return event;
}

/* Sample the x and Q2 parameters in a DIS event using a bilinear interpolation
 * of the PDF and rejection sampling.
 */
static enum ent_return transport_sample_yQ2(struct ent_physics * physics,
    struct ent_context * context, struct ent_state * neutrino, int proget,
    double * y_, double * Q2_)
{
/* Energy threshold between brute force rejection sampling or a
 * matched enveloppe. Caution : changing xmin has a strong impact on this
 * optimisation.
 */
#define DIS_BRUTE_FORCE_THRESHOLD 1E+08

        /* Unpack the tables, ect ... */
        const double * const pdf =
            physics->dis_pdf + proget * DIS_Q2_N * DIS_X_N;
        const double * const cdf =
            physics->dis_cdf + proget * (DIS_Q2_N - 1) * (DIS_X_N - 1);
        const double Q2max = 2. * ENT_MASS_NUCLEON * neutrino->energy;
        int iQ2max = ceil(log(Q2max / physics->dis_Q2min) / physics->dis_dlQ2);
        if (iQ2max >= DIS_Q2_N)
                iQ2max = DIS_Q2_N - 1;
        else if (iQ2max == 0)
                iQ2max = 1;
        else if (iQ2max < 0) {
                /* This should not occur since the total cross-section would
                 * be null in this case.
                 */
                return ENT_RETURN_DOMAIN_ERROR;
        }
        const double ly0 =
            log(physics->dis_Q2min / (physics->dis_xmin * Q2max)) -
            physics->dis_dlx;

        double y = 0., Q2 = 0.;
        for (;;) {
                double x = 0.;

                /* Sample (x, Q2) over a bilinear interpolation of the
                 * asymptotic pdf in (ln(x), ln(Q2)).
                 */
                int ix, iQ2;
                double r;
                if (neutrino->energy >= DIS_BRUTE_FORCE_THRESHOLD) {
                        /* Map the table's cell index (ix, iQ2) using a
                         * dichotomy and rejection sampling. */
                        int i1 = iQ2max * (DIS_X_N - 1) - 1;
                        r = cdf[i1] * context->random(context);
                        if (r <= cdf[0]) {
                                i1 = ix = iQ2 = 0;
                        } else {
                                int i0 = 0;
                                while (i1 - i0 > 1) {
                                        int i2 = (i0 + i1) / 2;
                                        if (r > cdf[i2])
                                                i0 = i2;
                                        else
                                                i1 = i2;
                                }
                                ix = i1 % (DIS_X_N - 1);
                                iQ2 = i1 / (DIS_X_N - 1);
                                const double ly = ly0 +
                                    iQ2 * physics->dis_dlQ2 -
                                    ix * physics->dis_dlx;
                                if (ly >= 0.) continue;
                                r -= cdf[i0];
                        }
                } else {
                        /* At low energies it becomes more efficient to spend
                         * some time for computing a matched enveloppe rather
                         * than relying on brute force rejection sampling.
                         */
                        const double a = physics->dis_dlQ2 / physics->dis_dlx;
                        const double b = ly0 / physics->dis_dlx;
                        double d = 0.;
                        int i;
                        for (i = 0; i <= iQ2max; i++) {
                                int jmin = floor(a * i + b);
                                if (jmin < 0)
                                        jmin = 0;
                                else if (jmin > DIS_X_N - 2)
                                        jmin = DIS_X_N - 2;
                                int j;
                                for (j = jmin; j < DIS_X_N - 1; j++) {
                                        const int k = i * (DIS_X_N - 1) + j;
                                        const double delta = (k == 0) ?
                                            cdf[k] :
                                            cdf[k] - cdf[k - 1];
                                        d += delta;
                                }
                        }

                        iQ2 = ix = 0;
                        r = d * context->random(context);
                        d = 0.;
                        for (i = 0; i <= iQ2max; i++) {
                                int jmin = floor(a * i + b);
                                if (jmin < 0)
                                        jmin = 0;
                                else if (jmin > DIS_X_N - 2)
                                        jmin = DIS_X_N - 2;
                                int j;
                                for (j = jmin; j < DIS_X_N - 1; j++) {
                                        const int k = i * (DIS_X_N - 1) + j;
                                        const double delta = (k == 0) ?
                                            cdf[k] :
                                            cdf[k] - cdf[k - 1];
                                        if (r <= d + delta) {
                                                ix = j;
                                                iQ2 = i;
                                                goto cell_found;
                                        }
                                        d += delta;
                                }
                        }
                cell_found:
                        r -= d;
                }

                /* Sample over the cell. First sample over the marginal
                 * CDF for x, given the bilinear model.
                 */
                const double f00 = pdf[iQ2 * DIS_X_N + ix];
                const double f01 = pdf[iQ2 * DIS_X_N + ix + 1];
                const double f10 = pdf[(iQ2 + 1) * DIS_X_N + ix];
                const double f11 = pdf[(iQ2 + 1) * DIS_X_N + ix + 1];
                const double ax = f10 + f11 - f00 - f01;
                const double bx = f00 + f01;
                const double cx =
                    -4. * r / (physics->dis_dlQ2 * physics->dis_dlx);
                double dx = bx * bx - ax * cx;
                if (dx < 0.)
                        dx = 0.;
                else
                        dx = sqrt(dx);
                const double hx = (dx - bx) / ax;
                x = physics->dis_xmin * exp((ix + hx) * physics->dis_dlx);
                if (x <= 0.) continue;

                /* Then sample over the conditional CDF for Q2. */
                const double f0 = f00 * (1. - hx) + f10 * hx;
                const double f1 = f01 * (1. - hx) + f11 * hx;
                const double ay = f1 - f0;
                const double by = f0;
                const double cy = -context->random(context) * (f0 + f1);
                double dy = by * by - ay * cy;
                if (dy < 0.)
                        dy = 0.;
                else
                        dy = sqrt(dy);
                const double hy = (dy - by) / ay;
                Q2 = physics->dis_Q2min * exp((iQ2 + hy) * physics->dis_dlQ2);
                if (Q2 >= Q2max) continue;

                /* Reject the sampled value if it violates the kinematic
                 * limit.
                 */
                y = Q2 / (x * Q2max);
                if (y > 1.) continue;

                /* Reject according to the true PDF, i.e. including the
                 * (1 - y)**2 factor(s). First let us compute the relevant
                 * structure functions.
                 */
                const int nf = physics->pdf->nf;
                const int eps = (proget / 2) % 2 ? -1 : 1;
                const double Z = (proget / 4);
                const double N = 1. - Z;
                float xfx[LHAPDF_NF_MAX];
                lha_pdf_compute(physics->pdf, x, Q2, xfx);

                double factor0, factor1;
                if ((proget % 2) == 0) {
                        /* Charged current DIS process. */
                        const double d = (Z <= 0.) ? 0. : xfx[1 * eps + nf];
                        const double u = (N <= 0.) ? 0. : xfx[2 * eps + nf];
                        const double s = xfx[3 * eps + nf];
                        const double b = xfx[5 * eps + nf];
                        const double F1 = Z * d + N * u + s + b;

                        const double dbar = (N <= 0.) ? 0. : xfx[-1 * eps + nf];
                        const double ubar = (Z <= 0.) ? 0. : xfx[-2 * eps + nf];
                        const double cbar = xfx[-4 * eps + nf];
                        const double F2 = N * dbar + Z * ubar + cbar;

                        factor0 = F1 + F2;
                        factor1 = F1 + F2 * (1. - y) * (1. - y);
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
                        factor0 = (gp2 + gm2) * (Z * u + N * d + c);
                        factor1 = (gp2 + gm2 * y12) * (Z * u + N * d + c);
                        factor0 += (gpp2 + gpm2) * (Z * d + N * u + s + b);
                        factor1 +=
                            (gpp2 + gpm2 * y12) * (Z * d + N * u + s + b);

                        const double dbar = xfx[-1 * eps + nf];
                        const double ubar = xfx[-2 * eps + nf];
                        const double sbar = xfx[-3 * eps + nf];
                        const double cbar = xfx[-4 * eps + nf];
                        const double bbar = xfx[-5 * eps + nf];
                        factor0 += (gm2 + gp2) * (Z * ubar + N * dbar + cbar);
                        factor1 +=
                            (gm2 + gp2 * y12) * (Z * ubar + N * dbar + cbar);
                        factor0 +=
                            (gpm2 + gpp2) * (Z * dbar + N * ubar + sbar + bbar);
                        factor1 += (gpm2 + gpp2 * y12) *
                            (Z * dbar + N * ubar + sbar + bbar);
                }

                /* Then, reject sample accordingly. */
                if (context->random(context) * factor0 <= factor1) break;
        }

        *y_ = y;
        *Q2_ = Q2;
        return ENT_RETURN_SUCCESS;
}

/* Sample the inelasticity, _y_, for an interaction with an electron. The
 * sampling is done with a combination of inverse and rejection sampling.
 */
static double transport_sample_y(struct ent_context * context,
    struct ent_state * state, int proget, enum ent_pid * ejectile,
    enum ent_pid * recoil, double * mu)
{
        const double yZ =
            0.5 * ENT_MASS_Z * ENT_MASS_Z / (ENT_MASS_ELECTRON * state->energy);
        const double yW =
            0.5 * ENT_MASS_W * ENT_MASS_W / (ENT_MASS_ELECTRON * state->energy);
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

        } else if (proget <= PROGET_ELASTIC_NU_BAR_X) {
                const double r1 = (proget == PROGET_ELASTIC_NU_X) ?
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
                const double u = context->random(context);
                y = (1. + yW) * u / (yW + u);
                *ejectile = ENT_PID_NU_E;
                if (proget == PROGET_INVERSE_NU_MU_MU) {
                        *recoil = ENT_PID_MUON;
                        *mu = ENT_MASS_MUON;
                } else {
                        *recoil = ENT_PID_TAU;
                        *mu = ENT_MASS_TAU;
                }
        } else {
                const double u = context->random(context);
                y = 1. - pow(u, 1. / 3.);
                if (proget == PROGET_INVERSE_NU_BAR_E_MU) {
                        *ejectile = ENT_PID_NU_BAR_MU;
                        *recoil = ENT_PID_MUON;
                        *mu = ENT_MASS_MUON;
                } else {
                        *ejectile = ENT_PID_NU_BAR_TAU;
                        *recoil = ENT_PID_TAU;
                        *mu = ENT_MASS_TAU;
                }
        }

        return y;
}

/* Rotate the state direction. */
static void transport_rotate(
    struct ent_state * state, double cos_theta, double cos_phi, double sin_phi)
{
        /* Unpack the direction. */
        double * const direction = state->direction;

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
                        const double nrm =
                            1. / sqrt(direction[0] * direction[0] +
                                     direction[2] * direction[2]);
                        u0x = -direction[2] * nrm, u0z = direction[0] * nrm;
                } else {
                        const double nrm =
                            1. / sqrt(direction[1] * direction[1] +
                                     direction[2] * direction[2]);
                        u0y = direction[2] * nrm, u0z = -direction[1] * nrm;
                }
        } else {
                if (a1 > a2) {
                        const double nrm =
                            1. / sqrt(direction[0] * direction[0] +
                                     direction[1] * direction[1]);
                        u0x = direction[1] * nrm, u0y = -direction[0] * nrm;
                } else {
                        const double nrm =
                            1. / sqrt(direction[1] * direction[1] +
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

/* Process an interaction vertex. */
static enum ent_return transport_vertex(struct ent_physics * physics,
    struct ent_context * context, int proget, struct ent_state * neutrino,
    struct ent_state * product)
{
        /* Process the corresponding vertex. */
        if (proget <= PROGET_NC_NU_BAR_PROTON) {
                double y, Q2;
                enum ent_return rc = transport_sample_yQ2(
                    physics, context, neutrino, proget, &y, &Q2);
                if (rc != ENT_RETURN_SUCCESS) return rc;
                struct ent_state product_;
                memcpy(&product_, neutrino, sizeof(product_));
                if ((proget % 2) == 0) {
                        /* Charged current event : compute the final states. */
                        double mu;
                        if (neutrino->pid == ENT_PID_NU_E) {
                                mu = ENT_MASS_ELECTRON;
                                product_.pid = ENT_PID_ELECTRON;
                        } else if (neutrino->pid == ENT_PID_NU_BAR_E) {
                                mu = ENT_MASS_ELECTRON;
                                product_.pid = ENT_PID_POSITRON;
                        } else if (neutrino->pid == ENT_PID_NU_MU) {
                                mu = ENT_MASS_MUON;
                                product_.pid = ENT_PID_MUON;
                        } else if (neutrino->pid == ENT_PID_NU_BAR_MU) {
                                mu = ENT_MASS_MUON;
                                product_.pid = ENT_PID_MUON_BAR;
                        } else if (neutrino->pid == ENT_PID_NU_TAU) {
                                mu = ENT_MASS_TAU;
                                product_.pid = ENT_PID_TAU;
                        } else if (neutrino->pid == ENT_PID_NU_BAR_TAU) {
                                mu = ENT_MASS_TAU;
                                product_.pid = ENT_PID_TAU_BAR;
                        } else {
                                return ENT_RETURN_DOMAIN_ERROR;
                        }

                        double Emu = neutrino->energy * (1. - y);
                        if (Emu > 0.) {
                                /* Compute the product lepton's energy and
                                 * direction.
                                 */
                                const double pmu = sqrt(Emu * (Emu + 2. * mu));
                                double ct = (Emu -
                                                0.5 * (Q2 + mu * mu) /
                                                    neutrino->energy) /
                                    pmu;
                                if (ct > 1.)
                                        ct = 1.;
                                else if (ct < -1.)
                                        ct = -1.;
                                const double phi =
                                    2. * M_PI * context->random(context);
                                const double cp = cos(phi);
                                const double sp = sin(phi);
                                transport_rotate(&product_, ct, cp, sp);

                                /* Compute the hadron's average direction
                                 * from momentum conservation.
                                 */
                                neutrino->direction[0] =
                                    neutrino->direction[0] * neutrino->energy -
                                    pmu * product_.direction[0];
                                neutrino->direction[1] =
                                    neutrino->direction[1] * neutrino->energy -
                                    pmu * product_.direction[1];
                                neutrino->direction[2] =
                                    neutrino->direction[2] * neutrino->energy -
                                    pmu * product_.direction[2];
                                double d = neutrino->direction[0] *
                                        neutrino->direction[0] +
                                    neutrino->direction[1] *
                                        neutrino->direction[1] +
                                    neutrino->direction[2] *
                                        neutrino->direction[2];
                                if (d > 0.) {
                                        d = 1. / sqrt(d);
                                        neutrino->direction[0] *= d;
                                        neutrino->direction[1] *= d;
                                        neutrino->direction[2] *= d;
                                }
                        } else
                                Emu = 0.;
                        product_.energy = Emu;
                        neutrino->energy *= y;
                        neutrino->pid = ENT_PID_HADRON;
                } else {
                        /* Neutral current event : compute the final states. */
                        product_.pid = ENT_PID_HADRON;
                        if (y < 1.) {
                                const double Ep = neutrino->energy * (1. - y);
                                double ct =
                                    1. - 0.5 * Q2 / (neutrino->energy * Ep);
                                if (ct < -1.) ct = -1.;
                                const double phi =
                                    2. * M_PI * context->random(context);
                                const double cp = cos(phi);
                                const double sp = sin(phi);
                                transport_rotate(neutrino, ct, cp, sp);

                                /* Compute the hadron's average direction
                                 * from momentum conservation.
                                 */
                                const double Eh = neutrino->energy * y;
                                product_.direction[0] =
                                    product_.direction[0] * neutrino->energy -
                                    Ep * neutrino->direction[0];
                                product_.direction[1] =
                                    product_.direction[1] * neutrino->energy -
                                    Ep * neutrino->direction[1];
                                product_.direction[2] =
                                    product_.direction[2] * neutrino->energy -
                                    Ep * neutrino->direction[2];
                                double d = product_.direction[0] *
                                        product_.direction[0] +
                                    product_.direction[1] *
                                        product_.direction[1] +
                                    product_.direction[2] *
                                        product_.direction[2];
                                if (d > 0.) {
                                        d = 1. / sqrt(d);
                                        product_.direction[0] *= d;
                                        product_.direction[1] *= d;
                                        product_.direction[2] *= d;
                                }
                                product_.energy = Eh;
                                neutrino->energy = Ep;
                        } else {
                                product_.energy = neutrino->energy;
                                neutrino->energy = 0.;
                        }
                }

                /* Copy back the product data if requested. */
                if (product != NULL)
                        memcpy(product, &product_, sizeof(product_));
        } else if (proget < PROGET_GLASHOW_HADRONS) {
                /* This is an interaction with an atomic electron and a
                 * neutrino in the final states.
                 */
                enum ent_pid ejectile = neutrino->pid;
                enum ent_pid recoil = ENT_PID_NONE;
                double y, ce, cr;
                for (;;) {
                        /* Let's first sample the inelasticity _y_. */
                        double mu;
                        y = transport_sample_y(
                            context, neutrino, proget, &ejectile, &recoil, &mu);
                        if (y < mu / neutrino->energy) continue;

                        /* Then, compute the cosines of the polar angles. */
                        const double pe =
                            neutrino->energy * (1. - y) + ENT_MASS_ELECTRON;
                        const double peL =
                            (neutrino->energy + ENT_MASS_ELECTRON) * (1. - y) +
                            0.5 * (ENT_MASS_ELECTRON * ENT_MASS_ELECTRON +
                                      mu * mu) /
                                neutrino->energy;
                        if (peL > pe) continue;
                        ce = peL / pe;
                        const double pr =
                            sqrt(y * y * neutrino->energy * neutrino->energy -
                                mu * mu);
                        const double prL = neutrino->energy - peL;
                        if (prL > pr) continue;
                        cr = prL / pr;
                        break;
                }

                /* Update the particles states. */
                const double phi = 2. * M_PI * context->random(context);
                const double cp = cos(phi);
                const double sp = sin(phi);
                if ((product != NULL) && (recoil != ENT_PID_NONE)) {
                        memcpy(product, neutrino, sizeof(*product));
                        product->pid = recoil;
                        product->energy = y * neutrino->energy;
                        transport_rotate(product, cr, -cp, -sp);
                }
                neutrino->pid = ejectile;
                neutrino->energy = neutrino->energy * (1. - y);
                transport_rotate(neutrino, ce, cp, sp);
        } else if (proget == PROGET_GLASHOW_HADRONS) {
                /* This is a total conversion of a anti nu_e neutrino on an
                 * atomic electron. */
                neutrino->pid = ENT_PID_HADRON;
        } else {
                return ENT_RETURN_DOMAIN_ERROR;
        }

        return ENT_RETURN_SUCCESS;
}

/* Compute the tranport cross-sections for a given projectile and
 * medium. */
static enum ent_return transport_cross_section(struct ent_physics * physics,
    enum ent_pid projectile, double energy, double Z, double A, double * cs)
{
        /* Build the interpolation or extrapolation factors. */
        double *cs0, *cs1;
        double p1, p2;
        int mode = cross_section_prepare(physics, energy, &cs0, &cs1, &p1, &p2);

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
        cs[0] = N0 * cross_section_compute(mode, 0, cs0, cs1, p1, p2);
        cs[1] = cs[0] + N0 * cross_section_compute(mode, 1, cs0, cs1, p1, p2);
        cs[2] = cs[1] + N1 * cross_section_compute(mode, 2, cs0, cs1, p1, p2);
        cs[3] = cs[2] + N1 * cross_section_compute(mode, 3, cs0, cs1, p1, p2);
        cs[4] = cs[3] + Z0 * cross_section_compute(mode, 4, cs0, cs1, p1, p2);
        cs[5] = cs[4] + Z0 * cross_section_compute(mode, 5, cs0, cs1, p1, p2);
        cs[6] = cs[5] + Z1 * cross_section_compute(mode, 6, cs0, cs1, p1, p2);
        cs[7] = cs[6] + Z1 * cross_section_compute(mode, 7, cs0, cs1, p1, p2);
        if (projectile == ENT_PID_NU_E)
                cs[8] = cs[7] +
                    Z * cross_section_compute(mode, 8, cs0, cs1, p1, p2);
        else
                cs[8] = cs[7];
        if (projectile == ENT_PID_NU_BAR_E)
                cs[9] = cs[8] +
                    Z * cross_section_compute(mode, 9, cs0, cs1, p1, p2);
        else
                cs[9] = cs[8];
        if ((projectile == ENT_PID_NU_MU) || (projectile == ENT_PID_NU_TAU))
                cs[10] = cs[9] +
                    Z * cross_section_compute(mode, 10, cs0, cs1, p1, p2);
        else
                cs[10] = cs[9];
        if ((projectile == ENT_PID_NU_BAR_MU) ||
            (projectile == ENT_PID_NU_BAR_TAU))
                cs[11] = cs[10] +
                    Z * cross_section_compute(mode, 11, cs0, cs1, p1, p2);
        else
                cs[11] = cs[10];
        if (projectile == ENT_PID_NU_MU)
                cs[12] = cs[11] +
                    Z * cross_section_compute(mode, 12, cs0, cs1, p1, p2);
        else
                cs[12] = cs[11];
        if (projectile == ENT_PID_NU_TAU)
                cs[13] =
                    cs[12] + cross_section_compute(mode, 13, cs0, cs1, p1, p2);
        else
                cs[13] = cs[12];
        if (projectile == ENT_PID_NU_BAR_E) {
                const double d =
                    Z * cross_section_compute(mode, 14, cs0, cs1, p1, p2);
                cs[14] = cs[13] + d;
                cs[15] = cs[14] +
                    Z * cross_section_compute(mode, 15, cs0, cs1, p1, p2);
                cs[16] = cs[15] + d * (ENT_WIDTH_W / ENT_WIDTH_W_TO_MUON - 1.);
        } else {
                cs[14] = cs[13];
                cs[15] = cs[13];
                cs[16] = cs[13];
        }

        /* Check the result and return. */
        if (cs[PROGET_N - 1] <= 0.) return ENT_RETURN_DOMAIN_ERROR;
        return ENT_RETURN_SUCCESS;
}

/* API interface for an interaction vertex. */
enum ent_return ent_vertex(struct ent_physics * physics,
    struct ent_context * context, enum ent_process process, enum ent_pid target,
    struct ent_state * neutrino, struct ent_state * product)
{
        ENT_ACKNOWLEDGE(ent_vertex);
        if (product != NULL) product->pid = ENT_PID_NONE;

        /* Check and format the inputs. */
        if ((physics == NULL) || (context == NULL) ||
            (context->random == NULL) || (neutrino == NULL))
                ENT_RETURN(ENT_RETURN_BAD_ADDRESS);

        /* Compute the process-target index. */
        int proget;
        enum ent_return rc;
        if ((rc = proget_compute(neutrino->pid, target, process, &proget)) !=
            ENT_RETURN_SUCCESS)
                ENT_RETURN(rc);

        /* Process the vertex. */
        ENT_RETURN(
            transport_vertex(physics, context, proget, neutrino, product));
}

/* API interface for the neutrino tranport. */
enum ent_return ent_transport(struct ent_physics * physics,
    struct ent_context * context, struct ent_state * neutrino,
    struct ent_state * product, enum ent_event * event)
{
        ENT_ACKNOWLEDGE(ent_transport);
        if (product != NULL) product->pid = ENT_PID_NONE;
        *event = ENT_EVENT_NONE;

        /* Check and format the inputs. */
        if ((physics == NULL) || (context == NULL) ||
            (context->medium == NULL) || (context->random == NULL) ||
            (neutrino == NULL))
                ENT_RETURN(ENT_RETURN_BAD_ADDRESS);

        /* Check for an initial limit violation. */
        enum ent_return rc = ENT_RETURN_SUCCESS;
        enum ent_event event_ = ENT_EVENT_NONE;
        if ((context->distance_max > 0.) &&
            (neutrino->distance >= context->distance_max)) {
                event_ = ENT_EVENT_LIMIT_DISTANCE;
                goto exit;
        }
        if ((context->grammage_max > 0.) &&
            (neutrino->grammage >= context->grammage_max)) {
                event_ = ENT_EVENT_LIMIT_GRAMMAGE;
                goto exit;
        }

        /* Initialise the transport. */
        double step = 0., density;
        struct ent_medium * medium = NULL;
        if ((rc = transport_step(context, neutrino, &medium, &step, &density,
                 context->grammage_max, &event_)) != ENT_RETURN_SUCCESS)
                goto exit;
        if (event_ != ENT_EVENT_NONE) goto exit;

        /* Compute the cross-sections. */
        double cs[PROGET_N];
        if ((rc = transport_cross_section(physics, neutrino->pid,
                 neutrino->energy, medium->Z, medium->A, cs)) !=
            ENT_RETURN_SUCCESS)
                goto exit;

        /* Randomise the depth of the next interaction. */
        enum ent_event foreseen = ENT_EVENT_NONE;
        const double Xint =
            medium->A * 1E-03 / (cs[PROGET_N - 1] * ENT_PHYS_NA);
        double Xlim = neutrino->grammage - Xint * log(context->random(context));
        if (context->grammage_max && (context->grammage_max < Xlim)) {
                Xlim = context->grammage_max;
                foreseen = ENT_EVENT_LIMIT_GRAMMAGE;
        }

        if (step)
                for (;;) {
                        /* Step until an event occurs. */
                        if ((rc = transport_step(context, neutrino, &medium,
                                 &step, &density, Xlim, &event_)) !=
                            ENT_RETURN_SUCCESS)
                                goto exit;
                        if (event_ != ENT_EVENT_NONE) break;
                }
        else {
                /* This is a uniform medium of infinite
                 * extension.
                 * Let's do a single straight step.
                 */
                event_ = transport_straight(context, neutrino, density, Xlim);
        }

        /* Process any interaction. */
        if ((event_ == ENT_EVENT_LIMIT_GRAMMAGE) &&
            (foreseen == ENT_EVENT_NONE)) {
                event_ = ENT_EVENT_INTERACTION;
                /* Randomise the interaction process and the target. */
                const double r = cs[PROGET_N - 1] * context->random(context);
                if (r < 0.)
                        rc = ENT_RETURN_DOMAIN_ERROR;
                else {
                        int proget;
                        for (proget = 0; proget < PROGET_N; proget++)
                                if (r <= cs[proget]) break;
                        rc = transport_vertex(
                            physics, context, proget, neutrino, product);
                }
        }

exit:
        if (event != NULL) *event = event_;
        ENT_RETURN(rc);
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
