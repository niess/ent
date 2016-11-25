/*
 *  An engine for Neutrinos Transport (ENT)
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
/* The Fermi coupling constant GF/(hbar*c)^3, in GeV^-2. */
#define ENT_PHYS_GF 1.1663787E-05
/* The Planck constant as hbar*c, in GeV * m. */
#define ENT_PHYS_HBC 1.97326978E-16
/* The Weinberg angle at MZ, as sin(theta_W)^2. */
#define ENT_PHYS_SIN_THETA_W_2 0.231295

#ifndef M_PI
/* Define pi, if unknown. */
#define M_PI 3.14159265358979323846
#endif

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
                            realloc(*buffer, (*buffer)->size + FILE_BLOCK_SIZE);
                        if (tmp == NULL) return ENT_RETURN_MEMORY_ERROR;
                        *buffer = tmp;
                        (*buffer)->size += FILE_BLOCK_SIZE;
                        size = FILE_BLOCK_SIZE + 1;
                        ptr += size - 1;
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

/* Locate the next word in the file buffer.
 *
 * __Warning__ : the word is stripped from trailling spaces. That for the
 * input buffer is modified.
 */
static void file_get_word(struct file_buffer * buffer, int skip)
{
        /* Scan the current buffer for a valid word. */
        char * p;
        for (p = buffer->cursor; *p != 0x0; p++) {
                if (*p == ' ') {
                        p++;
                        for (; *p == ' '; p++)
                                ;
                        if (*p == 0x0)
                                file_raise_error(
                                    buffer, ENT_RETURN_FORMAT_ERROR);
                        if (skip-- <= 0) {
                                char * q;
                                for (q = p; (*q != ' ') && (*q != 0x0); q++)
                                        ;
                                *q = 0x0;
                                buffer->cursor = p;
                                return;
                        }
                }
        }

        /* No word found. Raise an error. */
        file_raise_error(buffer, ENT_RETURN_FORMAT_ERROR);
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

static float file_get_float(struct file_buffer * buffer)
{
        enum ent_return rc;
        float f;
        if ((rc = file_get_float_(buffer, &f)) != ENT_RETURN_SUCCESS)
                file_raise_error(buffer, rc);
        return f;
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

/* Convert Lambda_{QCD} to alpha_S. */
static double qcd_lambda_to_alpha(
    unsigned int order, unsigned int nf, double q, double lambda)
{
        const double b0 = 11. - 2. / 3. * nf;
        const double t = log(q / lambda);
        double alpha = 1. / (b0 * t);

        if (order <= 1) return alpha; /* LO. */

        /* NLO. */
        const double b1 = 51. - 19. / 3. * nf;
        return alpha * (1. - b1 / b0 * alpha * log(2. * t));
}

/* Convert alpha_S to Lambda_{QCD}. */
static double qcd_alpha_to_lambda(
    unsigned int order, unsigned int nf, double alpha, double q)
{
        double b0 = 11. - 2. / 3. * nf;
        double t = 1. / (b0 * alpha);

        /* LO. */
        if (order <= 1) return q * exp(-t);

        /* NLO. */
        const double br = (51. - 19. / 3. * nf) / (b0 * b0);
        double as0, as1, ot, lt;

        do {
                lt = log(2. * t) / t;
                ot = t;
                as0 = (1. - br * lt) / (b0 * t);
                as1 = (-1. - br * (1. / t - 2. * lt)) / (b0 * t * t);
                t += (alpha - as0) / as1;
        } while (fabs(ot - t) / ot > 1E-05);
        return q * exp(-t);
}

/* Container for CTEQ PDFs in .tbl or .pds format.
 *
 * The `cteq_` data types and functions have been adapted from file cteqpdf.c
 * from Zoltan Nagy's CTEQ PDF library. Check the following for the original
 * code :
 *      http://www.desy.de/~znagy/Site/CTEQ_PDF.html.
 */
struct cteq_pdf {
        /* The name of the PDF table. */
        char * name;

        /* Alpha_s values. */
        unsigned int order, nf;
        /* Lambda(nf, MSbar) */
        double lambda[7];
        /* Quark masses. */
        double mass[7];
        /* Number of active flavours. */
        unsigned int nfmx;
        /* Maximum number of values? */
        unsigned int mxval;
        /*  Qini, Qmax, xmin  */
        double qini, qmax, xmin;
        /* xv, tv arrays and the pdf grid */
        unsigned int nx, nt;
        float *xv, *xvpow, *tv, *upd;

        /* Placeholder for variable size data. */
        char data[];
};

/* Exponent for CTEQ PDFs. */
#define CTEQ_XPOW 0.3

/* Load PDFs from a .tbl file. */
static enum ent_return tbl_load(FILE * stream, struct cteq_pdf ** pdf)
{
        *pdf = NULL;
        enum ent_return rc;
        struct cteq_pdf header;
        struct file_buffer * buffer = NULL;
        memset(&header, 0x0, sizeof(header));

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

        /* Parse the table name length. */
        int name_length = 0;
        file_get_line(&buffer, 0);
        file_get_word(buffer, 4);
        name_length = strlen(buffer->cursor);

        /* Parse the table format. */
        file_get_line(&buffer, 1);
        header.order = file_get_float(buffer);
        header.nf = file_get_float(buffer);
        header.mass[0] = 0.; /* Gluon mass. */
        header.lambda[header.nf] = file_get_float(buffer);
        int i;
        for (i = 1; i <= 6; i++) /* Quark masses. */
                header.mass[i] = file_get_float(buffer);
        file_get_line(&buffer, 1);
        header.nx = file_get_float(buffer);
        header.nt = file_get_float(buffer);
        header.nfmx = file_get_float(buffer);
        header.mxval = 2;

        /* Allocate and map the memory for the tables. */
        const unsigned int n_upd = (header.nx + 1) * (header.nt + 1) *
            (header.nfmx + 1 + header.mxval);
        const unsigned int size_upd =
            memory_padded_size(n_upd * sizeof(float), sizeof(double));
        const unsigned int size_x =
            memory_padded_size((header.nx + 1) * sizeof(float), sizeof(double));
        const unsigned int size_t =
            memory_padded_size((header.nt + 1) * sizeof(float), sizeof(double));
        const unsigned int size_name =
            memory_padded_size(name_length + 1, sizeof(double));
        const unsigned int size =
            sizeof(header) + size_name + 2 * size_x + size_t + size_upd;
        if ((*pdf = malloc(size)) == NULL) {
                rc = ENT_RETURN_MEMORY_ERROR;
                goto exit;
        }
        memcpy(*pdf, &header, sizeof(header));
        (*pdf)->name = (*pdf)->data;
        (*pdf)->xv = (float *)((char *)((*pdf)->name) + size_name);
        (*pdf)->xvpow = (float *)((char *)((*pdf)->xv) + size_x);
        (*pdf)->tv = (float *)((char *)((*pdf)->xvpow) + size_x);
        (*pdf)->upd = (float *)((char *)((*pdf)->tv) + size_t);

        /* Read Qinit and Qmax. */
        file_get_line(&buffer, 1);
        (*pdf)->qini = file_get_float(buffer);
        (*pdf)->qmax = file_get_float(buffer);

        /* Read the q values and convert them to t. */
        file_get_line(&buffer, 0);
        file_get_table(&buffer, (*pdf)->nt + 1, (*pdf)->tv);
        for (i = 0; i < (*pdf)->nt + 1; i++)
                (*pdf)->tv[i] =
                    log(log((*pdf)->tv[i] / (*pdf)->lambda[(*pdf)->nf]));

        /* Read xmin. */
        file_get_line(&buffer, 1);
        (*pdf)->xmin = file_get_float(buffer);
        file_get_line(&buffer, 0);
        file_get_float(buffer); /* Drop the 1st value. */

        /* Read the x table. */
        (*pdf)->xv[0] = 0.;
        file_get_table(&buffer, (*pdf)->nx, (*pdf)->xv + 1);

        /* Compute the xvpow values. */
        (*pdf)->xvpow[0] = 0.;
        for (i = 1; i <= (*pdf)->nx; i++)
                (*pdf)->xvpow[i] = pow((*pdf)->xv[i], CTEQ_XPOW);

        /* Read the grid. */
        file_get_line(&buffer, 1);
        file_get_table(&buffer, n_upd, (*pdf)->upd);

        /* Precompute the lambda values at the thresholds. */
        double as;
        unsigned int nf;

        for (nf = (*pdf)->nf + 1; nf <= 6; nf++) {
                as = qcd_lambda_to_alpha((*pdf)->order, nf - 1,
                    (*pdf)->mass[nf], (*pdf)->lambda[nf - 1]);
                (*pdf)->lambda[nf] = qcd_alpha_to_lambda(
                    (*pdf)->order, nf, as, (*pdf)->mass[nf]);
        }

        /* Under the charm mass every quark is considered as massless. */
        for (nf = (*pdf)->nf - 1; nf > 2; nf--) {
                as = qcd_lambda_to_alpha((*pdf)->order, nf + 1,
                    (*pdf)->mass[nf + 1], (*pdf)->lambda[nf + 1]);
                (*pdf)->lambda[nf] = qcd_alpha_to_lambda(
                    (*pdf)->order, nf, as, (*pdf)->mass[nf + 1]);
        }

        /* Roll back and copy the table name. */
        rewind(stream);
        file_get_line(&buffer, 0);
        file_get_word(buffer, 4);
        memcpy((*pdf)->name, buffer->cursor, name_length + 1);

exit:
        free(buffer);
        if (rc != ENT_RETURN_SUCCESS) {
                free(*pdf);
                *pdf = NULL;
        }
        return rc;
}

/* TODO: poorly optimised. */
static double cteq_polint4f(float * xa, float * ya, double x)
{
        const double h1 = xa[0] - x;
        const double h2 = xa[1] - x;
        const double h3 = xa[2] - x;
        const double h4 = xa[3] - x;

        const double tmp1 = (ya[1] - ya[0]) / (h1 - h2);
        const double d1 = h2 * tmp1;
        const double c1 = h1 * tmp1;

        const double tmp2 = (ya[2] - ya[1]) / (h2 - h3);
        const double d2 = h3 * tmp2;
        const double c2 = h2 * tmp2;

        const double tmp3 = (ya[3] - ya[2]) / (h3 - h4);
        const double d3 = h4 * tmp3;
        const double c3 = h3 * tmp3;

        const double tmp4 = (c2 - d1) / (h1 - h3);
        const double cd1 = h3 * tmp4;
        const double cc1 = h1 * tmp4;

        const double tmp5 = (c3 - d2) / (h2 - h4);
        const double cd2 = h4 * tmp5;
        const double cc2 = h2 * tmp5;

        const double tmp6 = (cc2 - cd1) / (h1 - h4);
        const double dd1 = h4 * tmp6;
        const double dc1 = h1 * tmp6;

        if (h3 + h4 < 0.0)
                return ya[3] + d3 + cd2 + dd1;
        else if (h2 + h3 < 0.0)
                return ya[2] + d2 + cd1 + dc1;
        else if (h1 + h2 < 0.0)
                return ya[1] + c2 + cd1 + dc1;
        else
                return ya[0] + c1 + cc1 + dc1;
}

static double cteq_pdf_compute(
    const struct cteq_pdf * pdf, int pid, double x, double q)
{
        const double onep = 1.00001;
        const unsigned int nx = pdf->nx, nq = pdf->nt;

        /* Check the inputs and locate the table indices by dichotomy. */
        if (q <= pdf->lambda[pdf->nf]) return 0.;
        if ((pid != 0) && (q <= pdf->mass[abs(pid)])) return 0.;

        int jlx = -1;
        int ju = nx + 1;
        while (ju - jlx > 1) {
                const int jm = (ju + jlx) / 2;
                if (x >= pdf->xv[jm])
                        jlx = jm;
                else
                        ju = jm;
        }

        int jx;
        if (jlx <= -1)
                return 0.; /* x < 0. */
        else if (jlx == 0)
                jx = 0;
        else if (jlx <= (int)nx - 2)
                jx = jlx - 1;
        else if (jlx == (int)nx - 1 || x < onep)
                jx = jlx - 2;
        else
                return 0.; /* x > 1. */

        const double tt = log(log(q / pdf->lambda[pdf->nf]));
        int jlq = -1;
        ju = nq + 1;
        while (ju - jlq > 1) {
                const int jm = (ju + jlq) / 2;
                if (tt >= pdf->tv[jm])
                        jlq = jm;
                else
                        ju = jm;
        }

        int jq;
        if (jlq <= 0)
                jq = 0;
        else if (jlq <= (int)nq - 2)
                jq = jlq - 1;
        else
                jq = nq - 3;

        /* Get the pdf function values at the lattice points. */
        const int ip = (pid > (int)pdf->mxval ? -pid : pid);
        const int jtmp =
            ((ip + pdf->nfmx) * (nq + 1) + (jq - 1)) * (nx + 1) + jx + 1;
        const double ss = pow(x, CTEQ_XPOW);

        float fvec[4];
        if (jx == 0) {
                float fij[4];
                int it;
                for (it = 0; it < 4; ++it) {
                        const int j1 = jtmp + (it + 1) * (nx + 1);
                        fij[0] = 0.;
                        fij[1] = (pdf->upd[j1]) * (pdf->xv[1]) * (pdf->xv[1]);
                        fij[2] =
                            (pdf->upd[j1 + 1]) * (pdf->xv[2]) * (pdf->xv[2]);
                        fij[3] =
                            (pdf->upd[j1 + 2]) * (pdf->xv[3]) * (pdf->xv[3]);

                        fvec[it] = cteq_polint4f(pdf->xvpow, fij, ss) / (x * x);
                }
        } else if (jlx == nx - 1) {
                int it;
                for (it = 0; it < 4; ++it)
                        fvec[it] = cteq_polint4f(pdf->xvpow + nx - 3,
                            pdf->upd + jtmp + (it + 1) * (nx + 1) - 1, ss);
        } else {
                const float * const svec = pdf->xvpow + jx - 1;
                const double s12 = svec[1] - svec[2];
                const double s13 = svec[1] - svec[3];
                const double s23 = svec[2] - svec[3];
                const double s24 = svec[2] - svec[4];
                const double s34 = svec[3] - svec[4];
                const double sy2 = ss - svec[2];
                const double sy3 = ss - svec[3];
                const double const1 = s13 / s23;
                const double const2 = s12 / s23;
                const double const3 = s34 / s23;
                const double const4 = s24 / s23;
                const double s1213 = s12 + s13;
                const double s2434 = s24 + s34;
                const double sdet = s12 * s34 - s1213 * s2434;
                const double tmp = sy2 * sy3 / sdet;
                const double const5 = (s34 * sy2 - s2434 * sy3) * tmp / s12;
                const double const6 = (s1213 * sy2 - s12 * sy3) * tmp / s34;

                int it;
                for (it = 0; it < 4; ++it) {
                        const int j1 = jtmp + (it + 1) * (nx + 1);
                        const double sf2 = pdf->upd[j1];
                        const double sf3 = pdf->upd[j1 + 1];
                        const double g1 = sf2 * const1 - sf3 * const2;
                        const double g4 = sf3 * const4 - sf2 * const3;
                        fvec[it] = (const5 * (pdf->upd[j1 - 1] - g1) +
                                       const6 * (pdf->upd[j1 + 2] - g4) +
                                       sf2 * sy3 - sf3 * sy2) /
                            s23;
                }
        }

        /* Interpolate in t. */
        if (jlq <= 0)
                return cteq_polint4f(pdf->tv, fvec, tt);
        else if (jlq >= nq - 1)
                return cteq_polint4f(pdf->tv + nq - 3, fvec, tt);

        const float * const tvec = pdf->tv + jq - 1;
        const double t12 = tvec[1] - tvec[2];
        const double t13 = tvec[1] - tvec[3];
        const double t23 = tvec[2] - tvec[3];
        const double t24 = tvec[2] - tvec[4];
        const double t34 = tvec[3] - tvec[4];
        const double ty2 = tt - tvec[2];
        const double ty3 = tt - tvec[3];
        const double tmp1 = t12 + t13;
        const double tmp2 = t24 + t34;
        const double tdet = t12 * t34 - tmp1 * tmp2;
        const double g1 = (fvec[1] * t13 - fvec[2] * t12) / t23;
        const double g4 = (fvec[2] * t24 - fvec[1] * t34) / t23;
        const double h00 = (t34 * ty2 - tmp2 * ty3) * (fvec[0] - g1) / t12 +
            (tmp1 * ty2 - t12 * ty3) * (fvec[3] - g4) / t34;

        return (h00 * ty2 * ty3 / tdet + fvec[1] * ty3 - fvec[2] * ty2) / t23;
}

#undef CTEQ_XPOW /* No more needed. */

enum ent_return ent_dcs_create(const char * data, struct ent_dcs ** dcs)
{
        enum ent_return rc;
        *dcs = NULL;

        FILE * stream;
        rc = ENT_RETURN_PATH_ERROR;
        if ((stream = fopen(data, "r")) == NULL) goto exit;

        struct cteq_pdf * pdf;
        if ((rc = tbl_load(stream, &pdf)) != ENT_RETURN_SUCCESS) goto exit;
        *dcs = (struct ent_dcs *)pdf;
        rc = ENT_RETURN_SUCCESS;

exit:
        if (stream != NULL) fclose(stream);
        return rc;
}

void ent_dcs_destroy(struct ent_dcs ** dcs)
{
        if ((dcs == NULL) || (*dcs == NULL)) return;

        free(*dcs);
        *dcs = NULL;
}

/* DCS for Deep Inelastic Scattering (DIS). */
static double dcs_dis(struct ent_dcs * dcs, enum ent_projectile projectile,
    double energy, double Z, double A, enum ent_process process, double x,
    double y)
{
        /* Compute the PDF. */
        const double Q2 = 2. * x * y * ENT_MASS_NUCLEON * energy;
        const double q = sqrt(Q2);
        struct cteq_pdf * pdf = (struct cteq_pdf *)dcs;
        int eps = (projectile > 0) ? 1 : -1; /* CP? */
        const double N = A - Z;              /* Number of neutrons. */
        double factor;
        if (process == ENT_PROCESS_DIS_CC) {
                /* Charged current DIS process. */
                const double d =
                    (Z <= 0.) ? 0. : cteq_pdf_compute(pdf, 1 * eps, x, q);
                const double u =
                    (N <= 0.) ? 0. : cteq_pdf_compute(pdf, 2 * eps, x, q);
                const double s = cteq_pdf_compute(pdf, 3 * eps, x, q);
                const double b = cteq_pdf_compute(pdf, 5 * eps, x, q);
                const double F1 = Z * d + N * u + A * (s + b);

                const double dbar =
                    (N <= 0.) ? 0. : cteq_pdf_compute(pdf, -1 * eps, x, q);
                const double ubar =
                    (Z <= 0.) ? 0. : cteq_pdf_compute(pdf, -2 * eps, x, q);
                const double cbar = cteq_pdf_compute(pdf, -4 * eps, x, q);
                const double F2 = N * dbar + Z * ubar + A * cbar;

                const double y1 = 1. - y;
                const double F = F1 + F2 * y1 * y1;
                const double MW2 = ENT_MASS_W * ENT_MASS_W;
                const double r = MW2 / (MW2 + Q2);
                factor = 2 * F * r * r;
        } else {
                /* Neutral current DIS process. */
                const double d = cteq_pdf_compute(pdf, 1 * eps, x, q);
                const double u = cteq_pdf_compute(pdf, 2 * eps, x, q);
                const double s = cteq_pdf_compute(pdf, 3 * eps, x, q);
                const double c = cteq_pdf_compute(pdf, 4 * eps, x, q);
                const double b = cteq_pdf_compute(pdf, 5 * eps, x, q);

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

                const double dbar = cteq_pdf_compute(pdf, -1 * eps, x, q);
                const double ubar = cteq_pdf_compute(pdf, -2 * eps, x, q);
                const double sbar = cteq_pdf_compute(pdf, -3 * eps, x, q);
                const double cbar = cteq_pdf_compute(pdf, -4 * eps, x, q);
                const double bbar = cteq_pdf_compute(pdf, -5 * eps, x, q);
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
        return energy * x * factor * (ENT_PHYS_GF * ENT_PHYS_HBC) *
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

/* Generic API function for computing the DCS. */
enum ent_return ent_dcs_compute(struct ent_dcs * dcs,
    enum ent_projectile projectile, double energy, double Z, double A,
    enum ent_process process, double x, double y, double * value)
{
        /* Check the inputs. */
        *value = 0.;
        if ((x > 1.) || (x < 0.) || (y > 1.) || (y < 0.))
                return ENT_RETURN_DOMAIN_ERROR;

        /* Compute the corresponding DCS. */
        if (process == ENT_PROCESS_ELASTIC)
                *value = dcs_elastic(projectile, energy, Z, y);
        else if ((process == ENT_PROCESS_DIS_CC) ||
            (process == ENT_PROCESS_DIS_NC))
                *value = dcs_dis(dcs, projectile, energy, Z, A, process, x, y);
        else if ((process == ENT_PROCESS_INVERSE_MUON) ||
            (process == ENT_PROCESS_INVERSE_TAU))
                *value = dcs_inverse(projectile, energy, process, Z, y);
        else
                return ENT_RETURN_DOMAIN_ERROR;

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
