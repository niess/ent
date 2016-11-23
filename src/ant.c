/*
 *  An engine for Neutrinos Transport (ANT)
 *  Copyright (C) 2016  Valentin Niess
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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

/* The ANT API. */
#include "ant.h"

/* Container for reading data from a text file. */
struct file_buffer {
        jmp_buf env;
        FILE * stream;
        enum ant_return code;
        unsigned int size;
        char * cursor;
        char line[];
};

/* Memory blocks size for file buffers. */
#define FILE_BLOCK_SIZE 2048

/* Create a new file buffer */
static enum ant_return file_buffer_create(
        struct file_buffer ** buffer, FILE * stream, jmp_buf env)
{
        if ((*buffer = malloc(FILE_BLOCK_SIZE)) == NULL)
                return ANT_RETURN_MEMORY_ERROR;
        memcpy((*buffer)->env, env, sizeof(jmp_buf));
        (*buffer)->stream = stream;
        (*buffer)->code = ANT_RETURN_SUCCESS;
        (*buffer)->size = FILE_BLOCK_SIZE-sizeof(**buffer);
        (*buffer)->cursor = NULL;

        return ANT_RETURN_SUCCESS;
}

/* Raise a file parsing error. */
static void file_raise_error(struct file_buffer * buffer, enum ant_return code)
{
        buffer->code = code;
        longjmp(buffer->env, 1);
}

/* Utility functions for reading a line from a file. */
static enum ant_return file_get_line_(struct file_buffer ** buffer, int skip)
{
        /* Get a new line. */
        char * ptr = (*buffer)->line;
        int size = (*buffer)->size;
        for (; skip >= 0; skip--) {
                for (;;) {
                        const char check = 0x1;
                        ptr[size-1] = check;
                        if (fgets(ptr, size, (*buffer)->stream) == NULL)
                                return ANT_RETURN_IO_ERROR;

                        /* Check that no overflow occured. */
                        if (ptr[size-1] == check) break;

                        /* Get more memory then ... */
                        void * tmp = realloc(*buffer, (*buffer)->size +
                            FILE_BLOCK_SIZE);
                        if (tmp == NULL) return ANT_RETURN_MEMORY_ERROR;
                        *buffer = tmp;
                        (*buffer)->size += FILE_BLOCK_SIZE;
                        size = FILE_BLOCK_SIZE+1;
                        ptr += size-1;
                }
        }

        /* Reset the line cursor and return. */
        (*buffer)->cursor = (*buffer)->line;
        return ANT_RETURN_SUCCESS;
}

#undef FILE_BLOCK_SIZE /* No more needed. */

static void file_get_line(struct file_buffer ** buffer, int skip)
{
        enum ant_return rc;
        if ((rc = file_get_line_(buffer, skip)) != ANT_RETURN_SUCCESS)
                file_raise_error(*buffer, rc);
}

/* Utility function for getting the next word in the buffer.
 *
 * __Warning__ : the word is stripped from trailling spaces. That for the
 * input buffer data are modified.
 */
static void file_get_word(struct file_buffer * buffer, int skip)
{
        /* Scan the current buffer for a valid word. */
        char * p;
        for (p = buffer->cursor; *p != 0x0; p++) {
                if (*p == ' ') {
                        p++;
                        for (; *p == ' '; p++) ;
                        if (*p == 0x0) file_raise_error(buffer,
                            ANT_RETURN_FORMAT_ERROR);
                        if (skip-- <= 0) {
                                char * q;
                                for (q = p; (*q != ' ') && (*q != 0x0); q++) ;
                                *q = 0x0;
                                buffer->cursor = p;
                                return;
                        }
                }
        }

        /* No word found. Raise an error. */
        file_raise_error(buffer, ANT_RETURN_FORMAT_ERROR);
}

/* Utility functions for parsing the next float in the buffer. */
static enum ant_return file_get_float_(struct file_buffer * buffer, float * f)
{
        char * endptr;
        *f = strtof(buffer->cursor, &endptr);
        if (buffer->cursor == endptr) return ANT_RETURN_FORMAT_ERROR;
        buffer->cursor = endptr;
        return ANT_RETURN_SUCCESS;
}

static float file_get_float(struct file_buffer * buffer)
{
        enum ant_return rc;
        float f;
        if ((rc = file_get_float_(buffer, &f)) != ANT_RETURN_SUCCESS)
                file_raise_error(buffer, rc);
        return f;
}

/* Utility function for parsing a table. */
static void file_get_table(
    struct file_buffer ** buffer, int size, float * table)
{
        enum ant_return rc;
        int failed = 0, n = 0;
        while (n < size) {
                if ((rc = file_get_float_(*buffer, table)) !=
                    ANT_RETURN_SUCCESS) {
                        if (failed) file_raise_error(*buffer, rc);
                        if ((rc = file_get_line_(buffer, 0))
                            != ANT_RETURN_SUCCESS)
                                file_raise_error(*buffer, rc);
                        failed = 1;
                }
                else {
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
        return i*pad_size;
}

static double qcd_lambda_to_alpha(
    unsigned int order, unsigned int nf, double q, double lambda)
{
        const double b0 = 11. - 2. / 3. * nf;
        const double t = log(q / lambda);
        double alpha = 1./(b0 * t);

        if(order <= 1) return alpha; /* LO. */

        /* NLO. */
        const double b1 = 51. - 19. / 3. * nf;
        return alpha*(1. - b1 / b0 * alpha * log(2. * t));
}

static double qcd_alpha_to_lambda(
        unsigned int order, unsigned int nf, double alpha, double q)
{
        double b0 = 11. - 2. / 3. * nf;
        double t = 1./(b0 * alpha);

        /* LO. */
        if(order <= 1) return q * exp(-t);

        /* NLO. */
        const double br = (51. - 19. / 3. * nf) / (b0 * b0);
        double as0, as1, ot, lt;

        do {
                lt = log(2. * t) / t;
                ot = t;
                as0 = (1. - br * lt) / (b0 * t);
                as1 = (-1. - br * (1. / t - 2. * lt)) / (b0 * t * t);
                t += (alpha - as0) / as1;
        }
        while(fabs(ot - t) / ot > 1E-05);
        return q * exp(-t);
}

/* Container for PDFs in .tbl format. */
struct tbl_pdf {
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
        /*  Qini, Qmax, xmin  */
        double qini, qmax, xmin;
        /* xv, tv arrays and the pdf grid */
        unsigned int nx, nt;
        float * xv, * xvpow, * tv, * upd;

        /* Placeholder for variable size data. */
        char data[];
};

/* Exponent for CTEQ PDFs. */
#define CTEQ_XPOW 0.3

/* Load PDFs from a .tbl file.
 *
 * Adapted from file cteqpdf.c from Zoltan Nagy's CTEQ PDF library. See :
 *      http://www.desy.de/~znagy/Site/CTEQ_PDF.html.
 */
static enum ant_return tbl_load(FILE * stream, struct tbl_pdf ** pdf)
{
        *pdf = NULL;
        enum ant_return rc;
        struct tbl_pdf header;
        struct file_buffer * buffer = NULL;
        memset(&header, 0x0, sizeof(header));

        /* Redirect any subsequent file parsing error. */
        jmp_buf env;
        if (setjmp(env) != 0) {
                rc = buffer->code;
                goto exit;
        }

        /* Create the temporary file buffer. */
        if ((rc = file_buffer_create(&buffer, stream, env))
            != ANT_RETURN_SUCCESS)
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

        /* Allocate and map the memory for the tables. */
        const unsigned int n_upd = (header.nx+1)*(header.nt+1)*(header.nfmx+3);
        const unsigned int size_upd = memory_padded_size(n_upd * sizeof(float),
            sizeof(double));
        const unsigned int size_x = memory_padded_size((header.nx+1) *
            sizeof(float), sizeof(double));
        const unsigned int size_t = memory_padded_size((header.nt+1) *
            sizeof(float), sizeof(double));
        const unsigned int size_name = memory_padded_size(name_length + 1,
            sizeof(double));
        const unsigned int size = sizeof(header)+size_name+2*size_x+size_t+
            size_upd;
        if ((*pdf = malloc(size)) == NULL) {
                rc = ANT_RETURN_MEMORY_ERROR;
                goto exit;
        }
        memcpy(*pdf, &header, sizeof(header));
        (*pdf)->name = (*pdf)->data;
        (*pdf)->xv = (float *)((char *)((*pdf)->name)+size_name);
        (*pdf)->xvpow = (float *)((char *)((*pdf)->xv)+size_x);
        (*pdf)->tv = (float *)((char *)((*pdf)->xvpow)+size_x);
        (*pdf)->upd = (float *)((char *)((*pdf)->tv)+size_t);

        /* Read Qinit and Qmax. */
        file_get_line(&buffer, 1);
        (*pdf)->qini = file_get_float(buffer);
        (*pdf)->qmax = file_get_float(buffer);

        /* Read the q values and convert them to t. */
        file_get_line(&buffer, 0);
        file_get_table(&buffer, (*pdf)->nt+1, (*pdf)->tv);
        for (i = 0; i < (*pdf)->nt+1; i++)
                (*pdf)->tv[i] = log(log((*pdf)->tv[i] /
                    (*pdf)->lambda[(*pdf)->nf]));

        /* Read xmin. */
        file_get_line(&buffer, 1);
        (*pdf)->xmin = file_get_float(buffer);
        file_get_line(&buffer, 0);
        file_get_float(buffer); /* Drop the 1st value. */

        /* Read the x table. */
        (*pdf)->xv[0] = 0.;
        file_get_table(&buffer, (*pdf)->nx, (*pdf)->xv+1);

        /* Compute the xvpow values. */
        (*pdf)->xvpow[0] = 0.;
        for(i = 1; i <= (*pdf)->nx; i++) (*pdf)->xvpow[i] = pow((*pdf)->xv[i],
            CTEQ_XPOW);

        /* Read the grid. */
        file_get_line(&buffer, 1);
        file_get_table(&buffer, n_upd, (*pdf)->upd);

        /* Precompute the lambda values at the thresholds. */
        double as;
        unsigned int nf;

        for (nf = (*pdf)->nf+1; nf <= 6; nf++) {
                as = qcd_lambda_to_alpha((*pdf)->order, nf-1,
                    (*pdf)->mass[nf], (*pdf)->lambda[nf-1]);
                (*pdf)->lambda[nf] = qcd_alpha_to_lambda((*pdf)->order, nf, as,
                    (*pdf)->mass[nf]);
        }

        /* Under the charm mass every quark is considered as massless. */
        for (nf = (*pdf)->nf-1; nf > 2; nf--) {
                as = qcd_lambda_to_alpha((*pdf)->order, nf+1,
                    (*pdf)->mass[nf+1], (*pdf)->lambda[nf+1]);
                (*pdf)->lambda[nf] = qcd_alpha_to_lambda((*pdf)->order, nf, as,
                    (*pdf)->mass[nf+1]);
        }

        /* Roll back and copy the table name. */
        rewind(stream);
        file_get_line(&buffer, 0);
        file_get_word(buffer, 4);
        memcpy((*pdf)->name, buffer->cursor, name_length+1);

exit:
        free(buffer);
        if (rc != ANT_RETURN_SUCCESS) {
                free(*pdf);
                *pdf = NULL;
        }
        return rc;
}

#undef CTEQ_XPOW /* No more needed. */

enum ant_return ant_dcs_create(const char * data, struct ant_dcs ** dcs)
{
        enum ant_return rc;
        *dcs = NULL;

        FILE * stream;
        rc = ANT_RETURN_PATH_ERROR;
        if ((stream = fopen(data, "r")) == NULL) goto exit;

        struct tbl_pdf * pdf;
        if ((rc = tbl_load(stream, &pdf)) != ANT_RETURN_SUCCESS)
                goto exit;
        *dcs = (struct ant_dcs *)pdf;
        rc = ANT_RETURN_SUCCESS;

exit:
        if (stream != NULL) fclose(stream);
        return rc;
}

void ant_dcs_destroy(struct ant_dcs ** dcs)
{
        if ((dcs == NULL) || (*dcs == NULL)) return;

        free(*dcs);
        *dcs = NULL;
}
