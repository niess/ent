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

#ifndef ANT_H
#define ANT_H
#ifdef __cplusplus
extern "C" {
#endif

/* For C standard streams. */
#ifndef FILE
#include <stdio.h>
#endif

/**
* Return codes used by ANT.
*/
enum ant_return {
        /** The operation succeeded. */
        ANT_RETURN_SUCCESS = 0,
        /** A wrong pointer address was provided, e.g. NULL. */
        ANT_RETURN_BAD_ADDRESS,
        /** Some input parameters are out of their validity range. */
        ANT_RETURN_DOMAIN_ERROR,
        /** Some input file has a wrong format. */
        ANT_RETURN_FORMAT_ERROR,
        /** Some read /write error occured. */
        ANT_RETURN_IO_ERROR,
        /** Some memory couldn't be allocated. */
        ANT_RETURN_MEMORY_ERROR,
        /** Some file couldn't be found. */
        ANT_RETURN_PATH_ERROR,
        /** Reference to an unknown model. */
        ANT_N_RETURNS
};

/**
* Generic function pointer.
*
* This is a generic function pointer used to identify the library functions,
* e.g. for error handling.
*/
typedef void ant_function_t(void);

/**
 * Callback for error handling.
 *
 * @param code     The `ant_return` error code.
 * @param caller   The library calling function.
 *
 * The user might provide its own error handler. It will be called at the
 * return of any ANT library function providing an error code.
 */
typedef void ant_handler_cb(enum ant_return code, ant_function_t caller);

/**
 * Opaque structure for handling neutrinos DCS data.
 */
struct ant_dcs;

/**
 * Create a new DCS set.
 *
 * @param data    The data file to create the set from.
 * @param dcs     A handle for the DCS set.
 *
 * Create a new DCS set from tabulated parton distributions data.
 */
enum ant_return ant_dcs_create(const char * data, struct ant_dcs ** dcs);

/**
 * Destroy a DCS set.
 *
 * @param dcs     A handle for the DCS set.
 *
 * Release the allocated memory. On exit _dcs_ is set to `NULL`.
 */
void ant_dcs_destroy(struct ant_dcs ** dcs);

/**
 * Return a string describing a `ant_return` code.
 *
 * @param rc    The return code.
 * @return A static string.
 *
 * This function is analog to the standard C `strerror` function but specific
 * to ANT return codes. It is thread safe.
 */
const char * ant_error_string(enum ant_return code);

/**
 * Return a string describing a ANT library function.
 *
 * @param function    The library function.
 * @return a static string.
 *
 * This function is meant for verbosing when handling errors. It is thread
 * safe.
 */
const char * ant_error_function(ant_function_t * function);

/**
 * Print a formated summary of an error.
 *
 * @param stream        The output stream where to print.
 * @param code          The `ant_return` error code.
 * @param caller        The library calling function.
 * @param tabulation    The tabulation separator or `NULL`.
 * @param newline       The newline separator or `NULL`.
 *
 * The output summary is formated in JSON. The *tabulation* and *newline*
 * parameters allow to control the output's rendering.
 */
void ant_error_print(FILE * stream, enum ant_return code, ant_function_t caller,
        const char * tabulation, const char * newline);

/**
 * Set or clear the error handler.
 *
 * @param handler    The error handler to set or `NULL`.
 *
 * Set the error handler callback for ANT library functions. If *handler* is
 * set to `NULL` error callbacks are disabled.
 *
 * __Warnings__
 *
 * The error handler is thread local.
 */
void ant_error_handler_set(ant_handler_cb * handler);

/**
 * Get the current error handler.
 *
 * @return The current error handler or `NULL` if none.
 *
 * __Warnings__
 *
 * The error handler is thread local.
 */
ant_handler_cb * ant_error_handler_get();

#ifdef __cplusplus
}
#endif
#endif
