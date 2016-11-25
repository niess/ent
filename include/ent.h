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

#ifndef ENT_H
#define ENT_H
#ifdef __cplusplus
extern "C" {
#endif

/* For C standard streams. */
#ifndef FILE
#include <stdio.h>
#endif

/**
* Return codes used by ENT.
*/
enum ent_return {
        /** The operation succeeded. */
        ENT_RETURN_SUCCESS = 0,
        /** A wrong pointer address was provided, e.g. NULL. */
        ENT_RETURN_BAD_ADDRESS,
        /** Some input parameters are out of their validity range. */
        ENT_RETURN_DOMAIN_ERROR,
        /** Some input file has a wrong format. */
        ENT_RETURN_FORMAT_ERROR,
        /** Some read /write error occured. */
        ENT_RETURN_IO_ERROR,
        /** Some memory couldn't be allocated. */
        ENT_RETURN_MEMORY_ERROR,
        /** Some file couldn't be found. */
        ENT_RETURN_PATH_ERROR,
        /** The number of return codes. */
        ENT_N_RETURNS
};

/**
 *  Projectiles, i.e. neutrinos flavours.
 */
enum ent_projectile {
        /** The tau anti-neutrino. */
        ENT_PROJECTILE_NU_TAU_BAR = -3,
        /** The muon anti-neutrino. */
        ENT_PROJECTILE_NU_MU_BAR = -2,
        /** The electron anti-neutrino. */
        ENT_PROJECTILE_NU_E_BAR = -1,
        /** The electron neutrino. */
        ENT_PROJECTILE_NU_E = 1,
        /** The muon neutrino. */
        ENT_PROJECTILE_NU_MU = 2,
        /** The tau neutrino. */
        ENT_PROJECTILE_NU_TAU = 3,
        /* The number of projectiles. */
        ENT_N_PROJECTILES = 6
};

/**
 *  Neutrino interaction processes.
 */
enum ent_process {
        /** The elastic scattering on electrons, e.g. nu+e->nu+e. */
        ENT_PROCESS_ELASTIC = 0,
        /** The neutral current DIS process. */
        ENT_PROCESS_DIS_NC,
        /** The charged current DIS process. */
        ENT_PROCESS_DIS_CC,
        /** The inverse muon decay : nu_mu+e->nu_e+mu. */
        ENT_PROCESS_INVERSE_MUON,
        /** The inverse tau decay : nu_tau+e->nu_tau+mu. */
        ENT_PROCESS_INVERSE_TAU,
        /* The number of processes. */
        ENT_N_PROCESSES
};

/**
* Generic function pointer.
*
* This is a generic function pointer used to identify the library functions,
* e.g. for error handling.
*/
typedef void ent_function_t(void);

/**
 * Callback for error handling.
 *
 * @param code     The `ent_return` error code.
 * @param caller   The library calling function.
 *
 * The user might provide its own error handler. It will be called at the
 * return of any ENT library function providing an error code.
 */
typedef void ent_handler_cb(enum ent_return code, ent_function_t caller);

/**
 * Opaque structure for handling neutrinos DCS data.
 */
struct ent_dcs;

/**
 * Create a new DCS set.
 *
 * @param data    The data file to create the set from.
 * @param dcs     A handle for the DCS set.
 *
 * Create a new DCS set from tabulated parton distributions data.
 */
enum ent_return ent_dcs_create(const char * data, struct ent_dcs ** dcs);

/**
 * Destroy a DCS set.
 *
 * @param dcs     A handle for the DCS set.
 *
 * Release the allocated memory. On exit _dcs_ is set to `NULL`.
 */
void ent_dcs_destroy(struct ent_dcs ** dcs);

/**
 * Compute the DCS.
 *
 * @param dcs           A handle for the DCS set.
 * @param projectile    The incoming projectile.
 * @param energy        The projectile total energy.
 * @param Z             The target charge number.
 * @param A             The target atomic mass number.
 * @param process       The interaction process.
 * @param x             The Bjorken _x_ fraction.
 * @param y             The energy loss fraction.
 * @param value         The coresponding DCS.
 *
 * Compute the DCS for the given projectile incoming on a target (Z, A) at
 * rest. For an isoscalar nucleon set `Z = 0.5` and `A = 1`.
 */
enum ent_return ent_dcs_compute(struct ent_dcs * dcs,
    enum ent_projectile projectile, double energy, double Z, double A,
    enum ent_process process, double x, double y, double * value);

/**
 * Return a string describing a `ent_return` code.
 *
 * @param rc    The return code.
 * @return A static string.
 *
 * This function is analog to the standard C `strerror` function but specific
 * to ENT return codes. It is thread safe.
 */
const char * ent_error_string(enum ent_return code);

/**
 * Return a string describing a ENT library function.
 *
 * @param function    The library function.
 * @return a static string.
 *
 * This function is meant for verbosing when handling errors. It is thread
 * safe.
 */
const char * ent_error_function(ent_function_t * function);

/**
 * Print a formated summary of an error.
 *
 * @param stream        The output stream where to print.
 * @param code          The `ent_return` error code.
 * @param caller        The library calling function.
 * @param tabulation    The tabulation separator or `NULL`.
 * @param newline       The newline separator or `NULL`.
 *
 * The output summary is formated in JSON. The *tabulation* and *newline*
 * parameters allow to control the output's rendering.
 */
void ent_error_print(FILE * stream, enum ent_return code, ent_function_t caller,
    const char * tabulation, const char * newline);

/**
 * Set or clear the error handler.
 *
 * @param handler    The error handler to set or `NULL`.
 *
 * Set the error handler callback for ENT library functions. If *handler* is
 * set to `NULL` error callbacks are disabled.
 *
 * __Warnings__
 *
 * The error handler is thread local.
 */
void ent_error_handler_set(ent_handler_cb * handler);

/**
 * Get the current error handler.
 *
 * @return The current error handler or `NULL` if none.
 *
 * __Warnings__
 *
 * The error handler is thread local.
 */
ent_handler_cb * ent_error_handler_get();

#ifdef __cplusplus
}
#endif
#endif
