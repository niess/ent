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
        /** The elastic scattering on electrons, e.g. nu + e -> nu + e. */
        ENT_PROCESS_ELASTIC = 0,
        /** The neutral current DIS process. */
        ENT_PROCESS_DIS_NC,
        /** The charged current DIS process. */
        ENT_PROCESS_DIS_CC,
        /** The inverse muon decay, e.g. nu_mu + e- -> nu_e + mu-. */
        ENT_PROCESS_INVERSE_MUON,
        /** The inverse tau decay, e.g. nu_tau + e- -> nu_e + tau-. */
        ENT_PROCESS_INVERSE_TAU,
        /** The resonant hadron(s) production in nu_e~ + e- -> W -> h. */
        ENT_PROCESS_GLASHOW_HADRON,
        /* The number of processes. */
        ENT_N_PROCESSES
};

/**
 *  Parton codes.
 */
enum ent_parton {
        /** The top anti-quark. */
        ENT_PARTON_T_BAR = -6,
        /** The beauty anti-quark. */
        ENT_PARTON_B_BAR = -5,
        /** The charm anti-quark. */
        ENT_PARTON_C_BAR = -4,
        /** The strange anti-quark. */
        ENT_PARTON_S_BAR = -6,
        /** The up anti-quark. */
        ENT_PARTON_U_BAR = -5,
        /** The down anti-quark */
        ENT_PARTON_D_BAR = -4,
        /** The gluon. */
        ENT_PARTON_GLUON = 0,
        /** The down quark. */
        ENT_PARTON_D = 1,
        /** The up quark */
        ENT_PARTON_U = 2,
        /** The strange quark. */
        ENT_PARTON_S = 3,
        /** The charm quark. */
        ENT_PARTON_C = 4,
        /** The beauty quark. */
        ENT_PARTON_B = 5,
        /** The top quark. */
        ENT_PARTON_T = 6,
        /* The number of partons. */
        ENT_N_PARTONS = 13
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
struct ent_physics;

/**
 * Create a new Physics environment.
 *
 * @param physics    A handle for the Physics environment.
 * @param pdf        The PDF file(s) to create the set from.
 * @return On success `ENT_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Create a new Physics environment conforming to the Standard Model of
 * Particle Physics and according to the given tabulations of Parton
 * Distribution Functions (PDF).
 *
 * __Error codes__
 *
 *     ENT_RETURN_FORMAT_ERROR    The pdf file format is not valid / supported.
 *
 *     ENT_RETURN_MEMORY_ERROR    Couldn't allocate memory.
 *
 *     ENT_RETURN_PATH_ERROR      The pdf file couldn't be found/opened.
 */
enum ent_return ent_physics_create(
    struct ent_physics ** physics, const char * pdf);

/**
 * Destroy a Physics environment.
 *
 * @param physics    A handle for the Physics.
 *
 * Release the allocated memory. On exit _physics_ is set to `NULL`.
 */
void ent_physics_destroy(struct ent_physics ** physics);

/**
 * Compute a DCS.
 *
 * @param physics       A handle for the Physics.
 * @param projectile    The incoming projectile.
 * @param energy        The projectile total energy.
 * @param Z             The target charge number.
 * @param A             The target atomic mass number.
 * @param process       The interaction process.
 * @param x             The Bjorken _x_ fraction.
 * @param y             The energy loss fraction.
 * @param dcs           The coresponding DCS value.
 * @return On success `ENT_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Compute the Differential Cross-Section (DCS) in the Laboratory frame for
 * the given projectile incoming on a target (Z, A) at rest. For an isoscalar
 * nucleon set `Z = 0.5` and `A = 1`.
 *
 * __Error codes__
 *
 *     ENT_RETURN_DOMAIN_ERROR     Some input parameter is invalid.
 */
enum ent_return ent_physics_dcs(struct ent_physics * physics,
    enum ent_projectile projectile, double energy, double Z, double A,
    enum ent_process process, double x, double y, double * dcs);

/**
* Compute a PDF.
*
* @param physics    A handle for the Physics.
* @param parton     The parton of interest.
* @param x          The Bjorken _x_ fraction.
* @param q2         The squared momentum transfer.
* @param pdf        The coresponding PDF value.
* @return On success `ENT_RETURN_SUCCESS` is returned otherwise an error
* code is returned as detailed below.
*
* Compute the Parton Distribution Function (PDF) for the given _parton_ and
* kinematic parameters.
*
* __Error codes__
*
*     ENT_RETURN_DOMAIN_ERROR     Some input parameter is invalid.
*/
enum ent_return ent_physics_pdf(struct ent_physics * physics,
    enum ent_parton parton, double x, double q2, double * value);

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
 * @param function      The library calling function.
 * @param tabulation    The tabulation separator or `NULL`.
 * @param newline       The newline separator or `NULL`.
 *
 * The output summary is formated in JSON. The *tabulation* and *newline*
 * parameters allow to control the output's rendering.
 */
void ent_error_print(FILE * stream, enum ent_return code,
    ent_function_t function, const char * tabulation, const char * newline);

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
