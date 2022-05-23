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
        /** Some memory could not be allocated. */
        ENT_RETURN_MEMORY_ERROR,
        /** Some file could not be found. */
        ENT_RETURN_PATH_ERROR,
        /** The number of return codes. */
        ENT_N_RETURNS
};

/**
 *  Particles ID's, following the PDG numbering scheme.
 */
enum ent_pid {
        /** The tau anti-neutrino. */
        ENT_PID_NU_BAR_TAU = -16,
        /** The anti-tau. */
        ENT_PID_TAU_BAR = -15,
        /** The muon anti-neutrino. */
        ENT_PID_NU_BAR_MU = -14,
        /** The anti-muon. */
        ENT_PID_MUON_BAR = -13,
        /** The electron anti-neutrino. */
        ENT_PID_NU_BAR_E = -12,
        /** The positron. */
        ENT_PID_POSITRON = -11,
        /** The anti-top quark. */
        ENT_PID_TOP_BAR = -6,
        /** The anti-top quark. */
        ENT_PID_BOTTOM_BAR = -5,
        /** The anti-bottom quark. */
        ENT_PID_CHARM_BAR = -4,
        /** The anti-strange quark. */
        ENT_PID_STRANGE_BAR = -3,
        /** The anti-up quark. */
        ENT_PID_UP_BAR = -2,
        /** The anti-down quark. */
        ENT_PID_DOWN_BAR = -1,
        /** Tag for none product. */
        ENT_PID_NONE = 0,
        /* The down quark. */
        ENT_PID_DOWN = 1,
        /* The up quark. */
        ENT_PID_UP = 2,
        /* The strange quark. */
        ENT_PID_STRANGE = 3,
        /* The charm quark. */
        ENT_PID_CHARM = 4,
        /* The bottom quark. */
        ENT_PID_BOTTOM = 5,
        /* The top quark. */
        ENT_PID_TOP = 6,
        /** The electron. */
        ENT_PID_ELECTRON = 11,
        /** The electron neutrino. */
        ENT_PID_NU_E = 12,
        /** The muon. */
        ENT_PID_MUON = 13,
        /** The muon neutrino. */
        ENT_PID_NU_MU = 14,
        /** The tau. */
        ENT_PID_TAU = 15,
        /** The tau neutrino. */
        ENT_PID_NU_TAU = 16,
        /** Unresolved hadron(s). */
        ENT_PID_HADRON = 100,
        /** A neutron. */
        ENT_PID_NEUTRON = 2112,
        /** A proton. */
        ENT_PID_PROTON = 2212,
        /** W- boson. */
        ENT_PID_W_MINUS = -24,
        /** Z boson. */
        ENT_PID_Z = 23,
        /** W+ boson. */
        ENT_PID_W_PLUS = 24
};

/**
 *  Neutrino interaction processes.
 */
enum ent_process {
        /** No specific process. */
        ENT_PROCESS_NONE = -1,
        /** The elastic scattering on electrons, e.g. nu + e -> nu + e. */
        ENT_PROCESS_ELASTIC = 0,
        /** The neutral current DIS process. */
        ENT_PROCESS_DIS_NC,
        /** The charged current DIS process (all contributions). */
        ENT_PROCESS_DIS_CC,
        /** The charged current DIS process (w/o top production). */
        ENT_PROCESS_DIS_CC_OTHER,
        /** The charged current DIS process (top production only). */
        ENT_PROCESS_DIS_CC_TOP,
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
        ENT_PARTON_S_BAR = -3,
        /** The up anti-quark. */
        ENT_PARTON_U_BAR = -2,
        /** The down anti-quark */
        ENT_PARTON_D_BAR = -1,
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

struct ent_medium;
struct ent_state;
/**
 * Density callback for a material medium.
 *
 * @param medium      A related `ent_medium` or `NULL`.
 * @param state       The state for which the _density_ is requested.
 * @param density     The corresponding density value in kg/m^3.
 * @return A proposed step limit distance, in m.
 *
 * The callback must return a proposed Monte-Carlo stepping distance, in m,
 * consistent with the size of the propagation medium density inhomogeneities,
 * e. g. 1 % of &rho; / |&nabla; &rho;|. Note that returning zero or less
 * signs that the corresponding medium is uniform.
 *
 * **Warning** : it is an error to return zero or less for any position of the
 * medium if at least one area is not uniform.
 */
typedef double ent_density_cb(
    struct ent_medium * medium, struct ent_state * state, double * density);

/**
 * Data for describing a propagation medium.
 *
 * These are the data required by _ENT_ for describing a propagation medium.
 * The user might implement his own data structure on top of it.
 */
struct ent_medium {
        /** The effective charge number of the medium's material. */
        double Z;
        /** The effective mass number of the medium's material. */
        double A;
        /** The medium density callback. */
        ent_density_cb * density;
};

struct ent_context;
/**
 * Medium callback for a Monte-Carlo state.
 *
 * @param context   The Monte-Carlo context requiring a medium.
 * @param state     The Monte-Carlo state for which the medium is requested.
 * @param medium    A pointer to the corresponding medium or `NULL` if the state
 * has exit the simulation area.
 * @return The proposed step size or `0` for an infinite medium.
 *
 * The callback must propose a Monte-Carlo stepping distance, in m,
 * consistent with the geometry. Note that returning zero or less
 * signs that the corresponding medium has no boundaries.
 *
 * **Warning** : it is an error to return zero or less for any state if the
 * initial medium is finite, i.e. returns a stricly positive step value.
 */
typedef double ent_medium_cb(struct ent_context * context,
    struct ent_state * state, struct ent_medium ** medium);

/**
 * Callback providing a stream of pseudo random numbers.
 *
 * @param context The Monte-Carlo context requiring a random number.
 * @return A uniform pseudo random number in [0;1].
 *
 * **Note** : this is the only random stream used by ENT. The user must unsure
 * proper behaviour, i.e. that a flat distribution in [0;1] is indeed returned.
 *
 * **Warning** : if multiple contexts are used the user must ensure that this
 * callback is thread safe, e.g. by using independant streams for each context
 * or a locking mechanism in order to share a single random stream.
 */
typedef double ent_random_cb(struct ent_context * context);

/**
 * Callback providing an a priori relative density of ancestors.
 *
 * @param context    The Monte-Carlo context requiring an ancestor.
 * @param ancestor   The PID of the required ancestor.
 * @param daughter   The daughter state.
 * @return The relative density of the ancestor.
 *
 * This callback **must** be provided for backward Monte-Carlo. It is expected
 * to return an _a priori_ density estimate for an ancestor at the daughter's
 * location. Returning `0` or less indicates a null density, which disables
 * the corresponding channel.
 *
 * **Warning** : if multiple contexts are used the user must ensure that this
 * callback is thread safe.
 */
typedef double ent_ancestor_cb(struct ent_context * context,
    enum ent_pid ancestor, struct ent_state * daughter);

/**
 * Callback for custom actions during the Monte-Carlo stepping.
 *
 * @param context    The operating Monte-Carlo context.
 * @param medium     The medium at the current step.
 * @param state      The particle state at the current step.
 * @return On success `ENT_RETURN_SUCCESS` must be returned otherwise any error
 * code can be returned.
 *
 * Providing this callback is optionnal. It is called at the end of major
 * Monte-Carlo steps when invoking `ent_transport`.
 *
 * **Warning** : if multiple contexts are used the user must ensure that this
 * callback is thread safe.
 */
typedef enum ent_return ent_stepping_cb(struct ent_context * context,
    struct ent_medium * medium, struct ent_state * state);

/**
 * Data for a Monte-Carlo context.
 *
 * These are the data required by _ENT_ for describing a Monte-Carlo context.
 * The user might implement his own data structure on top of it.
 */
struct ent_context {
        /** The medium callback. */
        ent_medium_cb * medium;
        /** The random engine callback. */
        ent_random_cb * random;
        /** The ancestor callback, or `NULL` for forward Monte-Carlo. */
        ent_ancestor_cb * ancestor;
        /** A custom stepping action callback, or `NULL` if none. */
        ent_stepping_cb * stepping_action;
        /** A user supplied distance limit for the transport, or `0`. */
        double distance_max;
        /** A user supplied grammage limit for the transport, or `0`. */
        double grammage_max;
};

/**
 * Data for a particle Monte-Carlo state.
 *
 * These are the data set required by _ENT_ for describing a particle
 * Monte-Carlo state. The user might implement his own data structure on top
 * of it.
 */
struct ent_state {
        /** The particle type. */
        enum ent_pid pid;
        /** The particle energy, in GeV. */
        double energy;
        /** The total distance travelled by the particle, in m. */
        double distance;
        /** The total grammage travelled by the particle, in kg/m^2. */
        double grammage;
        /** The Monte-Carlo weight. */
        double weight;
        /** The particle absolute position, in m. */
        double position[3];
        /** The particle momentum's direction. */
        double direction[3];
};

/**
 * Create a new physics environment.
 *
 * @param physics     A handle for the physics environment.
 * @param data        Physics data (PDFs, SFs for DIS or a physics dump).
 * @param cs          A cross-section file, for DIS processes, or `NULL`.
 * @return On success `ENT_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Create a new physics environment according to the given physics data. Those
 * can be Parton Distribution Functions (PDFs), DIS Structure Functions (SFs) or
 * a previous physics dump generated with `ent_physics_dump`. If PDFs are
 * provided, then DIS SFs are computed by ENT using LO expressions.
 *
 * __Note__: PDFs files must be in Les Houches Accord (LHA) format. Structure
 * functions are provided using an ENT specific binary format (.ent).
 *
 * If *cs* is not `NULL`, then it must point to file containing cross-section
 * values for DIS processes. ENT's cross-sections are re-scaled accordingly.
 * Otherwise, ENT's cross-sections for DIS processes are computed from physics
 * data.
 *
 * __Error codes__
 *
 *     ENT_RETURN_FORMAT_ERROR    An input file format is not valid / supported.
 *
 *     ENT_RETURN_MEMORY_ERROR    Could not allocate memory.
 *
 *     ENT_RETURN_PATH_ERROR      An input file could not be found/opened.
 */
enum ent_return ent_physics_create(
    struct ent_physics ** physics, const char * data, const char * cs);

/**
 * Destroy a physics environment.
 *
 * @param physics    A handle for the physics.
 *
 * Release the allocated memory. On exit _physics_ is set to `NULL`.
 */
void ent_physics_destroy(struct ent_physics ** physics);

/**
 * Dump physics data to a file.
 *
 * @param physics    A handle for the physics.
 * @param outfile    File where to dump the physics data.
 *
 * Dump physics data to a file using an ENT specific binary format (.ent). New
 * physics instances can be created directly from these data, using
 * `ent_physics_create`. This allows for faster initialisation.
 *
 * __Error codes__
 *
 *     ENT_RETURN_IO_ERROR        Could not write to output file.
 *
 *     ENT_RETURN_PATH_ERROR      The output file could not be found/opened.
 */
enum ent_return ent_physics_dump(
    const struct ent_physics * physics, const char * outfile);

/**
 * The total or process specific cross-section.
 *
 * @param physics          A handle for the physics.
 * @param projectile       The incoming projectile.
 * @param energy           The projectile total energy.
 * @param Z                The target charge number.
 * @param A                The target atomic mass number.
 * @param process          The interaction process.
 * @param cross_section    The coresponding DCS value.
 * @return On success `ENT_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Compute the cross-secction for the given projectile incoming on a
 * target (Z, A). Set *process* to `ENT_PROCESS_NONE` in order to Compute
 * the total cross-section.
 *
 * __Error codes__
 *
 *     ENT_RETURN_DOMAIN_ERROR     Some input parameter is invalid.
 */
enum ent_return ent_physics_cross_section(struct ent_physics * physics,
    enum ent_pid projectile, double energy, double Z, double A,
    enum ent_process process, double * cross_section);

/**
 * Compute a DCS.
 *
 * @param physics       A handle for the physics.
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
    enum ent_pid projectile, double energy, double Z, double A,
    enum ent_process process, double x, double y, double * dcs);

/**
* Compute a PDF.
*
* @param physics    A handle for the physics.
* @param parton     The parton of interest.
* @param x          The Bjorken _x_ fraction.
* @param q2         The squared momentum transfer.
* @param pdf        The coresponding PDF value.
* @return On success `ENT_RETURN_SUCCESS` is returned otherwise an error
* code is returned as detailed below.
*
* Compute the Parton Distribution Function (PDF) for the given _parton_ and
* kinematic parameters. If *parton* equals `ENT_N_PARTONS`, then the whole PDF
* set is returned, assuming that *pdf* is an array of at least `ENT_N_PARTONS`
* wide.
*
* __Error codes__
*
*     ENT_RETURN_DOMAIN_ERROR     Some input parameter is invalid.
*/
enum ent_return ent_physics_pdf(struct ent_physics * physics,
    enum ent_parton parton, double x, double q2, double * value);

/**
* Compute Structure Functions.
*
* @param physics    A handle for the physics.
* @param projectile The projectile PID.
* @param target     The target PID.
* @param process    The DIS process.
* @param x          The Bjorken _x_ fraction.
* @param q2         The squared momentum transfer.
* @param F2         The F2 SF, or `NULL`.
* @param F3         The F3 SF, or `NULL`.
* @param FL         The FL SF, or `NULL`.
* @return On success `ENT_RETURN_SUCCESS` is returned otherwise an error
* code is returned as detailed below.
*
* Compute the DIS Structure Functions (SFs) for the given projectile, target,
* process and kinematic parameters. Note that target must be one of
* `ENT_PID_NEUTRON` or `ENT_PID_PROTON`.
*
* __Error codes__
*
*     ENT_RETURN_DOMAIN_ERROR     Some input parameter is invalid.
*/
enum ent_return ent_physics_sf(struct ent_physics * physics,
    enum ent_pid projectile, enum ent_pid target, enum ent_process process,
    double x, double Q2, double * F2, double * F3, double * FL);

/**
* Exit events for a neutrino Monte-Carlo transport.
*/
enum ent_event {
        /** No event, e.g. exit when an error occured. */
        ENT_EVENT_NONE = 0,
        /** The neutrino has triggered a backward decay from a muon. */
        ENT_EVENT_DECAY_MUON,
        /** The neutrino has triggered a backward decay from a tau. */
        ENT_EVENT_DECAY_TAU,
        /** The neutrino has exit the simulation area. */
        ENT_EVENT_EXIT,
        /** The neutrino travelled distance has reached a user limit. */
        ENT_EVENT_LIMIT_DISTANCE,
        /** The neutrino travelled grammage has reached a user limit. */
        ENT_EVENT_LIMIT_GRAMMAGE,
        /** The neutrino has interacted. */
        ENT_EVENT_INTERACTION
};

/** Maximum number of collision products. */
#define ENT_PRODUCTS_SIZE 3

/**
 * Simulate a Monte-Carlo collision.
 *
 * @param physics    A handle for the physics.
 * @param context    The Monte-Carlo simulation context.
 * @param state      The initial / final state of the tracked particle.
 * @param medium     The target medium.
 * @param process    The interaction process, if specified.
 * @param products   Additional interaction product(s), or `NULL` if not
 *                   required.
 * @return On success `ENT_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Simulated a Monte-Carlo collision for the given projectile and target. If
 * _process_ is negative the interaction process is randomly selected according
 * to the cross-sections of all possible processes.
 *
 * __Error codes__
 *
 *     ENT_RETURN_DOMAIN_ERROR     Some input parameter is inconsistent.
 */
enum ent_return ent_collide(struct ent_physics * physics,
    struct ent_context * context, struct ent_state * state,
    struct ent_medium * medium, enum ent_process process,
    struct ent_state * products);

/**
 * Perform a Monte-Carlo transport.
 *
 * @param physics    A handle for the physics or `NULL`.
 * @param context    The Monte-Carlo simulation context.
 * @param state      The initial / final state of the tracked particle.
 * @param products   Additional interaction product(s), or `NULL` if not
 *                   required.
 * @param event      The event that ended the Monte-Carlo transport.
 * @return On success `ENT_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Transport the given particle in a medium. The Monte-Carlo transport ends
 * whenever an interaction occurs, or the neutrino escapes, or a user
 * specified limit is reached. Note that if *physics* is `NULL` neutrino
 * interactions are disabled.
 *
 * __Error codes__
 *
 *     ENT_RETURN_DOMAIN_ERROR     Some input parameter is inconsistent or an
 * inconsistent step value was returned.
 */
enum ent_return ent_transport(struct ent_physics * physics,
    struct ent_context * context, struct ent_state * state,
    struct ent_state * products, enum ent_event * event);

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
