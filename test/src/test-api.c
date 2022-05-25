/* Unit testing framework */
#include "criterion/criterion.h"

/* The ENT library */
#include "ent.h"


/* Mock function for memory allocations. */
static void * malloc_fails(size_t size) { return NULL; }
static void setup_malloc_fails(void) { ent_malloc = &malloc_fails; }
static void teardown_malloc_fails(void) { ent_malloc = &malloc; }

/* Fixture(s) for API test. */
static ent_handler_cb * default_handler = NULL;

static void setup_api(void)
{
        default_handler = ent_error_handler_get();
        ent_error_handler_set(NULL);
}

static void teardown_api(void)
{
        ent_error_handler_set(default_handler);
        default_handler = NULL;
}

TestSuite(api, .init = setup_api, .fini = teardown_api);

/* Ancestor callback for backward case */
static enum ent_pid expected_ancestor = ENT_PID_NONE;

static double ancestor(struct ent_context * context,
    enum ent_pid ancestor, struct ent_state * daughter)
{
        if (ancestor == expected_ancestor) {
                return 1.;
        } else {
                return 0.;
        }
}

/* API tests */
Test(api, ent_collide)
{
        /* Check collider API */
        struct ent_physics * physics;
        const char * data = "share/ent/CMS11-physics.ent";
        ent_physics_create(&physics, data);

        struct ent_context * context;
        ent_context_create(&context);

        struct ent_state state = {
            .pid = ENT_PID_NU_E,
            .energy = 1E+09,
            .direction = {0., 0., 1.},
            .weight = 1.
        };
        struct ent_medium target = {.Z = 1., .A = 1.};
        struct ent_state products[ENT_PRODUCTS_SIZE] = {0};

        cr_expect(ent_collide(physics, context, &state, &target,
            ENT_PROCESS_DIS_CC_OTHER, products) == ENT_RETURN_SUCCESS);
        cr_expect((state.energy < 1E+09) && (state.energy > 0.));
        cr_expect(state.pid == ENT_PID_ELECTRON);
        cr_expect(products->pid == ENT_PID_HADRON);

        state.pid = ENT_PID_NU_E;
        state.energy = 1E+09;
        cr_expect(ent_collide(physics, context, &state, &target,
            ENT_PROCESS_ELASTIC, NULL) == ENT_RETURN_SUCCESS);
        cr_expect((state.energy < 1E+09) && (state.energy > 0.));
        cr_expect(state.pid == ENT_PID_NU_E);

        context->ancestor = &ancestor;
        expected_ancestor = ENT_PID_NU_E;
        state.pid = ENT_PID_ELECTRON;
        state.energy = 1E+07;
        cr_expect(ent_collide(physics, context, &state, &target,
            ENT_PROCESS_DIS_CC_OTHER, products) == ENT_RETURN_SUCCESS);
        cr_expect(state.energy > 1E+07);
        cr_expect(state.weight != 1.);
        cr_expect(state.pid == expected_ancestor);
        cr_expect(products->pid == ENT_PID_HADRON);

        expected_ancestor = ENT_PID_NU_MU;
        state.pid = ENT_PID_POSITRON;
        state.energy = 1E+07;
        state.weight = 1.;
        cr_expect(ent_collide(physics, context, &state, &target,
            ENT_PROCESS_DIS_CC_TOP, products) == ENT_RETURN_SUCCESS);
        cr_expect(state.energy > 1E+07);
        cr_expect(state.weight != 1.);
        cr_expect(state.pid == expected_ancestor);
        cr_expect(products[0].pid == ENT_PID_MUON);
        cr_expect(products[1].pid == ENT_PID_NU_E);
        cr_expect(products[2].pid == ENT_PID_HADRON);

        expected_ancestor = ENT_PID_NU_E;
        state.pid = ENT_PID_NU_E;
        state.energy = 1E+07;
        state.weight = 1.;
        cr_expect(ent_collide(physics, context, &state, &target,
            ENT_PROCESS_ELASTIC, NULL) == ENT_RETURN_SUCCESS);
        cr_expect(state.energy > 1E+07);
        cr_expect(state.weight != 1.);
        cr_expect(state.pid == expected_ancestor);

        /* Check wrong inputs */
        target.A = 0.;
        cr_expect(ent_collide(physics, context, &state, &target,
            ENT_PROCESS_ELASTIC, NULL) == ENT_RETURN_DOMAIN_ERROR);

        ent_context_destroy(&context);
        ent_physics_destroy(&physics);

        target.A = 1.;
        cr_expect(ent_collide(physics, context, &state, &target,
            ENT_PROCESS_ELASTIC, NULL) == ENT_RETURN_BAD_ADDRESS);
}

Test(api, ent_context)
{
        /* Check initialisation */
        struct ent_context * context = NULL;
        cr_expect(ent_context_create(&context) == ENT_RETURN_SUCCESS);
        cr_expect(context != NULL);
        cr_expect(context->medium == NULL);
        cr_expect(context->random != NULL);
        cr_expect(context->ancestor == NULL);
        cr_expect(context->stepping_action == NULL);
        cr_expect(context->distance_max == 0.);
        cr_expect(context->grammage_max == 0.);

        /* Check random getter and setter */
        unsigned long seed;
        seed = ent_context_random_seed(context);
        cr_expect(ent_context_random_seed(context) == seed);

        seed = 1;
        ent_context_random_set(context, &seed);
        cr_expect(ent_context_random_seed(context) == seed);

        ent_context_random_set(context, NULL);

        /* Check random function */
        const double u = context->random(context);
        cr_expect((u > 0.) && (u < 1.));

        /* Check destructor */
        ent_context_destroy(&context);
        cr_expect(context == NULL);

        /* Check extra destructor call */
        ent_context_destroy(&context);

        /* Check user context case */
        struct ent_context user_context_data = { 0 };
        struct ent_context * user_context = &user_context_data;
        ent_context_destroy(&user_context);

        /* Check the NULL case */
        cr_expect(ent_context_create(NULL) == ENT_RETURN_BAD_ADDRESS);
}

Test(api, ent_context_memory, .init = setup_malloc_fails,
    .fini = teardown_malloc_fails)
{
        /* Check memory error */
        struct ent_context * context;
        cr_expect(ent_context_create(&context) == ENT_RETURN_MEMORY_ERROR);
}

Test(api, ent_error)
{
        /* Check functions getter */
#define CHECK_FUNCTION(FUNC)                                                   \
        cr_expect(strcmp(                                                      \
            ent_error_function((ent_function_t *)&FUNC), #FUNC) == 0)

        CHECK_FUNCTION(ent_collide);
        CHECK_FUNCTION(ent_context_create);
        CHECK_FUNCTION(ent_context_destroy);
        CHECK_FUNCTION(ent_context_random_seed);
        CHECK_FUNCTION(ent_context_random_set);
        CHECK_FUNCTION(ent_error_function);
        CHECK_FUNCTION(ent_error_handler_get);
        CHECK_FUNCTION(ent_error_handler_set);
        CHECK_FUNCTION(ent_error_print);
        CHECK_FUNCTION(ent_error_string);
        CHECK_FUNCTION(ent_physics_create);
        CHECK_FUNCTION(ent_physics_cross_section);
        CHECK_FUNCTION(ent_physics_dcs);
        CHECK_FUNCTION(ent_physics_ddcs);
        CHECK_FUNCTION(ent_physics_destroy);
        CHECK_FUNCTION(ent_physics_dump);
        CHECK_FUNCTION(ent_physics_metadata);
        CHECK_FUNCTION(ent_physics_pdf);
        CHECK_FUNCTION(ent_physics_rescale);
        CHECK_FUNCTION(ent_physics_sf);
        CHECK_FUNCTION(ent_transport);
        CHECK_FUNCTION(ent_version);

        cr_expect(ent_error_function((ent_function_t *)&setup_api) == NULL);

#undef CHECK_FUNCTION

        /* Check getter and setter */
        cr_expect(ent_error_handler_get() == NULL);
        ent_error_handler_set(default_handler);
        cr_expect(ent_error_handler_get() == default_handler);
        ent_error_handler_set(NULL);
        cr_expect(ent_error_handler_get() == NULL);

        /* Check error messages */
        cr_expect(ent_error_string(ENT_RETURN_SUCCESS) != NULL);
        cr_expect(ent_error_string(ENT_RETURN_BAD_ADDRESS) != NULL);
        cr_expect(ent_error_string(ENT_RETURN_DOMAIN_ERROR) != NULL);
        cr_expect(ent_error_string(ENT_RETURN_FORMAT_ERROR) != NULL);
        cr_expect(ent_error_string(ENT_RETURN_IO_ERROR) != NULL);
        cr_expect(ent_error_string(ENT_RETURN_MEMORY_ERROR) != NULL);
        cr_expect(ent_error_string(ENT_RETURN_PATH_ERROR) != NULL);
        cr_expect(ent_error_string(-1) == NULL);

        /* Check printer */
        ent_error_print(NULL, ENT_RETURN_SUCCESS,
            (ent_function_t *)&ent_collide, NULL, NULL);

        FILE * stream = fopen("/dev/null", "w");
        ent_error_print(stream, ENT_RETURN_SUCCESS,
            (ent_function_t *)&ent_collide, NULL, NULL);
        fclose(stream);
}

Test(api, ent_physics)
{
        /* Check constructor */
        const char * data = "share/ent/CMS11-physics.ent";
        struct ent_physics * physics = NULL;
        cr_expect(ent_physics_create(&physics, data) == ENT_RETURN_SUCCESS);
        cr_expect(physics != NULL);

        /* Check meta-data */
        cr_expect(ent_physics_metadata(physics) != NULL);
        cr_expect(ent_physics_metadata(NULL) == NULL);

        /* Check cross-section getter */
        double Z = 0.5, A = 1.;
        double energy = 1E+09; /* GeV */
        double cs;

        cs = 0.;
        cr_expect(ent_physics_cross_section(physics, ENT_PID_NU_E, energy, Z, A,
            ENT_PROCESS_DIS_CC, &cs) == ENT_RETURN_SUCCESS);
        cr_expect(cs > 0.);

        cs = 0.;
        cr_expect(ent_physics_cross_section(physics, ENT_PID_NU_E, energy, Z, A,
            ENT_PROCESS_DIS_NC, &cs) == ENT_RETURN_SUCCESS);
        cr_expect(cs > 0.);

        cs = 0.;
        Z = 1.;
        cr_expect(ent_physics_cross_section(physics, ENT_PID_NU_E, energy, Z, A,
            ENT_PROCESS_ELASTIC, &cs) == ENT_RETURN_SUCCESS);
        cr_expect(cs > 0.);

        cs = 0.;
        cr_expect(ent_physics_cross_section(physics, ENT_PID_NU_MU, energy, Z,
            A, ENT_PROCESS_INVERSE_MUON, &cs) == ENT_RETURN_SUCCESS);
        cr_expect(cs > 0.);

        cs = 0.;
        cr_expect(ent_physics_cross_section(physics, ENT_PID_NU_TAU, energy, Z,
            A, ENT_PROCESS_INVERSE_TAU, &cs) == ENT_RETURN_SUCCESS);
        cr_expect(cs > 0.);

        cs = 0.;
        cr_expect(ent_physics_cross_section(physics, ENT_PID_NU_BAR_E, energy,
            Z, A, ENT_PROCESS_GLASHOW_HADRON, &cs) == ENT_RETURN_SUCCESS);
        cr_expect(cs > 0.);

        cs = 1.;
        cr_expect(ent_physics_cross_section(physics, ENT_PID_MUON, energy,
            Z, A, ENT_PROCESS_DIS_CC, &cs) == ENT_RETURN_DOMAIN_ERROR);
        cr_expect(cs == 0.);

        /* Check rescaling */
        double cs0 = 0., cs1 = 0.;
        Z = 0.;
        cr_expect(ent_physics_cross_section(physics, ENT_PID_NU_E, energy, Z, A,
            ENT_PROCESS_DIS_NC, &cs0) == ENT_RETURN_SUCCESS);
        cr_expect(ent_physics_rescale(physics, NULL) == ENT_RETURN_SUCCESS);
        cr_expect(ent_physics_cross_section(physics, ENT_PID_NU_E, energy, Z, A,
            ENT_PROCESS_DIS_NC, &cs1) == ENT_RETURN_SUCCESS);
        cr_expect(cs0 != cs1);

        const char * cs_file = "share/ent/CMS11-cross-section.txt";
        cr_expect(ent_physics_rescale(physics, cs_file) == ENT_RETURN_SUCCESS);
        cr_expect(ent_physics_cross_section(physics, ENT_PID_NU_E, energy, Z, A,
            ENT_PROCESS_DIS_NC, &cs1) == ENT_RETURN_SUCCESS);
        cr_expect(cs0 == cs1);

        cs_file = "share/ent/TOTO21-cross-section.txt";
        cr_expect(
            ent_physics_rescale(physics, cs_file) == ENT_RETURN_PATH_ERROR);
        cr_expect(ent_physics_cross_section(physics, ENT_PID_NU_E, energy, Z, A,
            ENT_PROCESS_DIS_NC, &cs1) == ENT_RETURN_SUCCESS);
        cr_expect(cs0 == cs1);

        cs_file = "share/ent/CMS11-physics.ent";
        cr_expect(
            ent_physics_rescale(physics, cs_file) == ENT_RETURN_FORMAT_ERROR);
        cr_expect(ent_physics_cross_section(physics, ENT_PID_NU_E, energy, Z, A,
            ENT_PROCESS_DIS_NC, &cs1) == ENT_RETURN_SUCCESS);
        cr_expect(cs0 == cs1);

        /* Check DCS getter */
        double dcs, y;
        Z = 1., y = 1E-03;

        dcs = 0.;
        cr_expect(ent_physics_dcs(physics, ENT_PID_NU_E, energy, Z, A,
            ENT_PROCESS_ELASTIC, y, &dcs) == ENT_RETURN_SUCCESS);
        cr_expect(dcs > 0.);

        dcs = 0.;
        cr_expect(ent_physics_dcs(physics, ENT_PID_NU_MU, energy, Z, A,
            ENT_PROCESS_INVERSE_MUON, y, &dcs) == ENT_RETURN_SUCCESS);
        cr_expect(dcs > 0.);

        dcs = 0.;
        cr_expect(ent_physics_dcs(physics, ENT_PID_NU_TAU, energy, Z, A,
            ENT_PROCESS_INVERSE_TAU, y, &dcs) == ENT_RETURN_SUCCESS);
        cr_expect(dcs > 0.);

        dcs = 0.;
        cr_expect(ent_physics_dcs(physics, ENT_PID_NU_BAR_E, energy, Z, A,
            ENT_PROCESS_GLASHOW_HADRON, y, &dcs) == ENT_RETURN_SUCCESS);
        cr_expect(dcs > 0.);

        dcs = 1.;
        cr_expect(ent_physics_dcs(physics, ENT_PID_NU_BAR_E, energy, Z, A,
            ENT_PROCESS_NONE, y, &dcs) == ENT_RETURN_DOMAIN_ERROR);
        cr_expect(dcs == 0.);

        dcs = 1.;
        cr_expect(ent_physics_dcs(physics, ENT_PID_MUON, energy, Z, A,
            ENT_PROCESS_ELASTIC, y, &dcs) == ENT_RETURN_DOMAIN_ERROR);
        cr_expect(dcs == 0.);

        dcs = 1.;
        cr_expect(ent_physics_dcs(physics, ENT_PID_NU_E, energy, Z, A,
            ENT_PROCESS_ELASTIC, -1., &dcs) == ENT_RETURN_DOMAIN_ERROR);
        cr_expect(dcs == 0.);

        /* Check DIS DDCS getter */
        double ddcs, x = 1E-01;
        Z = 0.5, y = 1E-02;

        ddcs = 0.;
        cr_expect(ent_physics_ddcs(physics, ENT_PID_NU_E, energy, Z, A,
            ENT_PROCESS_DIS_CC, x, y, &ddcs) == ENT_RETURN_SUCCESS);
        cr_expect(ddcs > 0.);

        ddcs = 0.;
        cr_expect(ent_physics_ddcs(physics, ENT_PID_NU_E, energy, Z, A,
            ENT_PROCESS_DIS_CC_OTHER, x, y, &ddcs) == ENT_RETURN_SUCCESS);
        cr_expect(ddcs > 0.);

        ddcs = 0.;
        cr_expect(ent_physics_ddcs(physics, ENT_PID_NU_E, energy, Z, A,
            ENT_PROCESS_DIS_CC_TOP, x, y, &ddcs) == ENT_RETURN_SUCCESS);
        cr_expect(ddcs > 0.);

        ddcs = 0.;
        cr_expect(ent_physics_ddcs(physics, ENT_PID_NU_E, energy, Z, A,
            ENT_PROCESS_DIS_NC, x, y, &ddcs) == ENT_RETURN_SUCCESS);
        cr_expect(ddcs > 0.);

        ddcs = 1.;
        cr_expect(ent_physics_ddcs(physics, ENT_PID_NU_E, energy, Z, A,
            ENT_PROCESS_NONE, x, y, &ddcs) == ENT_RETURN_DOMAIN_ERROR);
        cr_expect(ddcs == 0.);

        ddcs = 1.;
        cr_expect(ent_physics_ddcs(physics, ENT_PID_MUON, energy, Z, A,
            ENT_PROCESS_DIS_NC, x, y, &ddcs) == ENT_RETURN_DOMAIN_ERROR);
        cr_expect(ddcs == 0.);

        ddcs = 1.;
        cr_expect(ent_physics_ddcs(physics, ENT_PID_NU_E, energy, Z, A,
            ENT_PROCESS_DIS_NC, -1., y, &ddcs) == ENT_RETURN_DOMAIN_ERROR);
        cr_expect(ddcs == 0.);

        /* Check SF getter */
        double Q2 = 1E+04;
        double sf[3] = {0., 0., 0.};

        sf[0] = sf[1] = sf[2] = 0.;
        cr_expect(ent_physics_sf(physics, ENT_PID_NU_E, ENT_PID_NEUTRON,
            ENT_PROCESS_DIS_CC, x, Q2, sf, sf + 1, sf + 2) ==
            ENT_RETURN_SUCCESS);
        cr_expect((sf[0] > 0.) && (sf[1] > 0.) && (sf[2] > 0.));

        sf[0] = sf[1] = sf[2] = 0.;
        cr_expect(ent_physics_sf(physics, ENT_PID_NU_E, ENT_PID_PROTON,
            ENT_PROCESS_DIS_CC_OTHER, x, Q2, sf, sf + 1, sf + 2) ==
            ENT_RETURN_SUCCESS);
        cr_expect((sf[0] > 0.) && (sf[1] > 0.) && (sf[2] > 0.));

        sf[0] = sf[1] = sf[2] = 0.;
        cr_expect(ent_physics_sf(physics, ENT_PID_NU_E, ENT_PID_NEUTRON,
            ENT_PROCESS_DIS_CC_TOP, x, Q2, sf, sf + 1, sf + 2) ==
            ENT_RETURN_SUCCESS);
        cr_expect((sf[0] > 0.) && (sf[1] > 0.) && (sf[2] > 0.));

        sf[0] = sf[1] = sf[2] = 0.;
        cr_expect(ent_physics_sf(physics, ENT_PID_NU_E, ENT_PID_PROTON,
            ENT_PROCESS_DIS_NC, x, Q2, sf, sf + 1, sf + 2) ==
            ENT_RETURN_SUCCESS);
        cr_expect((sf[0] > 0.) && (sf[1] > 0.) && (sf[2] > 0.));

        sf[0] = sf[1] = sf[2] = 1.;
        cr_expect(ent_physics_sf(physics, ENT_PID_NU_E, ENT_PID_PROTON,
            ENT_PROCESS_DIS_NC, x, -1., sf, sf + 1, sf + 2) ==
            ENT_RETURN_DOMAIN_ERROR);
        cr_expect((sf[0] == 0.) && (sf[1] == 0.) && (sf[2] == 0.));

        sf[0] = sf[1] = sf[2] = 1.;
        cr_expect(ent_physics_sf(physics, ENT_PID_NU_E, ENT_PID_ELECTRON,
            ENT_PROCESS_DIS_NC, x, Q2, sf, sf + 1, sf + 2) ==
            ENT_RETURN_DOMAIN_ERROR);
        cr_expect((sf[0] == 0.) && (sf[1] == 0.) && (sf[2] == 0.));

        sf[0] = sf[1] = sf[2] = 1.;
        cr_expect(ent_physics_sf(physics, ENT_PID_MUON, ENT_PID_NEUTRON,
            ENT_PROCESS_DIS_CC, x, Q2, sf, sf + 1, sf + 2) ==
            ENT_RETURN_DOMAIN_ERROR);
        cr_expect((sf[0] == 0.) && (sf[1] == 0.) && (sf[2] == 0.));

        /* Check the PDF getter */
        double pdf[ENT_N_PARTONS];

        *pdf = 1.;
        cr_expect(ent_physics_pdf(
            physics, ENT_PARTON_U, x, Q2, pdf) == ENT_RETURN_SUCCESS);
        cr_expect(*pdf == 0.);

        int i;
        for (i = 0; i < ENT_N_PARTONS; i++) pdf[i] = 1.;
        cr_expect(ent_physics_pdf(
            physics, ENT_N_PARTONS, x, Q2, pdf) == ENT_RETURN_SUCCESS);
        for (i = 0; i < ENT_N_PARTONS; i++) {
                cr_expect(pdf[i] == 0.);
        }

        cr_expect(ent_physics_pdf(
            physics, ENT_PARTON_U, x, -1., pdf) == ENT_RETURN_DOMAIN_ERROR);

        cr_expect(ent_physics_pdf(
            physics, 100, x, Q2, pdf) == ENT_RETURN_DOMAIN_ERROR);

        /* Check dump */
        cr_expect(
            ent_physics_dump(physics, NULL, "/dev/null") == ENT_RETURN_SUCCESS);

        const char * meta = ent_physics_metadata(physics);
        cr_expect(
            ent_physics_dump(physics, meta, "/dev/null") == ENT_RETURN_SUCCESS);

        /* Check destructor */
        ent_physics_destroy(&physics);
        cr_expect(physics == NULL);

        /* Check extra destructor call */
        ent_physics_destroy(&physics);

        /* Check the NULL case */
        cr_expect(ent_physics_create(NULL, data) == ENT_RETURN_BAD_ADDRESS);

        /* Check missing path */
        physics = (void *)0x1;
        const char * missing_data = "share/ent/TOTO21-physics.ent";
        cr_expect(ent_physics_create(
            &physics, missing_data) == ENT_RETURN_PATH_ERROR);
        cr_expect(physics == NULL);

        /* Check wrong data */
        physics = (void *)0x1;
        const char * wrong_data = "share/ent/CMS11-cross-section.txt";
        cr_expect(ent_physics_create(
            &physics, wrong_data) == ENT_RETURN_FORMAT_ERROR);
        cr_expect(physics == NULL);
}

Test(api, ent_physics_memory, .init = setup_malloc_fails,
    .fini = teardown_malloc_fails)
{
        /* Check memory error */
        struct ent_physics * physics = NULL;
        const char * data = "share/ent/CMS11-physics.ent";
        cr_expect(
            ent_physics_create(&physics, data) == ENT_RETURN_MEMORY_ERROR);
}

Test(api, ent_version)
{
        /* Test getter */
        int major, minor, patch;
        ent_version(&major, &minor, &patch);

        cr_expect(major == ENT_VERSION_MAJOR);
        cr_expect(minor == ENT_VERSION_MINOR);
        cr_expect(patch == ENT_VERSION_PATCH);

        /* Test the NULL case */
        ent_version(NULL, NULL, NULL);
}
