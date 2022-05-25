PREFIX= $(PWD)

CFLAGS=  -O2 -std=c99 -pedantic -fPIC -Wall
INCLUDE= -Iinclude
LIBS=    -lm


# ENT library
VERSION_MAJOR= $(shell grep ENT_VERSION_MAJOR include/ent.h | cut -d' ' -f3)
VERSION_MINOR= $(shell grep ENT_VERSION_MINOR include/ent.h | cut -d' ' -f3)
VERSION_PATCH= $(shell grep ENT_VERSION_PATCH include/ent.h | cut -d' ' -f3)

LIB_NAME=      libent.so
LIB_SHORTNAME= $(LIB_NAME).$(VERSION_MAJOR)
LIB_FULLNAME=  $(LIB_SHORTNAME).$(VERSION_MINOR).$(VERSION_PATCH)

.PHONY: lib
lib: $(PREFIX)/lib/$(LIB_FULLNAME) \
     $(PREFIX)/lib/$(LIB_SHORTNAME) \
     $(PREFIX)/lib/$(LIB_NAME)

$(PREFIX)/lib/$(LIB_FULLNAME): src/ent.c include/ent.h | $(PREFIX)/lib
	$(CC) -o $@ $(CFLAGS) $(INCLUDE) -shared $< $(LIBS)

$(PREFIX)/lib/$(LIB_SHORTNAME): $(PREFIX)/lib/$(LIB_FULLNAME)
	@ln -fs $(LIB_FULLNAME) $@

$(PREFIX)/lib/$(LIB_NAME): $(PREFIX)/lib/$(LIB_SHORTNAME)
	@ln -fs $(LIB_SHORTNAME) $@

$(PREFIX)/lib:
	@mkdir -p $@


# Examples
EXAMPLE_LIBS= $(LIBS) -L$(PREFIX)/lib -Wl,-rpath,$(PREFIX)/lib -lent

.PHONY: examples
examples: $(PREFIX)bin/example-collide \
          $(PREFIX)/bin/example-physics \
          $(PREFIX)/bin/example-transport

$(PREFIX)/bin/example-%: examples/example-%.c | $(PREFIX)/bin lib
	$(CC) -o $@ $(CFLAGS) $(INCLUDE) $< $(EXAMPLE_LIBS)

$(PREFIX)/bin:
	@mkdir -p $@


# Unit tests
CR_PREFIX=test/criterion
CR_INCLUDE=$(CR_PREFIX)/include
CR_LIB=$(CR_PREFIX)/lib

TEST_CFLAGS= -O0 -g3 -std=c99 -pedantic -Wall -D ENT_MALLOC=mock_malloc
TEST_INCLUDE= $(INCLUDE) -I$(CR_INCLUDE)
TEST_LIBS= $(LIBS) -L$(CR_LIB) -Wl,-rpath,$(CR_LIB) -lcriterion

.PHONY: test
test: $(PREFIX)/bin/test-api
	@$<

$(PREFIX)/bin/test-api: test/src/test-api.c src/ent.c include/ent.h | $(PREFIX)/bin
	$(CC) -o $@ $(TEST_CFLAGS) $(TEST_INCLUDE) test/src/test-api.c src/ent.c $(TEST_LIBS)


# Coverage
COVERAGE_CFLAGS= $(TEST_CFLAGS) -fprofile-arcs -ftest-coverage

.PHONY: coverage
coverage: $(PREFIX)/bin/coverage-api
	@lcov --directory . --zerocounters
	@$<
	@lcov --directory . --capture --output-file coverage.info
	@lcov --remove coverage.info '*/test/*' '*/include/*' '/usr/*' --output-file coverage.info
	@lcov --list coverage.info
	@genhtml -o docs/coverage coverage.info
	@rm -f *.gcda *.gcno coverage.info

$(PREFIX)/bin/coverage-api: test/src/test-api.c src/ent.c include/ent.h | $(PREFIX)/bin
	$(CC) -o $@ $(COVERAGE_CFLAGS) $(TEST_INCLUDE) test/src/test-api.c src/ent.c $(TEST_LIBS)


# Cleanup
.PHONY: clean
clean:
	@rm -rf $(PREFIX)/bin $(PREFIX)/lib *.o *.gcda *.gcno coverage.info docs/coverage
