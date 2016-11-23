CFLAGS := -O2 -std=c99 -pedantic -fPIC -Wall
INCLUDE := -Iinclude
LIBS := -lm

.PHONY: lib examples clean

lib: lib/libant.so
	@rm -f *.o

examples: bin/example-basic

clean:
	@rm -rf bin lib *.o

lib/lib%.so: src/%.c include/%.h
	@mkdir -p lib
	@gcc -o $@ $(CFLAGS) $(INCLUDE) -shared $< $(LIBS)

bin/example-%: examples/example-%.c lib
	@mkdir -p bin
	@gcc -o $@ $(CFLAGS) $(INCLUDE) $< -Llib -lant
