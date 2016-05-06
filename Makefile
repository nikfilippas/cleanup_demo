CC=gcc
CFLAGS=-Wall -g -O0 -Iinclude -std=c99 -fPIC
CFLAGS+=-I/home/damonge/include
LDFLAGS=-lgsl -lgslcblas   -lm


OBJECTS=src/ccl_core.o src/ccl_utils.o src/ccl_power.o src/ccl_placeholder.o src/ccl_background.o
LIB=lib/libccl.a
DYLIB=lib/libccl.so
INC_CCL=


$(LIB): $(OBJECTS)
	ar rc $(LIB) $(OBJECTS)

$(DYLIB): $(OBJECTS)
	$(CC) -shared -o $(DYLIB) $(OBJECTS) $(LDFLAGS)


tests: tests/ccl_test_utils tests/ccl_test_power

# tests/ccl_test_utils: tests/ccl_test_utils.c $(LIB)
# 	$(CC) -o tests/ccl_test_utils  $(CFLAGS) $(LDFLAGS) -lccl tests/ccl_test_utils.c

tests/% : tests/%.c $(LIB)
	$(CC)  $(CFLAGS) $< -o $@ -Llib -lccl $(LDFLAGS) 

src/%.o: src/%.c
	$(CC) -c $(CFLAGS) $< -o $@


clean:
	rm -rf *.dSYM *.o *.a  test_core_cosmo src/*.o lib/*.a