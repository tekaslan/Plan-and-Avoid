SRC := $(filter-out src/dubins_example.c $(wildcard src/*_cgi.c), $(wildcard src/*.c))
OBJ := $(addprefix obj/,$(notdir $(SRC:.c=.o)))
DEPS := geo traj
HDR = include/dubins.h $(foreach dep, $(DEPS), ../$(dep)/include/$(dep).h)
DEPS_SRC := $(foreach dep, $(DEPS), $(wildcard ../$(dep)/src/*.c))
DEPS_OBJ := $(foreach dep, $(DEPS), $(addprefix ../$(dep)/obj/,$(notdir $(patsubst %.c,%.o,$(filter-out $(wildcard ../$(dep)/src/*_test.c) $(wildcard ../$(dep)/src/*_example.c),$(wildcard ../$(dep)/src/*.c))))))

EXAMPLE := bin/dubins_example

LIB := include/ $(foreach dep, $(DEPS), ../$(dep)/include)

CC=gcc
CFLAGS = -c -std=c99 -Wall -Wextra -g $(foreach lib, $(LIB), -I$(lib))
CLINK = -lm -lcurl -lgsl -lgslcblas

.PRECIOUS: obj/%.o

all: chk_deps $(OBJ)

example: $(EXAMPLE)

chk_deps:
	@$(foreach dep, $(DEPS), test -s ../$(dep) || { echo "Needs $(dep) library for compiling. Run GetDeps.sh! Aborting..."; exit; };)

$(DEPS_OBJ): $(DEPS_SRC)
	$(foreach dep, $(DEPS), make -C ../$(dep);)

bin/%: obj/%.o $(OBJ) $(DEPS_OBJ)
	$(CC) -o $@ $^ $(CLINK)

obj/%.o: src/%.c $(HDR)
	$(CC) $(CFLAGS) -o $@ $<

clean:
	$(foreach dep, $(DEPS), rm -rf ../$(dep)/obj/*)
	rm -rf obj/* bin/*
