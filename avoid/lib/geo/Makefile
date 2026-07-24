SRC := $(filter-out $(wildcard src/*_test.c),\
	   $(filter-out $(wildcard src/*_example.c), $(wildcard src/*.c)))
OBJ := $(addprefix obj/,$(notdir $(SRC:.c=.o)))
HDR = include/geo.h 
EXAMPLE := bin/geo_example
TEST:= bin/geo_test

LIB := include/ 

CC=gcc
CFLAGS = -c -Wall -Wextra -O2 -std=c99 $(foreach lib, $(LIB), -I$(lib))
CLINK = -lm 

.PRECIOUS: obj/%.o

all: $(OBJ)

example: $(EXAMPLE)

test: $(TEST)

bin/%: obj/%.o $(OBJ) 
	$(CC) -o $@ $^ $(CLINK)

obj/%.o: src/%.c $(HDR)
	$(CC) $(CFLAGS) -o $@ $<

clean:
	$(foreach dep, $(DEPS), rm -rf ../$(dep)/obj/*)
	rm -rf obj/* bin/*
