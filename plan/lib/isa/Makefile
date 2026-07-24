SRC := $(filter-out $(wildcard src/*_test.c),\
	   $(filter-out $(wildcard src/*_example.c), $(wildcard src/*.c)))
OBJ := $(addprefix obj/,$(notdir $(SRC:.c=.o)))
HDR = include/isa.h 
EXAMPLE := bin/isa_example

LIB := include/ 

CC=gcc
CFLAGS = -c -Wall -Wextra -g -std=c99 $(foreach lib, $(LIB), -I$(lib))
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
