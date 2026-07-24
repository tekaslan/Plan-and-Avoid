DEP_DIR := geo
DEP_FILES := $(foreach d, $(DEP_DIR),\
             $(filter-out %standalone.o,\
             $(addsuffix .o,\
             $(addprefix ../$d/obj/,\
             $(basename $(notdir $(wildcard ../$d/src/*.c)))))))

C_FILES := $(wildcard src/*.c)
OBJ_FILES := $(addprefix obj/,$(notdir $(C_FILES:.c=.o)))


ILIBS := -I./include $(foreach d, $(DEP_DIR), -I../$d/include)
CCFLAGS := -g -std=c99 -c -Wall $(ILIBS)
LDFLAGS := -lm

FP := bin/fp

CC = gcc

all: $(OBJ_FILES)

fp:  $(OBJ_FILES)
	$(foreach d, $(DEP_DIR), make -C ../$d;)
	$(CC) -o $(FP) $(DEP_FILES) $^ $(LDFLAGS)

obj/%.o: src/%.c
	$(CC) $(CCFLAGS) -c -o $@ $<

clean:
	$(foreach d, $(DEP_DIR), make clean -C ../$d;)
	rm -f obj/*.o bin/*.o 
