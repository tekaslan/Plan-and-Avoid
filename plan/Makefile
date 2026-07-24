# Compiler variables
CC = /usr/bin/clang
CXX = /usr/bin/clang++

# Locate all directories inside lib/*/include
INCLUDE_DIRS := $(wildcard lib/*/include)

# Flags
CFLAGS = -g -fdiagnostics-color=always \
    -I/opt/homebrew/Cellar/proj/9.5.1/include \
    -I/opt/homebrew/Cellar/spatialindex/2.1.0/include/ \
	-I/opt/homebrew/Cellar/nlopt/2.10.0/include/ \
	$(patsubst %, -I%, $(INCLUDE_DIRS))

LDFLAGS = -L/usr/local/lib -lshp \
	-L/opt/homebrew/Cellar/proj/9.5.1/lib -lproj \
	-Wl,-rpath,/usr/local/lib \
	-L/opt/homebrew/Cellar/spatialindex/2.1.0/lib -lspatialindex \
	-L/opt/homebrew/Cellar/nlopt/2.10.0/lib -lnlopt \
	-stdlib=libc++ \
	-L/usr/lib -lstdc++

CXXFLAGS = -std=c++11 -stdlib=libc++ $(CFLAGS)

# Library source files
LIB_SRC = \
	lib/census/src/census.c \
	lib/dubins/src/dubins.c \
	lib/geo/src/geo.c \
	lib/isa/src/isa.c \
	lib/math_utils/src/math_utils.c \
	lib/node/src/node.c \
	lib/priorityQ/src/priorityQ.c \
    lib/search/src/search.c \
	lib/traj/src/traj.c \
	lib/aclm/src/aclm.c \
	lib/airtraffic/src/airtraffic.c
    
# C++ source for the wrapper
CXX_SRC = lib/census/src/rtree_wrapper.cpp

# Target executable name
TARGET = plan

# Object files
OBJ = obj/plan.o \
	lib/census/obj/census.o \
	lib/dubins/obj/dubins.o \
	lib/geo/obj/geo.o \
	lib/isa/obj/isa.o \
	lib/math_utils/obj/math_utils.o \
	lib/node/obj/node.o \
	lib/priorityQ/obj/priorityQ.o \
	lib/search/obj/search.o \
    lib/traj/obj/traj.o \
	lib/aclm/obj/aclm.o \
	lib/census/obj/rtree_wrapper.o \
	lib/airtraffic/obj/airtraffic.o
    
all: $(TARGET)

# Compile the main file
obj/plan.o: src/plan.c
	$(CC) $(CFLAGS) -c $< -o $@
	
lib/census/obj/census.o: lib/census/src/census.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/dubins/obj/dubins.o: lib/dubins/src/dubins.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/geo/obj/geo.o: lib/geo/src/geo.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/isa/obj/isa.o: lib/isa/src/isa.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/math_utils/obj/math_utils.o: lib/math_utils/src/math_utils.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/node/obj/node.o: lib/node/src/node.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/priorityQ/obj/priorityQ.o: lib/priorityQ/src/priorityQ.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/search/obj/search.o: lib/search/src/search.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/traj/obj/traj.o: lib/traj/src/traj.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/aclm/obj/aclm.o: lib/aclm/src/aclm.c
	$(CC) $(CFLAGS) -c $< -o $@

lib/airtraffic/obj/airtraffic.o: lib/airtraffic/src/airtraffic.c
	$(CC) $(CFLAGS) -c $< -o $@

# Compile C++ wrapper file to lib/census/obj/rtree_wrapper.o
lib/census/obj/rtree_wrapper.o: lib/census/src/rtree_wrapper.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link all object files into the final executable
$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $(TARGET) $(LDFLAGS)

# Clean the build
clean:
	rm -f $(OBJ) $(TARGET)