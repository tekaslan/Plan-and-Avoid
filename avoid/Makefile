CXX = /usr/bin/clang++
CC  = /usr/bin/clang

CPPFLAGS = -Iinclude -Ilib/resolution/include \
		   -I/opt/homebrew/Cellar/nlopt/2.10.0/include \
		   -I/opt/homebrew/Cellar/spatialindex/2.1.0/include/ \
           -Ilib/geo/include \
		   -Ilib/dubins/include \
		   -Ilib/traj/include \
		   -Ilib/math_utils/include \
		   -Ilib/airtraffic/include \
		   -Ilib/census/include \
		   -Ilib/resolution/include \
		   -Ilib/aclm/include \
		   -Ilib/node/include \
		   -Ilib/search/include \
		   -Ilib/priorityQ/include \

CXXFLAGS = -std=c++17 -Wall -Wextra $(SANFLAGS)
CFLAGS   = -std=c11  -Wall -Wextra $(SANFLAGS)

LDFLAGS  = $(SANFLAGS) \
		   -L/opt/homebrew/Cellar/nlopt/2.10.0/lib -lnlopt \
		   -L/opt/homebrew/Cellar/spatialindex/2.1.0/lib -lspatialindex \
		   -L/usr/local/lib -lshp \
		   -Wl,-rpath,/usr/local/lib

TARGET = avoid
OBJDIR = obj

CPP_SRC = src/avoid.cpp \
          lib/resolution/src/advisory.cpp \
		  lib/resolution/src/conflict.cpp \
		  lib/resolution/src/traffic.cpp \
		  lib/census/src/rtree_wrapper.cpp
		  
C_SRC   = lib/geo/src/geo.c \
		  lib/dubins/src/dubins.c \
		  lib/traj/src/traj.c \
		  lib/math_utils/src/math_utils.c \
		  lib/airtraffic/src/airtraffic.c \

CPP_OBJ = $(patsubst %.cpp,$(OBJDIR)/%.o,$(CPP_SRC))
C_OBJ   = $(patsubst %.c,$(OBJDIR)/%.o,$(C_SRC))

all: $(TARGET)

$(TARGET): $(CPP_OBJ) $(C_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJDIR) $(TARGET)

run: $(TARGET)
	./$(TARGET)