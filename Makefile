# Compiler (command-line overrides remain supported: make CXX=...)
CXX = mpiicpc
MKLFLAG ?= -mkl

# Flags
CXXFLAGS = -O2 -std=c++17 -I $(EIGEN3_INCLUDE_DIR) -DPFQMC_SCALE_SAFE_UDT -w
PFAPACK_ROOT ?= ./inc/pfapack
LFLAGS = $(PFAPACK_ROOT)/c_interface/libcpfapack.a $(PFAPACK_ROOT)/fortran/libpfapack.a

# Source directories
SRC_DIR = src
INC_DIR = inc
OBJ_DIR = obj
BDIR = bin

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)

# Object files
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

# Compile and link
$(BDIR)/main: $(OBJS) $(OBJ_DIR)/main.o | $(BDIR)
		$(CXX) $(MKLFLAG) $(CXXFLAGS) $(OBJS) $(OBJ_DIR)/main.o -o $@ $(LFLAGS)

# Compile C++ source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
		$(CXX) $(MKLFLAG) $(CXXFLAGS) -c $< -o $@ -I $(INC_DIR)

$(OBJ_DIR)/main.o: main.cpp | $(OBJ_DIR)
		$(CXX) $(MKLFLAG) $(CXXFLAGS) -c main.cpp -o $@ -I $(INC_DIR)

$(OBJ_DIR) $(BDIR):
		mkdir -p $@


.PHONY: clean

clean:
		rm -f $(OBJS) $(OBJ_DIR)/main.o $(BDIR)/main
