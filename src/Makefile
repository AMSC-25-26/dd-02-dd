# ======================================================
# Makefile - MPI Schwarz Solver
# ======================================================

CXX      := mpicxx
CXXFLAGS := -O2 -std=c++20 -Wall -Wextra -pedantic
LDFLAGS  :=

# =========================
# DIRECTORIES
# =========================
SRC_DIR    := .
INC_DIR    := ../include
BUILD_DIR  := ../build
OUTPUT_DIR := ../output

TARGET := $(BUILD_DIR)/schwarz_mpi

SRCS := main.cpp schwarz_solver.cpp sequential_solver.cpp
OBJS := $(SRCS:%.cpp=$(BUILD_DIR)/%.o)

# =========================
# DEFAULT RUN PARAMETERS
# (used ONLY in run_args)
# =========================
NNODES      ?= 200
OVERLAP     ?= 4
MU          ?= 0.01
C           ?= 5.0
A           ?= 0.0
B           ?= 1.0
UA          ?= 0.0
UB          ?= 0.0
MAX_IT      ?= 2000
TOL         ?= 1e-6
RUN_SEQ     ?= 1
FORCING     ?= 1
NPROC       ?= 4

ARGS := $(NNODES) $(OVERLAP) $(MU) $(C) $(A) $(B) \
        $(UA) $(UB) $(MAX_IT) $(TOL) $(RUN_SEQ) $(FORCING)

.PHONY: all clean run run_args run_seq help

# =========================
# BUILD
# =========================
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

# =========================
# RUN TARGETS
# =========================

# ---- INTERACTIVE RUN ----
# No command-line arguments -> terminal input is used
run: $(TARGET)
	@mkdir -p $(OUTPUT_DIR)
	cd $(OUTPUT_DIR) && mpirun -np $(NPROC) $(abspath $(TARGET))

# ---- BATCH RUN WITH ARGS ----
# Parameters passed from Makefile variables
run_args: $(TARGET)
	@mkdir -p $(OUTPUT_DIR)
	cd $(OUTPUT_DIR) && mpirun -np $(NPROC) $(abspath $(TARGET)) $(ARGS)

# ---- SEQUENTIAL ONLY ----
run_seq: $(TARGET)
	@mkdir -p $(OUTPUT_DIR)
	cd $(OUTPUT_DIR) && mpirun -np 1 $(abspath $(TARGET))

# =========================
# CLEAN
# =========================
clean:
	rm -rf $(BUILD_DIR) $(OUTPUT_DIR)

# =========================
# HELP
# =========================
help:
	@echo "Directories:"
	@echo "  include/             -> header files (.hpp)"
	@echo "  build/               -> object files and executable"
	@echo "  output/              -> runtime outputs"
	@echo ""
	@echo "Targets:"
	@echo "  make run                 -> interactive input"
	@echo "  make run_args            -> run with Makefile parameters"
	@echo "  make run_seq             -> sequential interactive run"
