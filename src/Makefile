# ======================================================
# Makefile - MPI Schwarz Solver
# ======================================================

CXX      := mpicxx
CXXFLAGS := -O2 -std=c++20 -Wall -Wextra -pedantic
LDFLAGS  :=
TARGET   := schwarz_mpi

SRCS := main.cpp schwarz_solver.cpp sequential_solver.cpp
OBJS := $(SRCS:.cpp=.o)

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

%.o: %.cpp schwarz_solver.hpp sequential_solver.hpp thomas.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# =========================
# RUN TARGETS
# =========================

# ---- INTERACTIVE RUN ----
# No command-line arguments -> terminal input is used
run: $(TARGET)
	mpirun -np $(NPROC) ./$(TARGET)

# ---- BATCH RUN WITH ARGS ----
# Parameters passed from Makefile variables
run_args: $(TARGET)
	mpirun -np $(NPROC) ./$(TARGET) $(ARGS)

# ---- SEQUENTIAL ONLY ----
run_seq: $(TARGET)
	mpirun -np 1 ./$(TARGET)

# =========================
# CLEAN
# =========================
clean:
	rm -f $(TARGET) $(OBJS) *.csv

# =========================
# HELP
# =========================
help:
	@echo "Targets:"
	@echo "  make run                 -> interactive input"
	@echo "  make run_args            -> run with Makefile parameters"
	@echo "  make run_seq             -> sequential interactive run"
	@echo ""
	@echo "Examples:"
	@echo "  make run"
	@echo "  make run_args NNODES=500 MU=0.05 FORCING=4"
	@echo "  make run NPROC=8"
