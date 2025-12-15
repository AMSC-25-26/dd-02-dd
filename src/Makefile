# ======================================================
# Makefile - MPI Schwarz Solver with Sequential Comparison
# ======================================================

CXX      := mpicxx
CXXFLAGS := -O2 -std=c++20 -Wall -Wextra -pedantic
LDFLAGS  :=
TARGET   := schwarz_mpi

SRCS := main.cpp schwarz_solver.cpp sequential_solver.cpp
OBJS := $(SRCS:.cpp=.o)

# =========================
# DEFAULT RUN PARAMETERS
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

.PHONY: all clean run run_seq help

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
run: $(TARGET)
	mpirun -np $(NPROC) ./$(TARGET) $(ARGS)

run_seq: $(TARGET)
	mpirun -np 1 ./$(TARGET) $(ARGS)

# =========================
# CLEAN
# =========================
clean:
	rm -f $(TARGET) $(OBJS) *.csv

# =========================
# HELP
# =========================
help:
	@echo "Usage examples:"
	@echo "  make run"
	@echo "  make run NPROC=8"
	@echo "  make run FORCING=4"
	@echo "  make run NNODES=500 OVERLAP=6 MU=0.05"
