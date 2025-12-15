# Makefile - MPI Schwarz Solver with Sequential Comparison

CXX      := mpicxx
CXXFLAGS := -O2 -std=c++20 -Wall -Wextra -pedantic
LDFLAGS  :=
TARGET   := schwarz_mpi

SRCS := main.cpp schwarz_solver.cpp sequential_solver.cpp
OBJS := $(SRCS:.cpp=.o)

.PHONY: all clean run run_seq

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp schwarz_solver.hpp sequential_solver.hpp thomas.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	mpirun -np 4 ./$(TARGET)

run_seq: $(TARGET)
	mpirun -np 1 ./$(TARGET) 200 0 0.01 5.0 2000 1e-6 1

clean:
	rm -f $(TARGET) $(OBJS) solution.csv solution_sequential.csv
