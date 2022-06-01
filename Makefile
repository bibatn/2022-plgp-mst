# defines
CC=gcc
MPICC=mpicc
CXX=g++
MPICXX=mpicxx
CFLAGS= -O3 -Wall -std=gnu99 -openmp
CXXFLAGS= -O3 -Wall
LDFLAGS= -O3 -lrt


TARGET = gen_valid_info validation gen_RMAT gen_RMAT_mpi mst_reference mst_reference_mpi mst gen_random  gen_random_mpi

all: $(TARGET)


# your own implementation, executable must called mst
mst: main.o mst.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

mst_reference_mpi: main.mpi.o graph_tools.mpi.o mst_reference.mpi.o
	$(MPICXX) $^ -o $@ $(LDFLAGS)

graph_tools.mpi.o: graph_tools.cpp
	$(MPICXX) -DUSE_MPI $(CXXFLAGS) -o $@ -c $<


# reference implementation
mst_reference: main.o mst_reference.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

gen_RMAT_mpi: gen_RMAT.mpi.o
	$(MPICXX) $^ -o $@ $(LDFLAGS)

gen_random_mpi: gen_random.mpi.o
	$(MPICXX) $^ -o $@ $(LDFLAGS)



gen_random: gen_random.o
	$(CXX) $^ -o $@ $(LDFLAGS)


gen_RMAT: gen_RMAT.o
	$(CXX) $^ -o $@ $(LDFLAGS)

validation: validation.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

gen_valid_info: graph_tools.o mst_reference.o gen_valid_info.o
	$(CXX) $^ -o $@ $(LDFLAGS)


%.mpi.o: %_mpi.cpp
	$(MPICXX) -DUSE_MPI $(CXXFLAGS) -o $@ -c $<

%.mpi.o: %.c
	$(MPICC) -DUSE_MPI $(CFLAGS) -o $@ -c $<

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -rf *.o $(TARGET)
