CC=gcc
MPICC=mpicc
CFLAGS=-std=c99 -Wall -fopenmp
LDFLAGS=-fopenmp
TARGETS=life_seq life_omp
MPI_TARGETS=  life_mpi_sync life_mpi_async life_mpi_pers
.PHONY: clean

all: $(TARGETS) $(MPI_TARGETS)


life_omp:life_omp.o

life_seq:life_seq.o

life_mpi_sync:
	$(MPICC) -std=c99 -Wall life_mpi_sync.c -o life_mpi_sync -lm
life_mpi_async:
	$(MPICC) -std=c99 -Wall life_mpi_async.c -o life_mpi_async -lm
life_mpi_pers:
	$(MPICC) -std=c99 -Wall life_mpi_pers.c -o life_mpi_pers -lm


$(TARGETS):
	$(CC) -o $@ $^ $(LDFLAGS)

%.o:%.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o *~ $(TARGETS) $(MPI_TARGETS)
