CC=gcc
Cflags=-std=c99 -Wall
LDFLAGS=-fopenmp
TARGETS=life_seq life_omp

.PHONY: clean

all: $(TARGETS)

life_omp:life_omp.o

life_seq:life_seq.o

$(TARGETS):
	$(CC) -o $@ $^ $(LDFLAGS)

%.o:%.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o *~ $(TARGETS)
