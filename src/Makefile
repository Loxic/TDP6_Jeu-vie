CC=gcc
Cflags=-std=c99 -Wall
LDFLAGS=
TARGETS=life_seq

.PHONY: clean

all: $(TARGETS)

life_seq:life_seq.o

$(TARGETS):
	$(CC) -o $@ $^ $(LDFLAGS)

%.o:%.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o *~ $(TARGETS)
