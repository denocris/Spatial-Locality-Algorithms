CC = g++
CFLAGS=-I.

.SUFFIXES:

all: main

main: main.x
		./main.x

main.x: main.o
	$(CC) $^ -o $@

%.o: %.cpp
	$(CC) -c $< $(CFLAGS) -o $@  # metti @ all'inizio per togliere warning


clean:
	rm -f *.o *~ *.x *.out *.txt

.PHONY: clean
