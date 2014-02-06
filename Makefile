CC = g++
CFLAGS = -O2 -g -Wall -lpthread
CORES = $(shell grep -c processor /proc/cpuinfo)

all: parprog

parprog: klmain.o
	$(CC) $(CFLAGS) klmain.o -o parprog

klmain.o: klmain.cpp
	$(CC) -DCORES=$(CORES) -c klmain.cpp

clean:
	rm -rf *.o parprog
