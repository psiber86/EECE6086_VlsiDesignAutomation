CC = g++
CFLAGS = -O2 -g -Wall

all: klalgo

klalgo: klmain.o
	$(CC) $(CFLAGS) klmain.o -o klalgo

klmain.o: klmain.cpp
	$(CC) -c klmain.cpp

clean:
	rm -rf *.o klalgo
