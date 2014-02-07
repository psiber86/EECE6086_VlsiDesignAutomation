CC = g++
CFLAGS = -O2 -g -Wall

all: parprog

parprog: klmain.o
	$(CC) $(CFLAGS) klmain.o -o parprog

klmain.o: klmain.cpp
	$(CC) $(CFLAGS) -c klmain.cpp

clean:
	rm -rf *.o parprog
