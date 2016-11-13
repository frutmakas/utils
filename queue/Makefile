OPTS = -O3 -Wall -Wno-deprecated


all: list.o fifo.o

fifo.o:  Makefile ../globaldef.h fifo.cpp fifo.h
	g++ -c $(OPTS) -I../ fifo.cpp

list.o:  Makefile ../globaldef.h list.cpp list.h
	g++ -c $(OPTS) -I../ list.cpp
