OPTS=-O3 -Wall -Wno-deprecated


all : nrutil.o

nrutil.o :  Makefile ../globaldef.h nrutil.cpp nrutil.h
	g++ -c  $(OPTS)  -I../ nrutil.cpp
