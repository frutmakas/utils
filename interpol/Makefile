NR_PATH = ../nr
NR=${NR_PATH}/nrutil.cpp ${NR_PATH}/nrutil.h
OPTS=-O3 -Wall -Wno-deprecated

all: interpol.o
	
interpol.o :   Makefile ../globaldef.h interpol.cpp interpol.h
	g++ -c $(OPTS) -I../ interpol.cpp
