UTILS_PATH = ../tools
UTILS = ${UTILS_PATH}/utilitis.cpp ${UTILS_PATH}/utilitis.h
TOOLS = ${UTILS_PATH}/tools.cpp ${UTILS_PATH}/tools.h
RANDGEN_PATH=../rand
RAND=${RANDGEN_PATH}/randgen.cpp ${RANDGEN_PATH}/randgen.h
OPTS=-O3 -Wall -Wno-deprecated

all : fadingtdma.o

fadingtdma.o :   Makefile ../globaldef.h fadingtdma.cpp fadingtdma.h ${TOOLS}  ${UTILS}
	g++ -c $(OPTS) -I../ fadingtdma.cpp


