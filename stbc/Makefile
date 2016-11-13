UTILS_PATH = ../tools
UTILS = ${UTILS_PATH}/utilitis.cpp ${UTILS_PATH}/utilitis.h
TOOLS = ${UTILS_PATH}/tools.cpp ${UTILS_PATH}/tools.h
RANDGEN_PATH=../rand
RAND=${RANDGEN_PATH}/randgen.cpp ${RANDGEN_PATH}/randgen.h

OPTS=-Wall -O3 -Wno-deprecated 

all : conversion.o

conversion.o :  Makefile ../globaldef.h conversion.cpp conversion.h  ${UTILS}
	g++ -c $(OPTS) -I../ conversion.cpp


