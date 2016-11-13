UTILS_PATH = ../tools
UTILS = ${UTILS_PATH}/utilitis.cpp ${UTILS_PATH}/utilitis.h
TOOLS = ${UTILS_PATH}/tools.cpp ${UTILS_PATH}/tools.h
RANDGEN_PATH=../rand
RAND=${RANDGEN_PATH}/randgen.cpp ${RANDGEN_PATH}/randgen.h
OPTS = -O3 -Wall -Wno-deprecated

all : fft.o 

fft.o:  Makefile ../globaldef.h fft.cpp fft.h ${UTILS} ${TOOLS}
	g++ -c $(OPTS) -I../ fft.cpp


