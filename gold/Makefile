UTILS_PATH = ../tools
UTILS = ${UTILS_PATH}/utilitis.cpp ${UTILS_PATH}/utilitis.h
TOOLS = ${UTILS_PATH}/tools.cpp ${UTILS_PATH}/tools.h
RANDGEN_PATH=../rand
RAND=${RANDGEN_PATH}/randgen.cpp ${RANDGEN_PATH}/randgen.h
OPTS=-O3 -Wall -Wno-deprecated

all : gold.o 

gold.o : ${UTILS} gold.cpp gold.h
	g++ $(OPTS) -c -I../ gold.cpp


