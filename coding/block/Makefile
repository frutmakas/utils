UTILS_PATH=../../tools
UTILS=$(UTILS_PATH)/utilitis.cpp $(UTILS_PATH)/utilitis.h
BINMATRIX=${UTILS_PATH}/binsparse.h ${UTILS_PATH}/binsparse.h
QUEUE_PATH=../../queue
LIST=${QUEUE_PATH}/list.h ${QUEUE_PATH}/list.cpp
OPTS=-O3 -Wall -Wno-deprecated

all: ldpc.o

ldpc.o : Makefile ../../globaldef.h ${UTILS} ${BINMATRIX} ${LIST} ldpc.cpp ldpc.h
	g++ -c $(OPTS) -I../../ ldpc.cpp

