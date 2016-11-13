#********************************************************************
#*                        File Information                         #*
#********************************************************************
#* $Name:  $
#* $Author: syed $
#* $Date: 2004/03/21 12:51:39 $
#* $Revision: 1.1.2.14 $
#* $Id: Makefile,v 1.1.2.14 2004/03/21 12:51:39 syed Exp $
#********************************************************************/

UTILS_PATH = ../tools
UTILS = ${UTILS_PATH}/utilitis.cpp ${UTILS_PATH}/utilitis.h
BINMATRIX=${UTILS_PATH}/binsparse.h ${UTILS_PATH}/binsparse.cpp
BINTYPES=${UTILS_PATH}/binmat.cpp ${UTILS_PATH}/binmat.h
TOOLS = ${UTILS_PATH}/tools.cpp ${UTILS_PATH}/tools.h
RANDGEN_PATH=../rand
RAND=${RANDGEN_PATH}/randgen.cpp ${RANDGEN_PATH}/randgen.h
PILOT = $(UTILS_PATH)/pilot.cpp $(UTILS_PATH)/pilot.h
EXT_UTILS=${UTILS_PATH}/ext_utilitis.h ${UTILS_PATH}/ext_utilitis.cpp
INTERLEAVER=${UTILS_PATH}/interleaver.h ${UTILS_PATH}/interleaver.cpp

OPTS=-O3 -Wall -Wno-deprecated

all : interleaver.o ext_utilitis.o utilitis.o tools.o binsparse.o binmat.o pilot.o

ext_utilitis.o: ${EXT_UTILS} Makefile ../globaldef.h
	g++ -c ${OPTS} -I../ ext_utilitis.cpp

interleaver.o: ${INTERLEAVER} Makefile ${TOOLS} ${UTILS} ../globaldef.h
	g++ -c ${OPTS} -I../ interleaver.cpp

binsparse.o : Makefile ../globaldef.h ${BINMATRIX}
	g++ -c $(OPTS) -I../ binsparse.cpp

utilitis.o : ${UTILS}  Makefile ../globaldef.h
	g++ -c $(OPTS) -I../ utilitis.cpp

tools.o:  ${TOOLS} ${UTILS}  Makefile ../globaldef.h
	g++ -c $(OPTS) -I../ tools.cpp

binmat.o: ${BINTYPES} ../globaldef.h  Makefile 
	g++ -c $(OPTS) -I../ binmat.cpp

pilot.o: $(UTILS) $(PILOT) Makefile ../globaldef.h
	g++ -c $(OPTS) -I../ pilot.cpp
