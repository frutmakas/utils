CANAL_PATH=../canal
CANAL=$(CANAL_PATH)/canal.cpp $(CANAL_PATH)/canal.cpp \
      $(CANAL_PATH)/awgn.cpp $(CANAL_PATH)/awgn.h
FADING_PATH=../fading
FADING=${FADING_PATH}/fadingtdma.cpp ${FADING_PATH}/fadingtdma.h
CDMA_PATH=../cdma
CDMA=$(CDMA_PATH)/gold.cpp $(CDMA_PATH)/gold.h \
     $(CDMA_PATH)/walsh.cpp $(CDMA_PATH)/walsh.h \
     $(CDMA_PATH)/spread.cpp $(CDMA_PATH)/spread.h
MODU_PATH=../modulation
MODU=$(MODU_PATH)/psk.cpp $(MODU_PATH)/psk.h
UTILS_PATH = ../tools
UTILS = $(UTILS_PATH)/utilitis.cpp $(UTILS_PATH)/utilitis.h
TOOLS = $(UTILS_PATH)/tools.cpp $(UTILS_PATH)/tools.h
INTERLEAVER = ${UTILS_PATH}/interleaver.cpp ${UTILS_PATH}/interleaver.h
EXT_UTILS=${UTILS_PATH}/ext_utilitis.h ${UTILS_PATH}/ext_utilitis.cpp
RANDGEN_PATH=../rand
RAND=$(RANDGEN_PATH)/randgen.cpp $(RANDGEN_PATH)/randgen.h
CODING_PATH=../coding
CODING=$(CODING_PATH)/coding.cpp $(CODING_PATH)/coding.h \
       $(CODING_PATH)/tcm_coding.cpp $(CODING_PATH)/tcm_coding.h
FFT_PATH=../fft
FFT=$(FFT_PATH)/fft.cpp $(FFT_PATH)/fft.h
INTERPOL_PATH=../interpol
INTERPOL=$(INTERPOL_PATH)/interpol.cpp $(INTERPOL_PATH)/interpol.h
BINMAT=${UTILS_PATH}/binmat.h ${UTILS_PATH}/binmat.cpp
NR_PATH=../nr
NR=$(NR_PATH)/nrutil.cpp $(NR_PATH)/nrutil.h
STBC_PATH=../stbc
STBC=$(STBC_PATH)/conversion.cpp $(STBC_PATH)/conversion.h
QUEUE_PATH=../queue
FIFO=$(QUEUE_PATH)/fifo.cpp $(QUEUE_PATH)/fifo.h
LIST=$(QUEUE_PATH)/list.cpp $(QUEUE_PATH)/list.h
BLOCK_PATH=$(CODING_PATH)/block
LDPC=$(BLOCK_PATH)/ldpc.cpp $(BLOCK_PATH)/ldpc.h
SORT_PATH=../sort
SORT=$(SORT_PATH)/quicksort.cpp $(SORT_PATH)/quicksort.h $(SORT_PATH)/heapsort.cpp $(SORT_PATH)/heapsort.h
BINMATRIX=${UTILS_PATH}/binsparse.cpp ${UTILS_PATH}/binsparse.h
PILOT = $(UTILS_PATH)/pilot.cpp $(UTILS_PATH)/pilot.h
OFDM = $(MODU_PATH)/ofdm.h $(MODU_PATH)/ofdm.cpp

SHOULDBE_LIB_FILES = $(CANAL_PATH)/canal.o $(CANAL_PATH)/awgn.o $(CDMA_PATH)/gold.o \
          $(MODU_PATH)/psk.o $(UTILS_PATH)/tools.o $(CODING_PATH)/coding.o \
          $(CODING_PATH)/tcm_coding.o $(FFT_PATH)/fft.o $(INTERPOL_PATH)/interpol.o \
          $(NR_PATH)/nrutil.o $(STBC_PATH)/conversion.o $(RANDGEN_PATH)/randgen.o $(UTILS_PATH)/utilitis.o \
	   $(QUEUE_PATH)/fifo.o $(QUEUE_PATH)/list.o $(BLOCK_PATH)/ldpc.o \
	    $(SORT_PATH)/quicksort.o $(SORT_PATH)/heapsort.o ${UTILS_PATH}/binsparse.o \
	    ${UTILS_PATH}/binmat.o ${UTILS_PATH}/ext_utilitis.o $(UTILS_PATH)/pilot.o $(MODU_PATH)/ofdm.o ${UTILS_PATH}/interleaver.o

LIB_FILES = $(CANAL_PATH)/canal.o \
           $(UTILS_PATH)/tools.o  \
          $(FFT_PATH)/fft.o $(INTERPOL_PATH)/interpol.o \
          $(NR_PATH)/nrutil.o $(STBC_PATH)/conversion.o \
	    $(RANDGEN_PATH)/randgen.o $(UTILS_PATH)/utilitis.o  \
	    $(QUEUE_PATH)/fifo.o $(QUEUE_PATH)/list.o $(BLOCK_PATH)/ldpc.o \
	    $(SORT_PATH)/quicksort.o $(SORT_PATH)/heapsort.o ${UTILS_PATH}/binsparse.o \
	   ${UTILS_PATH}/ext_utilitis.o ${UTILS_PATH}/binmat.o ${MODU_PATH}/psk.o \
	   ${FADING_PATH}/fadingtdma.o $(UTILS_PATH)/pilot.o $(MODU_PATH)/ofdm.o \
	   $(CDMA_PATH)/gold.o $(CDMA_PATH)/walsh.o $(CDMA_PATH)/spread.o ${UTILS_PATH}/interleaver.o

OPTS=-O3 -Wall -Wno-deprecated

#all : test.o test2.o emtest.o

#all : test.o  emtest.o em2.o

#all : mccdma2.o MIMOStbc.o cdmatest.o cdma.o sfbcofdm.o
all : sfbcofdmldpc.o stbcofdmldpc.o
	#g++ -lm -o../test_utilitis $(LIB_FILES) test.o
#	g++ -lm -o../test_ldpc $(LIB_FILES) test2.o
	#g++ -lm -o../emtest $(LIB_FILES) emtest.o
	#g++ -lm -o../em2  $(LIB_FILES) em2.o
	#g++ -lm -o../mimostbc  $(LIB_FILES) MIMOStbc.o
	#g++ -lm -o ../mccdma $(LIB_FILES) mccdma2.o
	#g++ -lm -o ../cdmatest $(LIB_FILES) cdmatest.o
	#g++ -lm -o ../scdma $(LIB_FILES) cdma.o
	g++ -lm -o ../sfbcofdmldpc ${LIB_FILES} sfbcofdmldpc.o
	g++ -lm -o ../stbcofdmldpc ${LIB_FILES} stbcofdmldpc.o

MIMOStbc.o  : Makefile ../globaldef.h MIMOStbc.cpp $(CANAL) $(CDMA) $(MODU) $(TOOLS) $(CODING) $(FFT) $(INTERPOL) $(NR) $(STBC) $(RAND) $(UTILS) $(FIFO) $(LIST) $(BINMATRIX) $(LDPC) $(SORT) ${BINMATH} $(MODU) $(PILOT) $(OFDM)
	g++ -c $(OPTS) -I../ -D__TIMESTAMP__="\" `date` \"" MIMOStbc.cpp

sfbcofdmldpc.o : Makefile ../globaldef.h sfbcofdmldpc.cpp ${CANAL} ${CDMA} ${MODU} ${INTERLEAVER} ${CODING} ${RAND} ${UTILS} ${FIFO} ${LIST} ${BINMATRIX} ${LDPC} ${SORT} ${BINMATH} 
	g++ -c ${OPTS} -I../ -D__TIMESTAMP__="\" `date` \"" sfbcofdmldpc.cpp

stbcofdmldpc.o : Makefile ../globaldef.h stbcofdmldpc.cpp ${CANAL} ${CDMA} ${MODU} ${INTERLEAVER} ${CODING} ${RAND} ${UTILS} ${FIFO} ${LIST} ${BINMATRIX} ${LDPC} ${SORT} ${BINMATH} 
	g++ -c ${OPTS} -I../ -D__TIMESTAMP__="\" `date` \"" stbcofdmldpc.cpp

mccdma2.o : Makefile ../globaldef.h mccdma2.cpp $(CANAL) $(CDMA) $(MODU) $(TOOLS) $(RAND) $(UTILS) $(EXT_UTILS)
	g++ -c $(OPTS) -I../ -D__TIMESTAMP__="\" `date` \"" mccdma2.cpp

cdmatest.o : Makefile ../globaldef.h cdmatest.cpp $(CANAL) $(CDMA) $(MODU) $(TOOLS) $(RAND) $(UTILS) $(EXT_UTILS)
	g++ -c $(OPTS) -I../ -D__TIMESTAMP__="\" `date` \"" cdmatest.cpp
cdma.o : Makefile ../globaldef.h cdma.cpp $(CANAL) $(CDMA) $(MODU) $(TOOLS) $(RAND) $(UTILS) $(EXT_UTILS) $(LDPC)
	g++ -c $(OPTS) -I../ -D__TIMESTAMP__="\" `date` \"" cdma.cpp
test.o : Makefile ../globaldef.h test.cpp $(CANAL) $(CDMA) $(MODU) $(TOOLS) $(CODING) $(FFT) $(INTERPOL) $(NR) $(STBC) $(RAND) $(UTILS) $(FIFO) $(LIST) $(BINMATRIX) $(LDPC) $(SORT) ${BINMATH} $(MODU) $(PILOT) $(OFDM)
	g++ -c $(OPTS) -I../ -D__TIMESTAMP__="\" `date` \"" test.cpp

test2.o : Makefile ../globaldef.h test2.cpp $(CANAL) $(CDMA) $(MODU) $(TOOLS) $(CODING) $(FFT) $(INTERPOL) $(NR) $(STBC) $(RAND) $(UTILS) $(FIFO) $(LIST) $(BINMATRIX) $(LDPC) $(SORT) ${BINMATH} $(MODU) $(PILOT) $(OFDM)
	g++ -c $(OPTS) -I../ test2.cpp

emtest.o :  Makefile ../globaldef.h  emtest.cpp $(CANAL) $(CDMA) $(MODU) $(TOOLS) $(CODING) $(FFT) $(INTERPOL) $(NR) $(STBC) $(RAND) $(UTILS) $(FIFO) $(LIST) $(BINMATRIX) $(LDPC) $(SORT) ${BINMATH} $(MODU) $(PILOT) $(OFDM)
	g++ -c $(OPTS) -I../ emtest.cpp


em2.o: Makefile ../globaldef.h em2.cpp $(CANAL) $(CDMA) $(MODU) $(TOOLS) $(CODING) $(FFT) $(INTERPOL) $(NR) $(STBC) $(RAND) $(UTILS) $(FIFO) $(LIST) $(BINMATRIX) $(LDPC) $(SORT) ${BINMATH} $(MODU) $(PILOT) $(OFDM)
	g++ -c $(OPTS) -I../ em2.cpp


