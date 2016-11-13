/* goldcode                       */
/* Nov-10-01, G.R. MOHAMMAD-KHANI, creation  */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "globaldef.h"

#define B_0001 1

#ifdef _NSP_INTEL_

#define nsp_UsesVector
#define nsp_UsesConversion
#define nsp_UsesSampleGen

#include <nsp.h>



double *dGoldCode(int InitialValue[2], int Polynomial[2], int polydeg, int idxwidth, int idx, int *code_size){

//int InitialValue[2]={31,31};
//int Polynomial[2]={18, 27};

	
	int seqsize=(B_0001 << polydeg) - 1;
	*code_size=seqsize;

	int *SequenceA=(int*)malloc(sizeof(int)*seqsize);
	int *SequenceB=(int*)malloc(sizeof(int)*seqsize);
	double *SequenceC=(double *)malloc(sizeof(double)*seqsize);
	
	int BitMask = (B_0001 << polydeg) - 1;
	unsigned int LFSR, Temp, Result;
	int BlockLength = (B_0001 << polydeg) - 1;
	int     *p, k=0;  
	double *pc;
	int    *pSeq0, *pSeq1;
	int    UIndex;
	int    IndexMask   = (B_0001 << idxwidth) - 1;
	int    *pLast = SequenceB + BlockLength;
	
	/* Error Messages of Parametrisation */
	if((polydeg < 2) || (polydeg > 30)) {
		fprintf(stderr, "\nERROR_Code_1: \n");
		free(SequenceA);
		free(SequenceB);
		free(SequenceC);
		return NULL;
	}

	if((idxwidth<polydeg)||(idxwidth>polydeg+1)||(idxwidth>30)) {
		fprintf(stderr, "\nERROR_Code_2: \n");
		free(SequenceA);
		free(SequenceB);
		free(SequenceC);
		return NULL;
	}
	
	if(((idxwidth == polydeg) & (idx < 0 || idx >=BlockLength)) || ((idx == (polydeg + 1)) & (idx < 0 || idx > BlockLength+1))){
		fprintf(stderr, "\nERROR_Code_3: \n");
		free(SequenceA);
		free(SequenceB);
		free(SequenceC);
		return NULL;
	}
	
	/* generate first m-sequence */
	LFSR  =  (unsigned int) InitialValue[0]; //<<<<------
	p     = SequenceA; 

	int i = BlockLength;
	
	while (i-->0){
		k++;
		Temp = (LFSR & Polynomial[0]); 
		Result = 0;
		int j=polydeg;
		while(j-->0){
			Result ^= (Temp & 1);
			Temp >>= 1;
		}
		*p++ = Result;
		LFSR = ((LFSR << 1) | Result) & BitMask;
		
	}

	/* generate second m-sequence */
	LFSR  = (unsigned int) InitialValue[1]; 
	p     = SequenceB; 
	i=BlockLength;
	while(i-->0){	  
		Temp = (LFSR & Polynomial[1]);
		Result = 0;
		int j=polydeg;
		while(j-->0){
			Result ^= (Temp & 1);
			Temp >>= 1;
		}
		*p++ = Result;
		LFSR = ((LFSR << 1) | Result) & BitMask;
	}

	/* generate gold_sequence */
	pc      = SequenceC;
	UIndex = idx & IndexMask;
	if (idx == BlockLength)
	{
		pSeq0 = SequenceA;
		i=BlockLength;
		while(i-->0) {
			*pc++ = (*pSeq0++)*2.0-1.0;
		}
	}
	else if (idx == BlockLength + 1)
	{
		pSeq1 = SequenceB;
		i=BlockLength;
		while(i-->0){
			*pc++ = (*pSeq1++)*2.0-1.0;
		}
	}
	else
	{
		pSeq0 = SequenceA;
		pSeq1 = SequenceB - UIndex + BlockLength;
		i=BlockLength;
		while(i-->0){
			if (pSeq1 == pLast) pSeq1 = SequenceB; 
			*pc++ = (*pSeq0++ ^ *pSeq1++)*2.0-1.0;
		}
	}
	double normalizer=1.0/sqrt(*code_size);
	nspdbMpy1(normalizer, SequenceC, *code_size);
	free(SequenceA);
	free(SequenceB);

	return SequenceC;	
}

DCplx *zGoldCode(int InitialValue[2], int Polynomial[2], int polydeg, int idxwidth, int idx, int *code_size){
	double *sequence = dGoldCode(InitialValue, Polynomial, polydeg, idxwidth, idx, code_size);
	double *ImZero = (double*)calloc(*code_size, sizeof(double));
	DCplx *zSeq = nspzMalloc(*code_size);
	nspzb2RealToCplx(sequence, ImZero,zSeq, *code_size);
	free(sequence);
	free(ImZero);
	return zSeq;
}

DCplx *zSpreadData(DCplx *data, int data_size, DCplx *spreadcode, int spread_size, int *spreaddata_size){
	*spreaddata_size=spread_size*data_size;
	DCplx *data_out=nspzMalloc(*spreaddata_size);
	DCplx *codecopy = nspzMalloc(spread_size);
	
	for(int i=0;i<data_size;i++) {
		nspzbCopy(spreadcode, codecopy, spread_size);
		nspzbMpy1(data[i], codecopy, spread_size);
		nspzbCopy(codecopy, &(data_out[i*spread_size]), spread_size);
	}

	nspFree(codecopy);
	return data_out;
}

DCplx *zSpreadData(DCplx *data, int data_size, double *spreadcode, int spread_size, int *spreaddata_size){
	double *ImZero = (double*)calloc(spread_size, sizeof(double));
	DCplx *zSeq = nspzMalloc(spread_size);
	nspzb2RealToCplx(spreadcode, ImZero, zSeq, spread_size);	
	DCplx *data_out=zSpreadData(data, data_size, spreadcode, spread_size, spreaddata_size);
	nspFree(ImZero);
	nspFree(zSeq);
	return data_out;
}

double *dSpreadData(double *data, int data_size, double *spreadcode, int spread_size, int *spreaddata_size){
	*spreaddata_size=spread_size*data_size;
	double *data_out=nspdMalloc(*spreaddata_size);
	double *codecopy = nspdMalloc(spread_size);
	
	for(int i=0;i<data_size;i++) {
		nspdbCopy(spreadcode, codecopy, spread_size);
		nspdbMpy1(data[i], codecopy, spread_size);
		nspdbCopy(codecopy, &(data_out[i*spread_size]), spread_size);
	}

	nspFree(codecopy);
	return data_out;
}

double *dUnspreadData(double *data, int data_size, double *spreadcode, int spread_size, int *unspreaded_size){
	*unspreaded_size = (int) ceil(data_size/spread_size);
	double *data_out = (double*)malloc(sizeof(double)*(*unspreaded_size));
	for(int i=0, k=0; i< *unspreaded_size;i++) {
		double sum = 0.0;
		for(int j=0;j<spread_size && k < data_size; j++, k++) {
			sum += data[k]*spreadcode[j];
		}
		data_out[i]=sum;
	}
	return data_out;
}

#ifdef GOLD_DEBUG

void main(int argc, char** argv){
	int code_size;

//	int InitialValue[2]={31,31};
//	int Polynomial[2]={18, 27};
	NSPWRandUniState statePtr;
	double *code=dGoldCode(InitialValue, Polynomial ,3,3,0, &code_size);
	int binary_len=32;
	short *wbinary=nspwMalloc(binary_len);
	double *binary=nspdMalloc(binary_len);

	int spreaddata_size=0;
	nspwRandUniInit(rand(), 0, 2, &statePtr);

	nspwbRandUni(&statePtr, wbinary, binary_len);
	nspdbIntToFloat(wbinary, binary, binary_len,sizeof(short)*8, NSP_Noflags);
	nspFree(wbinary);

	for(int i =0; i < binary_len;i++) binary[i]=binary[i]*2.0-1.0;
	double *coded=dSpreadData(binary, binary_len, code, code_size, &spreaddata_size);

	FILE *f=fopen("spread.txt", "wt");
	for(i=0;i<binary_len;i++){
		fprintf(f, "\n%1.1f X ", binary[i]);
		for(int j=0;j<code_size;j++) fprintf(f, "%1.1f ", code[j]);
		fprintf(f," --> ");
		for(j=0;j<code_size;j++) fprintf(f, "%1.1f ", coded[i*code_size+j]);
	}
	fclose(f);
	nspFree(coded);
	nspFree(code);
	nspFree(binary);
	
}
#endif


#else //_NSP_INTEL_
#include "utilitis.h"

DVector dGoldCode(int InitialValue[2], int Polynomial[2], int polydeg, int idxwidth, int idx){

//int InitialValue[2]={31,31};
//int Polynomial[2]={18, 27};

	
	int seqsize=(B_0001 << polydeg) - 1;
	DVector n;
	n.taille=seqsize;

	int *SequenceA= new int[seqsize];// (int*)malloc(sizeof(int)*seqsize);
	int *SequenceB= new int[seqsize];// (int*)malloc(sizeof(int)*seqsize);
	double *SequenceC= new double[seqsize];// (double *)malloc(sizeof(double)*seqsize);
	
	int BitMask = (B_0001 << polydeg) - 1;
	unsigned int LFSR, Temp, Result;
	int BlockLength = (B_0001 << polydeg) - 1;
	int     *p, k=0;  
	double *pc;
	int    *pSeq0, *pSeq1;
	int    UIndex;
	int    IndexMask   = (B_0001 << idxwidth) - 1;
	int    *pLast = SequenceB + BlockLength;
	
	/* Error Messages of Parametrisation */
	if((polydeg < 2) || (polydeg > 30)) {
		fprintf(stderr, "\nERROR_Code_1: \n");
		delete[] SequenceA;
		delete[] SequenceB;
		delete[] SequenceC;
		return DVector();
	}

	if((idxwidth<polydeg)||(idxwidth>polydeg+1)||(idxwidth>30)) {
		fprintf(stderr, "\nERROR_Code_2: \n");
		delete[] SequenceA;
		delete[] SequenceB;
		delete[] SequenceC;
		return DVector();
	}
	
	if(((idxwidth == polydeg) & (idx < 0 || idx >=BlockLength)) || ((idx == (polydeg + 1)) & (idx < 0 || idx > BlockLength+1))){
		fprintf(stderr, "\nERROR_Code_3: \n");
		delete[] SequenceA;
		delete[] SequenceB;
		delete[] SequenceC;
		return DVector();
	}
	
	/* generate first m-sequence */
	LFSR  =  (unsigned int) InitialValue[0]; //<<<<------
	p     = SequenceA; 

	int i = BlockLength;
	
	while (i-->0){
		k++;
		Temp = (LFSR & Polynomial[0]); 
		Result = 0;
		int j=polydeg;
		while(j-->0){
			Result ^= (Temp & 1);
			Temp >>= 1;
		}
		*p++ = Result;
		LFSR = ((LFSR << 1) | Result) & BitMask;
		
	}

	/* generate second m-sequence */
	LFSR  = (unsigned int) InitialValue[1]; 
	p     = SequenceB; 
	i=BlockLength;
	while(i-->0){	  
		Temp = (LFSR & Polynomial[1]);
		Result = 0;
		int j=polydeg;
		while(j-->0){
			Result ^= (Temp & 1);
			Temp >>= 1;
		}
		*p++ = Result;
		LFSR = ((LFSR << 1) | Result) & BitMask;
	}

	/* generate gold_sequence */
	pc      = SequenceC;
	UIndex = idx & IndexMask;
	if (idx == BlockLength)
	{
		pSeq0 = SequenceA;
		i=BlockLength;
		while(i-->0) {
			*pc++ = (*pSeq0++)*2.0-1.0;
		}
	}
	else if (idx == BlockLength + 1)
	{
		pSeq1 = SequenceB;
		i=BlockLength;
		while(i-->0){
			*pc++ = (*pSeq1++)*2.0-1.0;
		}
	}
	else
	{
		pSeq0 = SequenceA;
		pSeq1 = SequenceB - UIndex + BlockLength;
		i=BlockLength;
		while(i-->0){
			if (pSeq1 == pLast) pSeq1 = SequenceB; 
			*pc++ = (*pSeq0++ ^ *pSeq1++)*2.0-1.0;
		}
	}
	double normalizer=1.0/sqrt(n.taille);
	n.vect=SequenceC;
	DVector m=n*normalizer;
//	nspdbMpy1(normalizer, SequenceC, n.taille);
	delete[] SequenceA;
	delete[] SequenceB;

	return m;	
}

ZVector zGoldCode(int InitialValue[2], int Polynomial[2], int polydeg, int idxwidth, int idx){
	return ZVector(dGoldCode(InitialValue, Polynomial, polydeg, idxwidth, idx));
}

ZVector SpreadData(const ZVector &data, const ZVector &spreadcode){
	ZVector spreaded(data.taille*spreadcode.taille);

	for(int i=0;i<data.taille;i++) {
		spreaded.insert(spreadcode*data.vect[i], i*spreadcode.taille);
	}
	return spreaded;
}

ZVector SpreadData(const ZVector &data, const DVector &spreadcode){
	ZVector spreaded(data.taille*spreadcode.taille);

	for(int i=0;i<data.taille;i++) {
		spreaded.insert(spreadcode*data.vect[i], i*spreadcode.taille);
	}
	return spreaded;
}

ZVector SpreadData(const DVector &data, const ZVector &spreadcode){
	ZVector spreaded(data.taille*spreadcode.taille);

	for(int i=0;i<data.taille;i++) {
		spreaded.insert(spreadcode*data.vect[i], i*spreadcode.taille);
	}
	return spreaded;
}

DVector zSpreadData(const DVector &data, const DVector &spreadcode){
	DVector spreaded(data.taille*spreadcode.taille);

	for(int i=0;i<data.taille;i++) {
		spreaded.insert(spreadcode*data.vect[i], i*spreadcode.taille);
	}
	return spreaded;
}

DVector UnspreadData(const DVector &data, const DVector &spreadcode);


#endif