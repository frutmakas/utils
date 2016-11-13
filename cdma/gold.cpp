/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/04/05 18:49:15 $
 * $Revision: 1.1.2.3 $
 * $Id: gold.cpp,v 1.1.2.3 2004/04/05 18:49:15 syed Exp $
 ********************************************************************/

/* goldcode                       */
/* Nov-10-01, G.R. MOHAMMAD-KHANI, creation  */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "globaldef.h"

#define B_0001 1

#include "tools/utilitis.h"

DVector dGoldCode(int InitialValue[2], int Polynomial[2], int polydeg, int idxwidth, int idx, bool normalize=true){

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
		fprintf(stderr, "\nERROR_Code_2: (idxwidth<polydeg)||(idxwidth>polydeg+1)||(idxwidth>30) \n");
		fprintf(stderr, "\n+ idxwidth : %d\n+ polydeg : %d \n", idxwidth, polydeg);
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
	DVector m = normalize ? n*normalizer : n;
//	nspdbMpy1(normalizer, SequenceC, n.taille);
	delete[] SequenceA;
	delete[] SequenceB;

	return m;	
}

ZVector zGoldCode(int InitialValue[2], int Polynomial[2], int polydeg, int idxwidth, int idx , bool normalize=true){
	return ZVector(dGoldCode(InitialValue, Polynomial, polydeg, idxwidth, idx, normalize));
}
