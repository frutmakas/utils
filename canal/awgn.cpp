// qpskawgn.cpp : Defines the entry point for the console application.
//

#include  <math.h>

#include "globaldef.h"
#include "canal/awgn.h"
#include "tools/utilitis.h"
#include "rand/randgen.h"

DCplx zAWGN_nip(DCplx data) {
	return data+zRandGaus(&zgausStatePtr);
}

void zAWGN(DCplx *data){
	*data+=zRandGaus(&zgausStatePtr);	
}

ZVector zbAWGN_nip(const ZVector &data){
	ZVector noise(data.taille);
	zbRandGaus(&zgausStatePtr, noise);
	return data+noise;
}

void zbAWGN(const ZVector &data){
	ZVector noise(data.taille);
	zbRandGaus(&zgausStatePtr, noise);
	for(register int i=0;i<data.taille;i++)
		data.vect[i]+=noise.vect[i];
}

ZVector zbAWGN_nip(const DVector &data){
	ZVector noise(data.taille);
	zbRandGaus(&zgausStatePtr, noise);
	return noise+data;
}

double dAWGN_nip(double data) {
	return data+dRandGaus(&dgausStatePtr);
}

void dAWGN(double *data){
	(*data)+=dRandGaus(&dgausStatePtr);	
}

DVector dbAWGN_nip(const DVector &data){
	DVector noise(data.taille);
	dbRandGaus(&dgausStatePtr, noise);
	return noise+data;
}

ZVector dbAWGN_nip(const ZVector &data){
	ZVector noise(data.taille);
	zbRandGaus(&zgausStatePtr, noise);
	return data+noise;
}


void dbAWGN(const DVector &data){
	DVector noise(data.taille);
	dbRandGaus(&dgausStatePtr, noise);
	for(register int i=0;i<data.taille;i++)
		data.vect[i]+=noise.vect[i];
}

void dbAWGN(const ZVector &data){
	DVector noise(data.taille);
	dbRandGaus(&dgausStatePtr, noise);
	for(register int i=0;i<data.taille;i++)
		data.vect[i]+=noise.vect[i];
}
