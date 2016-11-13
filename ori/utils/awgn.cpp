// qpskawgn.cpp : Defines the entry point for the console application.
//

#include  <math.h>

#include "globaldef.h"
#include "awgn.h"

#ifdef _NSP_INTEL_

DCplx zAWGN_nip(DCplx data) {
	return nspzAdd(data, nspzRandGaus(&zgausStatePtr) );
}

void zAWGN(DCplx *data){
	*data=nspzAdd(*data ,nspzRandGaus(&zgausStatePtr));	
}

DCplx *zbAWGN_nip(DCplx *data, int size){
	DCplx *noise=nspzMalloc(size);
	nspzbRandGaus(&zgausStatePtr, noise, size);
	nspzbAdd2(data, noise, size);
	return noise;
}

void zbAWGN(DCplx *data, int size){
	DCplx *noise=nspzMalloc(size);
	nspzbRandGaus(&zgausStatePtr, noise, size);
	nspzbAdd2(noise, data, size);
	nspFree(noise);
}


double dAWGN_nip(double data) {
	return data+nspdRandGaus(&dgausStatePtr);
}

void dAWGN(double *data){
	(*data)+=nspdRandGaus(&dgausStatePtr);	
}

double *dbAWGN_nip(double *data, int size){
	double *noise=nspdMalloc(size);
	nspdbRandGaus(&dgausStatePtr, noise, size);
	nspdbAdd2(data, noise, size);
	return noise;
}

void dbAWGN(double *data, int size){
	double *noise=nspdMalloc(size);
	nspdbRandGaus(&dgausStatePtr, noise, size);
	nspdbAdd2(noise, data, size);
	nspFree(noise);
}

#else 
#include "utilitis.h"
#include "randgen.h"

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

#endif