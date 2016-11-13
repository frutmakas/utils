#ifndef _PSK_H_
#define _PSK_H_

#include "globaldef.h"
#ifdef _NSP_INTEL_

#define nsp_UsesSampleGen
#define nsp_UsesVector

#include "nsp.h"
DCplx *nspzGRAY_QPSK_mod(short *data, int data_size, int *qpsk_size) ;
short *nspzGRAY_QPSK_demod(DCplx *qpsk_data, int qpsk_size, int *binary_size);

double* nspdBPSK_mod(short *data, int data_size); 
double* nspdBPSK_mod(double *data, int data_size);

short* nspwBPSK_demod(double* data, int data_size);
double* nspdBPSK_demod(double* data, int data_size);

DCplx* nspzBPSK_mod(short *data, int data_size) ;
short* nspzBPSK_demod(DCplx* data, int data_size) ;

DCplx *nspz8PSK_mod(short *data, int data_size, int *psk8_size) ;
short *nspz8PSK_demod(DCplx *psk8_data, int psk8_size, int *binary_size);

#else
#include "utilitis.h"

#endif

void psk8init();

typedef struct {
	DCplx sym[8];
	bool initialized;
} TPSK8 ;

extern TPSK8 psk8;

ZVector z8PSK_mod(const WVector &data);
WVector z8PSK_demod(const ZVector &data);


#endif