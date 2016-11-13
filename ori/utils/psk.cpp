// qpskaw gn.cpp : Defines the entry point for the console application.
//
#include  <math.h>
#include  <stdio.h>
#include "psk.h"
#include "nspdcplx.h"
#include "nspmatrice.h"
#include "tools.h"

#ifdef _NSP_INTEL_
#define nsp_UsesConversion
#include "nsp.h"

#define INV_SQRT2 0,70711



#define DEBUG
DCplx *nspz8PSK_mod(short *data, int data_size, int *psk8_size){
	if (!psk8.initialized) {
		fprintf(stderr , "INITIALIZE PSK8 FIRST !!! call psk8init()");
		return NULL;
	}
	*psk8_size = (int)ceil(data_size/3.0);
	DCplx *data_out=nspzMalloc(*psk8_size);
	for(register int i = 0, j=0; i < data_size; i+=3, j++) {
#ifndef DEBUG
		data_out[j]=psk8.sym[bin2oct(&data[i], data_size-i)];
#else	 //DEBUG 
		int oct = bin2oct(&data[i], data_size-i);
		data_out[j]=psk8.sym[oct];
#endif
	}
	return data_out;
}

short *nspz8PSK_demod(DCplx *psk8_data, int psk8_size, int *binary_size){
	if (!psk8.initialized) {
		fprintf(stderr , "INITIALIZE PSK8 FIRST !!! call psk8init()");
		return NULL;
	}
	*binary_size=psk8_size*3;
	short *data_out=nspwMalloc(*binary_size);
	for(int i  =0 ,j=0; i < psk8_size; i++, j+=3) {
		int min_pos=0;
		double min_distance = dcplxmag(psk8_data[i]-psk8.sym[0]);
		double dist = 0.0;
		for(int k = 1; k < 8 ;k++) {
			dist = dcplxmag(psk8_data[i]-psk8.sym[k]);
			if (dist < min_distance) {
				min_distance=dist; min_pos=k;
			}
		}
		oct2bin(min_pos, &(data_out[j]));

	}
	return data_out;
}


/*  Modulation QPSK shifted, codage Gray  */
DCplx *nspzGRAY_QPSK_mod(short *data, int data_size, int *qpsk_size) {
	int i;
	double amp=sqrt(2.0)/2.0;
	short *cdata=data;
	*qpsk_size=data_size/2+data_size%2;
	
	DCplx *qpsk = nspzMalloc(*qpsk_size);
	
	for (i=0; i< *qpsk_size; i++){
		qpsk[i].re=(*(cdata++)==0) ? amp : -amp;
		qpsk[i].im=(*(cdata++)==0) ? amp : -amp;
	}
	return qpsk;
}

/*-------------------------------------------------------*/
/*  Demodulation QPSK shifted, codage Gray  */
short *nspzGRAY_QPSK_demod(DCplx *qpsk_data, int qpsk_size, int *binary_size){
	int i;
	*binary_size=qpsk_size*2;
	short *binary=nspwMalloc(*binary_size);
	short *xbin=binary;

	for (i=0; i< *binary_size; i++){
		*(xbin++) = (qpsk_data[i].re >= 0.0) ? 0 : 1;
		*(xbin++) = (qpsk_data[i].im >= 0.0) ? 0 : 1;
	}
	return binary;
}

/*----------------------------------------*/
double* nspdBPSK_mod(short *data, int data_size) {
	double *b=nspdMalloc(data_size);
	for (int i=0 ; i<data_size; i++)
		b[i]=data[i]*2.0-1.0;
		
	return b;
}


double* nspdBPSK_mod(double *data, int data_size) {
	double *b=nspdMalloc(data_size);
	for (int i=0 ; i<data_size; i++)
		b[i]=data[i]<0?-1.0:1.0;
	return b;
}

/*----------------------------------------*/
short* nspwBPSK_demod(double* data, int data_size) {
	short *b=nspwMalloc(data_size);
   	for (int i=0 ; i<data_size ; i++)
		b[i]=(data[i]>0?1:0);
	return b;
}

/*----------------------------------------*/
double* nspdBPSK_demod(double* data, int data_size) {
	double *b=nspdMalloc(data_size);
   	for (int i=0 ; i<data_size ; i++)
		b[i]=(data[i]>0?1.0:0.0);
	return b;
}
/*----------------------------------------*/
DCplx* nspzBPSK_mod(short *data, int data_size) {
	DCplx *b=nspzMalloc(data_size);
	for (int i=0 ; i<data_size; i++){
		b[i].re=data[i]*2.0-1.0;
		b[i].im=0.0;
	}
	return b;
}

/*----------------------------------------*/
short* nspzBPSK_demod(DCplx* data, int data_size) {
	short *b=nspwMalloc(data_size);
	double *d =nspdMalloc(data_size);
	nspzbPhase(data, d, data_size);
   	for (int i=0 ; i<data_size ; i++)
		b[i]=d[i]>0.0?1:0;
	nspFree(d);
	return b;
}

#else //_NSP_INTEL_

ZVector z8PSK_mod(const WVector &data){
	if (!psk8.initialized) {
		fprintf(stderr , "INITIALIZE PSK8 FIRST !!! call psk8init()");
		return NULL;
	}
	ZVector data_out((int)ceil(data.taille/3.0));
//	*psk8_size = (int)ceil(data.taille/3.0);
//	DCplx *data_out=nspzMalloc(*psk8_size);
	for(register int i = 0, j=0; i < data.taille; i+=3, j++) {
#ifndef DEBUG
		data_out.vect[j]=psk8.sym[bin2oct(&data.vect[i], data.taille-i)];
#else	
		int oct = bin2oct(&data.vect[i], data.taille-i);
		data_out.vect[j]=psk8.sym[oct];
#endif
	}
	return data_out;
}

WVector z8PSK_demod(const ZVector &data){
	if (!psk8.initialized) {
		fprintf(stderr , "INITIALIZE PSK8 FIRST !!! call psk8init()");
		return NULL;
	}

	WVector data_out(data.taille*3);

//	*binary_size=data.taille*3;
//	short *data_out=nspwMalloc(*binary_size);

	for(int i  =0 ,j=0; i < data.taille; i++, j+=3) {
		int min_pos=0;
		double min_distance = (data.vect[i]-psk8.sym[0]).mag();
		double dist = 0.0;
		for(int k = 1; k < 8 ;k++) {
			dist = (data.vect[i]-psk8.sym[k]).mag();
			if (dist < min_distance) {
				min_distance=dist; min_pos=k;
			}
		}
		oct2bin(min_pos, &(data_out.vect[j]));
	}

	return data_out;
}

#endif

void psk8init(){
	psk8.initialized = true;

    psk8.sym[0].re = 1.0;
	psk8.sym[0].im = 0.0;
    psk8.sym[1].re = 0.70711;
	psk8.sym[1].im = 0.70711;
    psk8.sym[2].re = 0.0;
	psk8.sym[2].im = 1.0;
    psk8.sym[3].re = -0.70711;
	psk8.sym[3].im = 0.70711;

    psk8.sym[4].re = -1.0;
	psk8.sym[4].im = 0.0;
    psk8.sym[5].re = -0.70711;
	psk8.sym[5].im = -0.70711;
    psk8.sym[6].re = 0.0;
	psk8.sym[6].im = -1.0;
    psk8.sym[7].re = 0.70711;
	psk8.sym[7].im = -0.70711;

#ifdef DEBUG
	for(int i = 0; i < 8; i++) 
		print(psk8.sym[i]);
#endif
}