#ifndef _FADING_CDMA_H_
#define _FADING_CDMA_H_

#include <nsp.h>


struct TZRayleighMultipath {
	DCplx *taps; 
	int taplen;
};

struct TDRayleighMultipath {
	double *taps; 
	int taplen;
};

/* 
	- taplist size = delayunit size = nbOfTaps; 
	- delayunit is referenced to T=0 in multiple of Tc 
	- delayunit is sorted ascending
	- delayunit[0]=0;
	- ex : delayunits = {0, 2, 3, 5, 10 } ; taplist = { 0.9, 0.7, 0.19, 0.02, 0.009 }; nbOfTaps = 5
		   the resulting TXRayleighMultipath will be 
	       TXRayleighMultipath::taps = { 0.9, 0, 0.7, 0.19, 0, 0.02, 0, 0, 0, 0, 0.009 }
		   TXRayleighMultipath::taps = 10

*/
TZRayleighMultipath zInitRayleighMultipath(short *delayunits, DCplx *taplist, int nbOfTaps);

/* 
	- taplist size = delayunit size = nbOfTaps; 
	- delayunit is referenced to T=0 in multiple of Tc 
	- delayunit is sorted ascending
	- delayunit[0]=0;
	- ex : delayunits = {0, 2, 3, 5, 10 } ; taplist = { 0.9, 0.7, 0.19, 0.02, 0.009 }; nbOfTaps = 5
		   the resulting TXRayleighMultipath will be 
	       TXRayleighMultipath::taps = { 0.9, 0, 0.7, 0.19, 0, 0.02, 0, 0, 0, 0, 0.009 }
		   TXRayleighMultipath::taps = 10
*/

TDRayleighMultipath dInitRayleighMultipath(short *delayunits, double *taplist, int nbOfTaps);


void dFreeRayleighMultipath(TDRayleighMultipath ray);

void zFreeRayleighMultipath(TZRayleighMultipath ray);


DCplx *zRayleighMultipath_nip(TZRayleighMultipath ray, DCplx *data_in, int data_size);

void zRayleighMultipath(TZRayleighMultipath ray, DCplx *data_in, int data_size);

double *dRayleighMultipath_nip(TDRayleighMultipath ray, double *data_in, int data_size);

void dRayleighMultipath(TDRayleighMultipath ray, double *data_in, int data_size);

DCplx *dRayleighMultipath_nip(TDRayleighMultipath ray, DCplx *data_in, int data_size);

void dRayleighMultipath(TDRayleighMultipath ray, DCplx *data_in, int data_size);

#endif