#include <math.h>

#include <stdlib.h>
#include <string.h>

#include "globaldef.h"

#ifdef _NSP_INTEL_

#define nsp_UsesVector
#define nsp_UsesSampleGen
#define nsp_UsesConversion
#define nsp_UsesConvolution

#include <nsp.h>

#else 
#error "THIS LIBRARY NEEDS INTEL SIGNAL PROCESSING PERFORMANCE LIBRARY"
#endif

#include "fadingcdma.h"


// taplist size = delayunit size = nbOfTaps; 
// delayunit is referenced to T=0 in multiple of Tc 
// delayunit is sorted ascending
// delayunit[0]=0;

TZRayleighMultipath zInitRayleighMultipath(short *delayunits, DCplx *taplist, int nbOfTaps){
	
	TZRayleighMultipath ray;
	ray.taplen = delayunits[nbOfTaps-1]+1;
	ray.taps = nspzMalloc(ray.taplen);
	nspzbZero(ray.taps, ray.taplen);

	for(int i=0, j=0;i<nbOfTaps;i++){
		ray.taps[delayunits[i]].re = taplist[i].re;
		ray.taps[delayunits[i]].im = taplist[i].im;		
	}

	return ray;
}

TDRayleighMultipath dInitRayleighMultipath(short *delayunits, double *taplist, int nbOfTaps){
	
	TDRayleighMultipath ray;
	ray.taplen = delayunits[nbOfTaps-1]+1;
	ray.taps = nspdMalloc(ray.taplen);
	nspdbZero(ray.taps, ray.taplen);

	for(int i=0, j=0;i<nbOfTaps;i++){
		ray.taps[delayunits[i]] = taplist[i];
		ray.taps[delayunits[i]] = taplist[i];		
	}

	return ray;
}

void dFreeRayleighMultipath(TDRayleighMultipath ray) {
	free(ray.taps);
}

void zFreeRayleighMultipath(TZRayleighMultipath ray){
	nspFree(ray.taps);
}

DCplx *zRayleighMultipath_nip(TZRayleighMultipath ray, DCplx *data_in, int data_size){
	DCplx *convres = nspzMalloc(ray.taplen+data_size-1);
	DCplx *convfin = nspzMalloc(data_size);
	nspzConv(data_in, data_size, ray.taps ,ray.taplen, convres);
	nspzbCopy(convres, convfin, data_size);
	nspFree(convres);
	return convfin;
}

void zRayleighMultipath(TZRayleighMultipath ray, DCplx *data_in, int data_size){
	DCplx *convres = nspzMalloc(ray.taplen+data_size-1);
	nspzConv(data_in, data_size, ray.taps ,ray.taplen, convres);
	nspzbCopy(convres, data_in, data_size);
	nspFree(convres);
}

double *dRayleighMultipath_nip(TDRayleighMultipath ray, double *data_in, int data_size){
	double *convres = nspdMalloc(ray.taplen+data_size-1);
	double *convfin = nspdMalloc(data_size);
	nspdConv(data_in, data_size, ray.taps ,ray.taplen, convres);
	nspdbCopy(convres, convfin, data_size);
	nspFree(convres);
	return convfin;
}

void dRayleighMultipath(TDRayleighMultipath ray, double *data_in, int data_size){
	double *convres = nspdMalloc(ray.taplen+data_size-1);
	nspdConv(data_in, data_size, ray.taps ,ray.taplen, convres);
	nspdbCopy(convres, data_in, data_size);
	nspFree(convres);
}

DCplx *dRayleighMultipath_nip(TDRayleighMultipath ray, DCplx *data_in, int data_size){
	DCplx *zray=nspzMalloc(ray.taplen);
	double *im = nspdMalloc(ray.taplen);
	nspdbZero(im, ray.taplen);
	nspzb2RealToCplx(ray.taps, im, zray, ray.taplen);
	nspFree(im);

	DCplx *convres = nspzMalloc(ray.taplen+data_size-1);
	
	nspzConv(data_in, data_size,  zray ,ray.taplen, convres);
	nspFree(zray);
	
	DCplx *convfin = nspzMalloc(data_size);
	nspzbCopy(convres, convfin, data_size);
	nspFree(convres);
	
	return convfin;
}

void dRayleighMultipath(TDRayleighMultipath ray, DCplx *data_in, int data_size){
	DCplx *zray=nspzMalloc(ray.taplen);
	double *im = nspdMalloc(ray.taplen);
	nspdbZero(im, ray.taplen);
	nspzb2RealToCplx(ray.taps, im, zray, ray.taplen);
	nspFree(im);

	DCplx *convres = nspzMalloc(ray.taplen+data_size-1);
	
	nspzConv(data_in, data_size,  zray ,ray.taplen, convres);
	nspFree(zray);
	
	nspzbCopy(convres, data_in, data_size);
	nspFree(convres);
	
}
