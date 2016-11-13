#ifndef _FADING_H_
#define _FADING_H_

#include "globaldef.h"

#ifdef _NSP_INTEL_

double *freqd(long nb_rays, double max_doppler);
double **fgwave(long nb_rays, long filter_length, double symbol_time) ;
void tapsgen(double  symbol_T, int filter_length, long nb_samples, long nb_rays, double fd[], double g[], DCplx **fadtaps);
DCplx **fadingtaps(long max_doppler, double symbol_time, long filter_length, long nb_samples, long nb_rays);

#endif

#ifndef _NSP_INTEL_

#include "utilitis.h"

vector freqd(long nb_rays, double max_doppler);
matrice fgwave(long nb_rays, long filter_length, double symbol_time);
cmplx_mat tapsgen(double  symbol_T, int filter_length, long nb_samples, long nb_rays, double *fd, double **g);
cmplx_mat fadingmain(double  symbol_time, long filter_length, long nb_samples, long nb_rays, double seed, double *fd, double **g);
cmplx_mat fadingtaps(long max_doppler, double symbol_time, long filter_length, long nb_samples, long nb_rays);


#endif

#endif