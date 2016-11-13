/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/06/06 09:23:20 $
 * $Revision: 1.1.2.21 $
 * $Id: fadingtdma.cpp,v 1.1.2.21 2006/06/06 09:23:20 syed Exp $
 ********************************************************************/
 
/******************************************************************

Fading channel generator, generating a number L of fading taps and 
outputing a mathematica list. The program generates complex taps 
for a fixed receiver filter model. It takes as inputs the

max_doppler		: maximum doppler frequency
symbol_T		: the symbol time in micro seconds
filter_length	: The filter tap number
nb_samples		: The number of samples to generate
nb_rays			: The number of simulated rays 
seed			: A random seed
fd				: Doppler frequency spectrum
g				: signal waveform

The current channel model is an exponential delay power spectrum
p_tau = EXP(-tau/b_tau); b_tau is a system constant with the Jakes 
doppler spectrum, and uniform angle distribution. 

Program is called in the matlab by

tap2=fadingmain(symbol_time,filter_length,nb_samples,nb_rays,seed,fd,g)';

The output is a matrix of complex sample values for all the taps.
The (2*k-1)th and the (2*k)th column is the real and imaginary part
for the kth tap respectively, where k=1,2....

Author: * Tao Wu       : Version 3.0 (August, 1998)
		* Syed Mohamad : 
			- (26 march  2002)  ported from Matlab source code
			- (29 August 2002)  optimized for ENSIL UTILITIS
			- (Dec 2003)  Optimized for speed
			- (March 2004) Fix multipath problem ..

******************************************************************/
//#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 


/*#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>*/
#include <iostream>
#include <fstream>
using namespace std;

#ifdef WIN_DOS
#pragma intrinsic(sin, cos, fabs, sqrt, pow, exp, log10)
#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ )
#endif
#include "globaldef.h"

#ifndef _NSP_INTEL_
#include "tools/utilitis.h"
#else 
#error "This library is meant to be used with the utils_nonsp branch"
#endif

#include "fadingtdma.h"
#include "rand/randgen.h"

#define pi (4 * atan(1.0))
#define max_rays 1000
#define max_taps 10
#define msymbol_T prhs[0]
#define mfilter_length prhs[1]
#define mnb_samples [2]
#define mnb_rays prhs[3]
#define mseed prhs[4]
#define mfd prhs[5]
#define mg prhs[6]
#define TAPS plhs[0]

#ifndef _NSP_INTEL_

dRandUniStatePtr tappy;

DVector freqd(long nb_rays, double max_doppler, dRandUniStatePtr *tapseed){
	DVector fd(nb_rays,0);
	dbRandUni(tapseed, fd);
	for(register int i = 0; i <nb_rays;i++)
		fd.vect[i]=max_doppler*cos(UTILITIS_M_2_PI*fd.vect[i]);
	return fd;
}

/* symbol_time is in µsec*/
DMatrix fgwave(long nb_rays, long filter_length, double symbol_time, dRandUniStatePtr *tapseed) {
	double symbol_T=symbol_time*1e-6;
	double symbol_T_1=1.0/symbol_T;
	double Tmax=10e-6 ,b_tau=1e-6;
	double *xtau = (double*)malloc(sizeof(double)*nb_rays);
	register double exp_tmax_btau = 1.0-exp(-Tmax/b_tau);
	register double tapfactor;
	DMatrix fg(nb_rays, filter_length);

	DVector randfactor(nb_rays);

	dbRandUni(tapseed, randfactor);

	for(register int i = 0; i < nb_rays;i++) {
		//double randfactor=dRandUni(tapseed); 
		xtau[i]=-b_tau*log10(1.0-randfactor.vect[i]*exp_tmax_btau);
		for(register int k=0;k<filter_length;k++){
			tapfactor=fabs(k-xtau[i]*symbol_T_1);
			fg.mat[i][k] = tapfactor < 1.0 ? 1.0-tapfactor : 0.0;
		}
	}
	free(xtau);
	return fg;
}

cmplx_mat tapsgen(double  symbol_T, int filter_length, long nb_samples, long nb_rays, const DVector &fd, const DMatrix &g, dRandUniStatePtr *tapseed){
	int i,j,k;
	cmplx sample=DCplx(0,0);
	symbol_T*=1e-6;
	cmplx_mat fadtaps(nb_samples, filter_length);

	/* Generate the coefficients of tap-delay-line */
	double coef=sqrt(2.0/nb_rays);
	DVector theta(nb_rays);

/*	for(int t=0;t<theta.taille;t++) {
		theta.vect[t]=M_2_PI*dRandUni(tapseed);
	}
	*/
	dbRandUni(tapseed, theta);

	DVector norm(filter_length);
	for(j=0; j<nb_samples; j++){
		double symbol_T_j=symbol_T*j;
		for(k=0; k<filter_length; k++){
			sample.re=sample.im=0.0;
			for(i=0;i<nb_rays;i++){				
				//double tau=theta.vect[i]+M_2_PI*fd.vect[i]*j*symbol_T; 
				double tau=UTILITIS_M_2_PI*(theta.vect[i]+fd.vect[i]*symbol_T_j); 
				sample.re+=cos(tau)*g.mat[i][k];
				sample.im+=sin(tau)*g.mat[i][k];
			}
			sample*=coef;
			fadtaps.mat[j][k]=sample;
			//norm.vect[k]+=sample.re*sample.re+sample.im*sample.im;
		}
	}

// normalization
/*	double nb_samples_1 = 1.0/nb_samples;
	for(register int kk=0;kk<filter_length;kk++) {
		norm.vect[kk]=1.0/sqrt(norm.vect[kk]*nb_samples_1);
		for(register int jj=0;jj<nb_samples;jj++) {
			fadtaps.mat[jj][kk] *= norm.vect[kk];
		}
	}
*/
	return fadtaps;
}

cmplx_mat fadingmain(double  symbol_time, long filter_length, long nb_samples, long nb_rays, double seed, const DVector &fd, const DMatrix &g, dRandUniStatePtr *tapseed){
	return tapsgen(symbol_time, filter_length, nb_samples, nb_rays, fd, g, tapseed);
}

cmplx_mat fadingtaps(double max_doppler, double symbol_time, long filter_length, long nb_samples, long nb_rays, dRandUniStatePtr *tapseed, double fd_fluct) { // tapseed must be initialized to uni[0;1]{
	if(fd_fluct>50 || fd_fluct< 0) {
		throw CFadingTapsException(__FILE__, __LINE__, FD_FLUCT_OUT_OF_RANGE);//'fd_fluctuation_percentage value must be between 0 and 50');
	}

	double delta=fd_fluct/100.0;
	int extra_samples = int(max_doppler*10);
	int samples = nb_samples+extra_samples;
	DVector doppler(filter_length);
	ZMatrix taps(nb_samples, filter_length);
	for (int l=0;l<filter_length;l++) { //ll=1:filter_length
		register int r;
		doppler.vect[l] = max_doppler +(dRandUni(tapseed)*max_doppler*delta*2-delta);
		//freqd(nb_rays,doppler.vect[l]);
		//fgwave(nb_rays,filter_length,symbol_time);
		ZMatrix tmp = tapsgen(symbol_time, 1, samples, nb_rays, freqd(nb_rays,doppler.vect[l], tapseed), 
							    fgwave(nb_rays,filter_length,symbol_time, tapseed),  tapseed);
		double sum=0;
		for (r=0;r<nb_samples;r++) {
			taps.mat[r][l]=tmp.mat[r+extra_samples][0];
			sum+=taps.mat[r][l].mod2();
		}
		double coef = sqrt(nb_samples/sum);
		for (r=0;r<nb_samples;r++) {
			taps.mat[r][l]=taps.mat[r][l]*coef;
		}
	}

	return taps;
}

cmplx_mat fadingtaps(double max_doppler, double symbol_time, long filter_length, long nb_samples, long nb_rays, dRandUniStatePtr *tapseed, const DVector &tapsamplitude, double fd_fluct) { // tapseed must be initialized to uni[0;1]{
	if(fd_fluct>50 || fd_fluct< 0) {
		throw CFadingTapsException(__FILE__, __LINE__, FD_FLUCT_OUT_OF_RANGE);//'fd_fluctuation_percentage value must be between 0 and 50');
	}
	if (tapsamplitude.taille!=filter_length) {
		throw CFadingTapsException(__FILE__, __LINE__, TAPS_AMPLITUDE_SIZE_ERROR);//'fd_fluctuation_percentage value must be between 0 and 50');
	}

	double delta=fd_fluct/100.0;
	int extra_samples = int(max_doppler*10);
	int samples = nb_samples+extra_samples;
	DVector doppler(filter_length);
	ZMatrix taps(nb_samples, filter_length);
	for (int l=0;l<filter_length;l++) { //ll=1:filter_length
		register int r;
		doppler.vect[l] = max_doppler +(dRandUni(tapseed)*max_doppler*delta*2.0-delta);
		//freqd(nb_rays,doppler.vect[l]);
		//fgwave(nb_rays,filter_length,symbol_time);
		ZMatrix tmp = tapsgen(symbol_time, 1, samples, nb_rays, freqd(nb_rays,doppler.vect[l], tapseed), 
							    fgwave(nb_rays,filter_length,symbol_time, tapseed),  tapseed);
		double sum=0;
		for (r=0;r<nb_samples;r++) {
			taps.mat[r][l]=tmp.mat[r+extra_samples][0];
			sum+=taps.mat[r][l].mod2();
		}
		double coef = sqrt(nb_samples/sum)*tapsamplitude.vect[l];
		for (r=0;r<nb_samples;r++) {
			taps.mat[r][l]=taps.mat[r][l]*coef;
		}
	}

	return taps;
}

#endif
