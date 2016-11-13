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
		* Syed Mohamad : (26 march 2002) 
			- ported from Matlab source code
			- optimized for INTEL SPLIB

******************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "globaldef.h"

#ifndef _NSP_INTEL_
#include "utilitis.h"
#endif


#ifdef _NSP_INTEL_

#include <nsp.h>
#include <nsprand.h>

#include  "nspconv.h"
#include  "nspcnv2d.h"
#include  "nspfir2.h"


#include  "nspcvrt.h"

#include "nspdct.h"

#include  "nspfirl.h"
#include  "nspfirh.h"

#include  "nspfirg.h"

#include  "nspiirl.h"
#include  "nspiirh.h"

#include  "nsplmsl.h"
#include  "nsplmsh.h"

#include  "nspmed.h"

#include  "nspmisc.h"

#include  "nsprand.h"
#include  "nsptone.h"
#include  "nsptrngl.h"

#include  "nspfft.h"
#include  "nspgrtzl.h"
#include  "nspgrtzw.h"

#include  "nsparith.h"
#include  "nspcorr.h"
#include  "nsplaw.h"
#include  "nsplnexp.h"
#include  "nspsampl.h"
#include  "nspfirh.h"
#include  "nspfirl.h"
#include  "nspfirg.h"
#include  "nsprsmpl.h"
#include  "nspvec.h"
#include  "nspdotp.h"
#include  "nspnorm.h"
#include  "nsplogic.h"
#include  "nspdiv.h"
#include  "nspatan.h"

#include  "nspwin.h"

#include  "nspwlt.h"


#include "nsperror.h"
#include "nspalloc.h"
#endif

#include "fadingtdma.h"

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

extern 	NSPDRandUniState tapstatePtr;


#ifndef _NSP_INTEL_
double uni_rand();

long int seed;

vector freqd(long nb_rays, double max_doppler){
	vector fd=init_vect(nb_rays,0);
//	double *fd=(double*)malloc(sizeof(double)*nb_rays);
	for(int i = 0; i <nb_rays;i++)
		fd.vect[i]=max_doppler*cos(uni_rand());
	return fd;

}

/* symbol_time is in µsec*/
matrice fgwave(long nb_rays, long filter_length, double symbol_time) {
	double symbol_T=symbol_time*1e-6;
	double Tmax=10e-6 ,b_tau=1e-6;
	double *xtau = (double*)malloc(sizeof(double)*nb_rays);
	double exp_tmax_btau = 1.0-exp(-Tmax/b_tau);
	double tapfactor;
	matrice fg;
	fg.mat=(double**)malloc(sizeof(double*)*nb_rays);
	fg.lin=nb_rays;
	fg.col=filter_length;

	srand(15424);

	for(int i = 0; i < nb_rays;i++) {
		double randfactor=((double)rand())/RAND_MAX;

		xtau[i]=-b_tau*log10(1-randfactor*exp_tmax_btau);
		fg.mat[i]=(double*)malloc(sizeof(double)*filter_length);
		for(int k=0;k<filter_length;k++){
			tapfactor=fabs((k*symbol_T-xtau[i])/symbol_T);
			if (tapfactor<1.0) 
				fg.mat[i][k]=1-tapfactor;
			else
				fg.mat[i][k]=0;
		}
	}
	free(xtau);
	return fg;
}

cmplx_mat tapsgen(double  symbol_T, int filter_length, long nb_samples, long nb_rays, vector fd, matrice g){
	int i,j,k;
	cmplx sample={0,0};
	symbol_T*=1e-6;

	cmplx_mat fadtaps;
	fadtaps.mat=(cmplx**)malloc(sizeof(cmplx*)*nb_samples);
	fadtaps.col=nb_samples;
	fadtaps.lin=filter_length;

	for(i = 0; i<nb_samples;i++)
		fadtaps.mat[i]=(cmplx*)malloc(sizeof(cmplx)*filter_length);

	/* Generate the coefficients of tap-delay-line */

	for(j=0; j<nb_samples; j++){
		for(k=0; k<filter_length; k++){
			sample.Re=0; sample.Im=0;
			for(i=0;i<nb_rays;i++){
				double tau=2*M_PI*fd.vect[i]*j*symbol_T+uni_rand();
				long phi=i+k*nb_rays;
				sample.Re+=cos(tau)*g.mat[i][k];
				sample.Im+=sin(tau)*g.mat[i][k];
			}
			double coef=sqrt(2)/sqrt(nb_rays);
			sample.Re*=coef;
			sample.Im*=coef;
			
			fadtaps.mat[j][k].Re=sample.Re;
			fadtaps.mat[j][k].Im=sample.Im;
		}
	}
	return fadtaps;
}

cmplx_mat fadingmain(double  symbol_time, long filter_length, long nb_samples, long nb_rays, double seed, vector fd, matrice g){
	return tapsgen(symbol_time, filter_length, nb_samples, nb_rays, fd, g);
}

cmplx_mat fadingtaps(long max_doppler, double symbol_time, long filter_length, long nb_samples, long nb_rays){
	vector fd=freqd(nb_rays, max_doppler);
	matrice g=fgwave(nb_rays,filter_length,symbol_time);
	cmplx_mat tap;
	tap=fadingmain(symbol_time, filter_length, nb_samples, nb_rays, 12345, fd, g);
	liberer(fd);
	/*for (int i=0;i<nb_rays;i++) 
		free(g[i]);*/
	liberer(g);
	return tap;
}
/* Uniformly distributed random variables generator */
double uni_rand() /*0->1*/
{
	double r;
	seed=(25733*seed+13849)%65536;
	r=(double)seed/65536;
	return r;
}

#endif

#ifdef _NSP_INTEL_

double *freqd(long nb_rays, double max_doppler){
	double *fd=nspdMalloc(nb_rays);
	nspdbRandUni(&tapstatePtr, fd, nb_rays);
	for(int i = 0; i <nb_rays;i++)
		fd[i]=max_doppler*cos(fd[i]);
	return fd;
}

/* symbol_time is in µsec*/
double **fgwave(long nb_rays, long filter_length, double symbol_time) {
	double symbol_T=symbol_time*1e-6;
	double Tmax=10e-6 ,b_tau=1e-6;
	double *xtau = nspdMalloc(nb_rays);
	double exp_tmax_btau = 1.0-exp(-Tmax/b_tau);
	double tapfactor;
	double **fg=(double**)malloc(sizeof(double*)*nb_rays);
	srand(15424);
	for(int i = 0; i < nb_rays;i++) {
		double randfactor=(((double)rand())/RAND_MAX);
		xtau[i]=(-b_tau*log10(1-randfactor*exp_tmax_btau))/symbol_T;
		fg[i]=nspdMalloc(filter_length);
		for(int k=0;k<filter_length;k++){
			tapfactor=fabs(k-xtau[i]);
			if (tapfactor<1.0) 
				fg[i][k]=1-tapfactor;
			else
				fg[i][k]=0;
		}
	}
	nspFree(xtau);
	return fg;
}


DCplx **tapsgen(double  symbol_T, int filter_length, long nb_samples, long nb_rays, double *fd, double **g){
	register int i,j,k;
	double theta[max_rays];

	DCplx sample={0.0,0.0};
	symbol_T*=1e-6;

	/* Generate the random variables */
	nspdbRandUni(&tapstatePtr, theta, max_rays);

	DCplx **fadtaps=(DCplx**)malloc(sizeof(DCplx*)*nb_samples);
	for(i = 0; i<nb_samples;i++)
		fadtaps[i]=nspzMalloc(filter_length);

	/* Generate the coefficients of tap-delay-line */
	double NSP_2PI_symbol_T=NSP_2PI*symbol_T; // optimisation
	
	DCplx coef={sqrt(2)/sqrt(nb_rays),0.0};
	for(j=0; j<nb_samples; j++){
		double xtmp=NSP_2PI_symbol_T*j; // optimisation
		for(k=0; k<filter_length; k++){
			sample.re=0.0; sample.im=0.0;
			for(i=0;i<nb_rays;i++){
				register double tau=xtmp*fd[i]+theta[i];
				sample.re+=cos(tau)*g[i][k];
				sample.im+=sin(tau)*g[i][k];
			
			}
			fadtaps[j][k]=nspzMpy(sample, coef);
		}
	}
	return fadtaps;
}


DCplx **fadingmain(double  symbol_time, long filter_length, long nb_samples, long nb_rays, double seed, double *fd, double **g){
	return tapsgen(symbol_time, filter_length, nb_samples, nb_rays, fd, g);
}

DCplx **fadingtaps(long max_doppler, double symbol_time, long filter_length, long nb_samples, long nb_rays){
	double *fd=freqd(nb_rays, max_doppler);
	double **g=fgwave(nb_rays,filter_length,symbol_time);
	DCplx **tap;
	tap=fadingmain(symbol_time, filter_length, nb_samples, nb_rays, 12345, fd, g);
	nspFree(fd);
	for (int i=0;i<nb_rays;i++) 
		nspFree(g[i]);
	free(g);
	return tap;
}

#endif


