// fadingtaps.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "globaldef.h"
#include <math.h>
#include <stdlib.h>


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
#include "fading.h"
#include "psk.h"
#include "awgn.h"


#define FILTER_LENGTH	3
#define NB_SAMPLES		100
#define MAX_DOPPLER		2700
#define NB_RAYS			500
#define SYMBOL_TIME		3.7


#define MOYENNE 0.0
//#define VARIANCE 0.3
#define IN_SIZE 100
#define VARIANCE_MIN 0
#define VARIANCE_MAX 1.0
//#define VARIANCE_STEP 0.1

#define NB_ITERATION 100
#define VARIANCE_STUDY
//#define TAPTEST

NSPZRandGausState zgausStatePtr;
NSPDRandGausState dgausStatePtr;
NSPWRandUniState  uniStatePtr;
NSPDRandUniState tapstatePtr;


int main(int argc, char* argv[])
{

#ifdef _NSP_INTEL_
	double EBN0=0.0001;
	double VARIANCE=0.0;
	

#ifndef VARIANCE_STUDY
#else 
	double VARIANCE_STEP=1.0;
	printf("\nEnter AWGN VARIANCE step : ");
	scanf("%lf", &VARIANCE_STEP );

	
#endif

//	VARIANCE = SYMBOL_TIME*1.0e-6/(2.0*pow(10.0,EBN0/10));
		
	nspzRandGausInit(314159, MOYENNE, VARIANCE, &zgausStatePtr);
	nspdRandGausInit(445755, MOYENNE, VARIANCE, &dgausStatePtr);
	nspdRandUniInit(1560, 0, NSP_2PI, &tapstatePtr);
	nspwRandUniInit (114521, 0, 2, &uniStatePtr);

#ifdef TAPTEST	
	FILE *f=fopen("diagram.txt", "wt");
	DCplx *tofft=nspzMalloc(512);
	int unlocked=1;

	for(int w=0;w<100;w++) {
		DCplx **taps = fadingtaps(MAX_DOPPLER,SYMBOL_TIME,FILTER_LENGTH,NB_SAMPLES,NB_RAYS);
		for(int x=0;x<NB_SAMPLES;x++){
			for(int y=0;y<FILTER_LENGTH;y++)
				fprintf(f,"%.3f;%.3f;;;", taps[x][y].re, taps[x][y].im);
			fprintf(f,"\n");
		}
		printf("# ");
	}

	fclose(f);
	double *d=freqd(512, MAX_DOPPLER);
	FILE *f3=fopen("raytime.txt", "wt");
	for(int n=0;n<512;n++) {
		fprintf(f3, "\n%5f", d[n]);
	}
	fclose(f3);
	
	nspdRealFft(d, 9,NSP_Forw);
	FILE *f2=fopen("rayfft.txt", "wt");
	for(n=0;n<512;n++) {
		fprintf(f2, "\n%5f", d[n]);
	}
	fclose(f2);
	nspFree(d);

#else
#ifndef QPSK
	
	printf("\n-**************** ENTERING BPSK TEST SYSTEM ********************-\n\n");
	printf("\n\a*** Multipath Rayleigh Fading Channel Parameter ***\n\nMAX_DOPPLER : %d\nSYMBOL_TIME : %f\nFILTER_LENGTH : %d\nNB_SAMPLES : %d\nNB_RAYS : %d", MAX_DOPPLER,SYMBOL_TIME,FILTER_LENGTH,NB_SAMPLES,NB_RAYS);
#ifndef VARIANCE_STUDY
	printf("\n\n\n*** Gaussian Random Parameter ***\n\nMOYENNE : %f\nVARIANCE : %f", MOYENNE, VARIANCE);
#else
	printf("\n\n\n*** Gaussian Random Parameter ***\n\nMOYENNE : %f\nCalculating TEB for VARIANCE varying from %f to %f with step %f", MOYENNE, VARIANCE_MIN, VARIANCE_MAX, VARIANCE_STEP);
#endif

	printf("\n\n*****************************************\n");
	printf("\n   Simulating BPSK.. Please wait until the computation terminate");
	int err=0;
#ifndef VARIANCE_STUDY	
	FILE *f=fopen("bpsktest.txt", "wt");
	fprintf(f, "\nBINARY_IN;BPSK;RAYLEIGH_RE;RAYLEIGH_IM;RAYLEIGH_MAG;BPSK_PERTUBE;BINARY_OUT");
#else
	FILE *f=fopen("tebbpsk.txt","wt");
	fprintf(f,"VARIANCE;TEB");
#endif
#ifdef VARIANCE_STUDY
 //for (EBN0=EBN0_MIN_DB;EBN0<EBN0_MAX_DB;EBN0+=EBN0_STEP) {
	for(VARIANCE=VARIANCE_MIN;VARIANCE<VARIANCE_MAX;VARIANCE+=VARIANCE_STEP) {
	double awgnseed=RAND_MAX/(rand()*3);
	double EBN0_LIN = pow(10,EBN0/10.0);
	//VARIANCE = SYMBOL_TIME*1.0e-6 /(2.0* EBN0_LIN);
	nspdRandGausInit(awgnseed, MOYENNE, VARIANCE, &dgausStatePtr);
	printf("\nComputing VARIANCE = %f...", VARIANCE);
#endif
	err=0;
	for(int w=0;w<100;w++) {
		DCplx **taps = fadingtaps(MAX_DOPPLER,SYMBOL_TIME,FILTER_LENGTH,NB_SAMPLES,NB_RAYS);
		
		short *binary_in = nspwMalloc(NB_SAMPLES);	
		nspwbRandUni(&uniStatePtr, binary_in, NB_SAMPLES);

		double *bpsk = nspdBPSK_mod(binary_in, NB_SAMPLES);
		//for(int u=0;u<NB_SAMPLES;u++) printf(" %d-->%.1f    ", binary_in[u], bpsk[u]);
		double *bpsk_pertube = nspdMalloc(NB_SAMPLES);
		nspdbCopy(bpsk, bpsk_pertube, NB_SAMPLES);
		
		//printf("\n");

		for(int t=0;t<NB_SAMPLES;t++){
			bpsk_pertube[t]=bpsk[t]*sqrt(taps[t][0].re*taps[t][0].re+taps[t][0].im*taps[t][0].im);
		//	printf("\t %f", bpsk_pertube[t]);
		}

		dbAWGN(bpsk_pertube,NB_SAMPLES);

		short *binary_out=nspdBPSK_demod(bpsk_pertube,NB_SAMPLES);
		
		//for(u=0;u<NB_SAMPLES;u++) printf(" \n%d-->%d    ", binary_in[u], binary_out[u]);

#ifndef VARIANCE_STUDY	
		for(t=0;t<NB_SAMPLES;t++)
			fprintf(f, "\n%d;%.3f;%.3f;%.3f;%.3f;%.3f;%d;",binary_in[t], bpsk[t], taps[t][0].re,taps[t][0].im,sqrt(taps[t][0].re*taps[t][0].re+taps[t][0].im*taps[t][0].im),bpsk_pertube[t], binary_out[t]);
#endif
		nspwbSub2(binary_in, binary_out, NB_SAMPLES, NSP_NO_SCALE,0);
		nspwbAbs1(binary_out, NB_SAMPLES);
		err+=nspwSum(binary_out, NB_SAMPLES, NSP_NO_SCALE,0);

		nspFree(bpsk_pertube);
		nspFree(bpsk);
		nspFree(binary_in);
		nspFree(binary_out);
		for(int g=0;g<NB_SAMPLES;g++)
			nspFree(taps[g]);
		free(taps);
	}
#ifdef VARIANCE_STUDY
	printf("..done...TEB = %f", (double)err/(double)(NB_ITERATION*NB_SAMPLES));
	fprintf(f,"\n%lf;%g", VARIANCE, (double)err/(double)(NB_ITERATION*NB_SAMPLES));
	}
#endif
#ifndef VARIANCE_STUDY
	fprintf(f,"TEB is %d", err);
	fclose(f);
#endif

#else // QPSK
	printf("\n-**************** ENTERING QPSK TEST SYSTEM ********************-\n\n");
	int err=0;
	FILE *f=fopen("qpskogo.txt", "wt");
	fprintf(f, "\nBINARY_IN;;QPSK.RE;QPSK.IM;RAYLEIGH_RE;RAYLEIGH_IM;RAYLEIGH_MAG;QPSK_PERTUBE.RE;QPSK_PERTUBE.IM;BINARY_OUT;;");
	for(int w=0;w<NB_ITERATION;w++) {
		short *binary_in = nspwMalloc(2*NB_SAMPLES);	
		nspwbRandUni(&uniStatePtr, binary_in, 2*NB_SAMPLES);
		int qpsk_size;
		DCplx *qpsk = GRAY_QPSK_mod(binary_in, 2*NB_SAMPLES,&qpsk_size);
		//for(int u=0;u<NB_SAMPLES;u++) printf(" %d-->%.1f    ", binary_in[u], bpsk[u]);
		DCplx *qpsk_pertube = nspzMalloc(qpsk_size);
		nspzbCopy(qpsk, qpsk_pertube, qpsk_size);
		DCplx **taps = fadingtaps(MAX_DOPPLER,SYMBOL_TIME,FILTER_LENGTH,qpsk_size,NB_RAYS);
		
		
		//printf("\n");
		double amp;
		DCplx compensation;
		for(int t=0;t<qpsk_size;t++){
			amp=sqrt(taps[t][0].re*taps[t][0].re+taps[t][0].im*taps[t][0].im);
			qpsk_pertube[t].re=amp;
			qpsk_pertube[t].im=amp;	
		}
		zbAWGN(qpsk_pertube,qpsk_size);
		short *binary_out=nspdBPSK_demod(bpsk_pertube,qpsk_size);
		
		//for(u=0;u<NB_SAMPLES;u++) printf(" \n%d-->%d    ", binary_in[u], binary_out[u]);

		
		for(t=0;t<qpsk_size;t++)
			fprintf(f, "\n%d;%d;%.3f;%.3f;%.3f;%.3f;%.3f;%.3f;%.3f;%d;%d",binary_in[2t], binary_in[2t+1], qpsk[t].re, qpsk[t].im, taps[t][0].re,taps[t][0].im,sqrt(taps[t][0].re*taps[t][0].re+taps[t][0].im*taps[t][0].im),qpsk_pertube[t].re,qpsk_pertube[t].im, binary_out[2t],binary_out[2t+1]);

		nspwbSub2(binary_in, binary_out, NB_SAMPLES, NSP_NO_SCALE,0);
		nspwbAbs1(binary_out, NB_SAMPLES);
		err+=nspwSum(binary_out, NB_SAMPLES, NSP_NO_SCALE,0);

		nspFree(bpsk_pertube);
		nspFree(bpsk);
		nspFree(binary_in);
		nspFree(binary_out);
		for(int g=0;g<NB_SAMPLES;g++)
			nspFree(taps[g]);
		free(taps);
	}
	fprintf(f,"Error is %d out of %d samples", err, NB_SAMPLES);
	fclose(f);

#endif // QPSK
	printf("\n\a\nComputation done...\nError is %d out of %d samples\nTEB = %lf\n", err, NB_SAMPLES*NB_ITERATION, (double)err/(double)(NB_ITERATION*NB_SAMPLES));

#endif //taptest

#endif

	return 0;
}

