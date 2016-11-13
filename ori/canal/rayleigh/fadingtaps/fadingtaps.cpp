// fadingtaps.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "globaldef.h"
#include <math.h>
#include <stdlib.h>


#ifdef _NSP_INTEL_
#define nsp_UsesVector
#define nsp_UsesConversion
#define nsp_UsesSampleGen
#define nsp_UsesConvolution

#include <nsp.h>

#endif
#include "fadingtdma.h"
#include "fadingcdma.h"
#include "psk.h"
#include "awgn.h"
#include "gold/gold.h"
#include "interpol/interpol.h"
#include "rayleigh.h"


#define FILTER_LENGTH	3
#define NB_SAMPLES		50
#define MAX_DOPPLER		100
#define NB_RAYS			500
#define SYMBOL_TIME		2


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

NSPDRandGausState  rayleighDGausStatePtr;

int main(int argc, char* argv[]){
	int N = 0;
	nspdRandGausInit(37, 0.0, 1.0, &rayleighDGausStatePtr);
	double *r = RayleighFading(95,1050 ,512, &N);
	FILE *f = fopen("goloki.txt", "wt");
	for(int i = 0; i<N;i++) 
		fprintf(f, "\n%d : %.6f", i, r[i]);
	fclose(f);

	return 0;
}

#ifdef bluurp

int main(int argc, char* argv[]){

	/*short delay[]    = {  0,   1,    4,    5,    8,   10,   11,   15};
	double raycoef[] = {0.8, 0.4, 0.25, 0.24, 0.19, 0.15, 0.14, 0.07};
	*/
	int InitialValue[2]={31,31};
	int Polynomial[2]={18, 27};

	short delay[]={  0, 1, 2, 3};
	double raycoef[]= { 0.8 , -0.5, 0.3, -0.2};
	short *wdata= nspwMalloc(NB_SAMPLES);

	double *bpskmultipath;
	
	nspwRandUniInit (1010100, 0, 2, &uniStatePtr);
	nspwbRandUni(&uniStatePtr, wdata, NB_SAMPLES);


	double data[NB_SAMPLES];
	nspdbIntToFloat(wdata, data, NB_SAMPLES, sizeof(short)*8, NSP_Noflags);
	nspFree(wdata);
	
	int code_size=0, spreadsize=0;

	double *gold=dGoldCode(InitialValue, Polynomial, 5, 5 , 0, &code_size);
	double *spreaded=dSpreadData(data, NB_SAMPLES, gold, code_size, &spreadsize);
	double *bpsk = nspdBPSK_mod(spreaded, spreadsize);

	TDRayleighMultipath ray = dInitRayleighMultipath(delay, raycoef, 4);
	bpskmultipath = dRayleighMultipath_nip(ray, bpsk, spreadsize);

	double *bpsk_demod = nspdBPSK_demod(bpskmultipath, spreadsize);

	int unspr_size=0;
	double *unspr = dUnspreadData(bpsk_demod, spreadsize, gold, code_size, &unspr_size);
	for(int k = 0; k< unspr_size;k++)
		unspr[k]=unspr[k]>0.0?1.0:0.0;

	printf("\n filter size  :%d", ray.taplen);
	for(int i = 0;i< NB_SAMPLES;i++) {
		printf("\n i= %d     data = %lf  res = %lf   status : %d", i, data[i], unspr[i], data[i]-unspr[i]);
	}

#ifdef louikjs
	double VARIANCE = 0.0;
	nspzRandGausInit(0, MOYENNE, VARIANCE, &zgausStatePtr);
	nspdRandGausInit(0, MOYENNE, VARIANCE, &dgausStatePtr);
	nspdRandUniInit(0, 0, NSP_2PI, &tapstatePtr);
	nspwRandUniInit (0, 0, 2, &uniStatePtr);
	DCplx **taps = fadingtaps(MAX_DOPPLER,SYMBOL_TIME,FILTER_LENGTH,NB_SAMPLES,NB_RAYS);

	for (int i = 0; i < 10; i++) 
		printf("\n(%.3lf, %.3lf) (%.3lf, %.3lf) (%.3lf, %.3lf)", taps[i][0].re,  taps[i][0].im, taps[i][1].re,  taps[i][1].im, taps[i][2].re,  taps[i][2].im); 
#endif

	return 0;
}

#endif
#ifdef HUTUIL

int bpsktest(int argc, char* argv[])
{

#ifdef _NSP_INTEL_
	double EBN0=0.0001;
	double VARIANCE=0.0;
	double VARIANCE_STEP=1.0;
	printf("\nEnter AWGN VARIANCE Value: ");
	scanf("%lf", &VARIANCE_STEP );

	nspzRandGausInit(314159, MOYENNE, VARIANCE, &zgausStatePtr);
	nspdRandGausInit(445755, MOYENNE, VARIANCE, &dgausStatePtr);
	nspdRandUniInit(1560, 0, NSP_2PI, &tapstatePtr);
	nspwRandUniInit (114521, 0, 2, &uniStatePtr);
#ifndef QPSK
	printf("\n-**************** ENTERING BPSK TEST SYSTEM ********************-\n\n");
	printf("\n\a*** Multipath Rayleigh Fading Channel Parameter ***\n\nMAX_DOPPLER : %d\nSYMBOL_TIME : %f\nFILTER_LENGTH : %d\nNB_SAMPLES : %d\nNB_RAYS : %d", MAX_DOPPLER,SYMBOL_TIME,FILTER_LENGTH,NB_SAMPLES,NB_RAYS);
	printf("\n\n\n*** Gaussian Random Parameter ***\n\nMOYENNE : %f\nCalculating TEB for VARIANCE varying from %f to %f with step %f", MOYENNE, VARIANCE_MIN, VARIANCE_MAX, VARIANCE_STEP);
	printf("\n\n*****************************************\n");
	printf("\n   Simulating BPSK.. Please wait until the computation terminate");
	int err=0;
	printf("\nComputing VARIANCE = %f...", VARIANCE);
	err=0;
	int InitialValue[2]={31,31};
	int Polynomial[2]={18, 27};
	for(int w=0;w<100;w++) {
		DCplx **taps = fadingtaps(MAX_DOPPLER,SYMBOL_TIME,FILTER_LENGTH,NB_SAMPLES,NB_RAYS);
		
		short *binary_in = nspwMalloc(NB_SAMPLES);	
		nspwbRandUni(&uniStatePtr, binary_in, NB_SAMPLES);
		int code_size=0, spreadsize=0;
		double *gold=dGoldCode(InitialValue, Polynomial, 5, 5 , 0, &code_size);
		double *spreaded=dSpreadData(binary_in, NB_SAMPLES, gold, code_size, &spreadsize);
		double *bpsk = nspdBPSK_mod(spreaded,spreadsize);
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
	printf("..done...TEB = %f", (double)err/(double)(NB_ITERATION*NB_SAMPLES));
	fprintf(f,"\n%lf;%g", VARIANCE, (double)err/(double)(NB_ITERATION*NB_SAMPLES));

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

#endif
	return 0;
}

#endif 