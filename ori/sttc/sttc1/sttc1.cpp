#include "globaldef.h"
#include <conio.h>

#ifdef _NSP_INTEL_

#define nsp_UsesVector
#define nsp_UsesConversion
#define nsp_UsesSampleGen

#include <nsp.h>
#include "conversion.h"
#include "psk.h"
#include <stdio.h>
#include "tools.h"
#include "awgn.h"


#include "nspmatrice.h"
#include "nspdcplx.h"



NSPDRandGausState  rayleighDGausStatePtr;
NSPZRandGausState  zgausStatePtr;
NSPZRandGausState  zgausnoiseStatePtr;
NSPDRandGausState  dgausStatePtr;
NSPWRandUniState  wuniStatePtr;

#define FILTER_LENGTH	1
//#define NB_SAMPLES		50
#define MAX_DOPPLER		1500
#define NB_RAYS			500
#define SYMBOL_TIME		2

const NB_BINARY=48; 

TPSK8 psk8;

int psk8testmain(int argc, char **argv) {
	psk8.initialized=false;
	short *binary_in=nspwMalloc(NB_BINARY);
	nspwRandUniInit(75152, 0, 2, &wuniStatePtr);
	nspwbRandUni(&wuniStatePtr, binary_in, NB_BINARY);
	nspzRandGausInit(17552, 0.0, 0.05, &zgausStatePtr);
	int psk8_size;
	psk8init(); 
	
	DCplx *psk8signal = nspz8PSK_mod(binary_in, NB_BINARY, &psk8_size);
	zbAWGN(psk8signal, psk8_size);
	int out_size=0;
	short *binary_out = nspz8PSK_demod(psk8signal, psk8_size, &out_size);
	int sum = 0;
	for(int i=0;i<out_size;i+=3) {
		printf("\n %d --> %d", bin2oct(&(binary_in[i]), NB_BINARY-i), bin2oct(&(binary_out[i]), NB_BINARY-i));
	}
	for( i=0;i<NB_BINARY;i++) 
		sum += (binary_in[i]==binary_out[i]) ? 0 : 1;

	
	printf("\n SUM ===> %d, NB_BINARY = %d outsize=%d", sum, NB_BINARY, out_size);
/*	for(int j=0,i=0;j<psk8_size;j++,i+=3) {
		printf("\n %d --> ", bin2oct(&(binary_in[i]), NB_BINARY-i));
		print(psk8signal[j]);
	}
	*/
	return 1;
}

int main(int argc, char **argv) {
	psk8init();
	short *binary_in=nspwMalloc(NB_BINARY);
	nspzRandGausInit(17552, 0.0, 0.5, &zgausStatePtr);

	nspdRandGausInit(71552, 0.0, 0.5, &dgausStatePtr);
	nspwRandUniInit(75152, 0, 2, &wuniStatePtr);
//	nspdRandUniInit(75512, 0, 1, &tapstatePtr);

	nspwbRandUni(&wuniStatePtr, binary_in, NB_BINARY);
	
	int NB_DATA = 0;
	DCplx *data = nspz8PSK_mod(binary_in, NB_BINARY, &NB_DATA);

	nspzRandGausInit(17552, 0.0, 0.1, &zgausnoiseStatePtr); // 0.712 a changer en 3/2SNR
	//nspzRandGausInit(17552, 0.0, 0.03 , &zgausnoiseStatePtr); // 0.712 a changer en 3/2SNR
	
	ZMatrix c_it_g2=G2(data, NB_DATA);

	
	// generating fading 
//	DCplx **taps = NULL;
	const NB_RECEIVE_ANTENNAS=1;
	const NB_TRANSMIT_ANTENNAS=2;

	ZMatrix alpha(NB_TRANSMIT_ANTENNAS, NB_RECEIVE_ANTENNAS);
	for(int k=0;k<NB_TRANSMIT_ANTENNAS;k++) {
		nspzbRandGaus(&zgausStatePtr ,alpha.mat[k], NB_RECEIVE_ANTENNAS);
//		alpha.mat[k][0].re=1.0;alpha.mat[k][0].im=0.0;
	}

	print("alpha ...", alpha);
	print("c_it_g2", c_it_g2);
	print("tmp", c_it_g2*alpha);

	ZMatrix sum=c_it_g2*alpha;
	ZMatrix noise(NB_RECEIVE_ANTENNAS, sum.line);
	for(k=0;k<NB_RECEIVE_ANTENNAS;k++) 
		nspzbRandGaus(&zgausnoiseStatePtr, noise.mat[k], sum.line);
		//nspzbZero(noise.mat[k], sum.line);
	noise.transpose();

	print("noise ...", noise);
	ZMatrix rtj=c_it_g2*alpha+noise;
	print("rtj", rtj);

	// decoding .... 

	DCplx *decode = nspzMalloc(rtj.line);
	for(int t = 0; t<rtj.line;t+=2) {
		DCplx x1={0.0, 0.0};
		DCplx x2={0.0, 0.0};
		double d1;
		double y=0.0;
		for (int j=0;j<NB_RECEIVE_ANTENNAS;j++) {
			x1=x1+rtj.mat[t][j]*conj(alpha.mat[0][j])+conj(rtj.mat[t+1][j])*alpha.mat[1][j];
			x2=x2+rtj.mat[t][j]*conj(alpha.mat[1][j])-conj(rtj.mat[t+1][j])*alpha.mat[0][j];
			for(int i=0;i<NB_TRANSMIT_ANTENNAS;i++) {
				y+=dcplxmod2(alpha.mat[i][j]);
			}
		}
		y=y-1;
		double min_metric1=dcplxmod2(x1-psk8.sym[0])+y*dcplxmod2(psk8.sym[0]);
		int min_idx1=0;
		double min_metric2=dcplxmod2(x2-psk8.sym[0])+y*dcplxmod2(psk8.sym[0]);
		int min_idx2=0;
		double tmp_metric=0.0;
		for(int k=1;k<8;k++) {
			tmp_metric=dcplxmod2(x1-psk8.sym[k])+y*dcplxmod2(psk8.sym[k]);
			if(tmp_metric<min_metric1) {
				min_idx1=k; min_metric1=tmp_metric;
			}
			tmp_metric=dcplxmod2(x2-psk8.sym[k])+y*dcplxmod2(psk8.sym[k]);
			if(tmp_metric<min_metric2) {
				min_idx2=k; min_metric2=tmp_metric;
			}
		}
		decode[t]=psk8.sym[min_idx1];
		decode[t+1]=psk8.sym[min_idx2];
	}
	int out_size=0;
	short *binary_out = nspz8PSK_demod(decode, rtj.line, &out_size);
	int sum1 = 0;
	for(int h=0;h<out_size;h++) {
		sum1 += (binary_in[h]==binary_out[h]) ? 0 : 1;
		printf("\n in : %d ----> out : %d", binary_in[h], binary_out[h]);
	}
	printf("\n Sum : %d", sum1);

	

	// print("rtj ...", rtj);




	return 1;
}

#else

#include "psk.h"
#include "utilitis.h"
#include "randgen.h"
#include "awgn.h"
#include "tools.h"
#include "conversion.h"

// NSPWRandUniState  wuniStatePtr;

#define FILTER_LENGTH	1
//#define NB_SAMPLES		50
#define MAX_DOPPLER		1500
#define NB_RAYS			500
#define SYMBOL_TIME		2

const NB_BINARY=10000; 

TPSK8 psk8;

dRandGausStatePtr  rayleighDGausStatePtr;
zRandGausStatePtr  zgausStatePtr;
zRandGausStatePtr  zgausnoiseStatePtr;
dRandGausStatePtr  dgausStatePtr;
TRandBitStatePtr wuniStatePtr;

int psk8testmain(int argc, char **argv) {

	psk8.initialized=false;

	wuniStatePtr.iseed = 5543;

	WVector binary_in(NB_BINARY);
	bRandBit1(&wuniStatePtr, binary_in);
	zRandGausInit(&zgausStatePtr, 0, 0.05);
//	int psk8_size;
	psk8init(); 
	
	ZVector psk8signal = z8PSK_mod(binary_in);
	zbAWGN(psk8signal);
	int out_size=0;
	WVector binary_out = z8PSK_demod(psk8signal);
	int sum = 0;
	for(int i=0;i<out_size;i+=3) {
		printf("\n %d --> %d", bin2oct(&(binary_in.vect[i]), NB_BINARY-i), bin2oct(&(binary_out.vect[i]), NB_BINARY-i));
	}
	for( i=0;i<NB_BINARY;i++) 
		sum += (binary_in.vect[i]==binary_out.vect[i]) ? 0 : 1;

	
	printf("\n SUM ===> %d, NB_BINARY = %d outsize=%d", sum, NB_BINARY, out_size);
/*	for(int j=0,i=0;j<psk8_size;j++,i+=3) {
		printf("\n %d --> ", bin2oct(&(binary_in[i]), NB_BINARY-i));
		print(psk8signal[j]);
	}
	*/
	return 1;
}
 
int main(int argc, char **argv) {
	psk8init();
	WVector binary_in(NB_BINARY);

	zRandGausInit(&zgausStatePtr, 0.5, 0.5);

	dRandGausInit(&dgausStatePtr, 0, 0.5);
	wuniStatePtr.iseed=6656;
//	nspdRandUniInit(75512, 0, 1, &tapstatePtr);

	bRandBit1(&wuniStatePtr, binary_in);
	cout << endl <<"input binaries : " << binary_in;
	
	int NB_DATA = 0;
	ZVector data = z8PSK_mod(binary_in);
	cout << endl <<"PSK8 modulated symbols : " ; // << data;

	zRandGausInit(&zgausnoiseStatePtr, 0.0,0.5); // 0.712 a changer en 3/2SNR
	//nspzRandGausInit(17552, 0.0, 0.03 , &zgausnoiseStatePtr); // 0.712 a changer en 3/2SNR
	
	ZMatrix c_it_g2=G2(data, data.taille);
	cout << endl <<"space time G2 coded output: "; // <<  c_it_g2;

	// generating fading 
//	DCplx **taps = NULL;
	const NB_RECEIVE_ANTENNAS=2;
	const NB_TRANSMIT_ANTENNAS=2;

	ZMatrix alpha(NB_TRANSMIT_ANTENNAS, NB_RECEIVE_ANTENNAS);
	for(int k=0;k<NB_TRANSMIT_ANTENNAS;k++) {
		zbRandGaus(&zgausStatePtr ,alpha.mat[k], NB_RECEIVE_ANTENNAS);
	}

	cout << endl << "alpha ..." << alpha.mag();
	getch();

	cout << endl << "ij path fading : "; // << alpha;

	ZMatrix sum=c_it_g2*alpha;
	ZMatrix noise(NB_RECEIVE_ANTENNAS, sum.line);
	for(k=0;k<NB_RECEIVE_ANTENNAS;k++) 
		zbRandGaus(&zgausnoiseStatePtr, noise.mat[k], sum.line);
		//nspzbZero(noise.mat[k], sum.line);
	noise.transpose();
	cout << endl <<"ij path noise over time: "; // << binary_in;
	
	//print("noise ...", noise);
	ZMatrix rtj=sum+noise;
	cout << endl <<"received signal : "; //<< rtj;

	//print("rtj", rtj);

	// decoding .... 

	ZVector decode(rtj.line);
	for(int t = 0; t<rtj.line;t+=2) {
		DCplx x1(0.0, 0.0);
		DCplx x2(0.0, 0.0);
		double d1;
		double y=0.0;
		for (int j=0;j<NB_RECEIVE_ANTENNAS;j++) {
			x1=x1+rtj.mat[t][j]*alpha.mat[0][j].conj()+rtj.mat[t+1][j].conj()*alpha.mat[1][j];
			x2=x2+rtj.mat[t][j]*alpha.mat[1][j].conj()-rtj.mat[t+1][j].conj()*alpha.mat[0][j];
			for(int i=0;i<NB_TRANSMIT_ANTENNAS;i++) {
				y+=alpha.mat[i][j].mod2();
			}
		}
		y=y-1;
		double min_metric1=(x1-psk8.sym[0]).mod2()+y*psk8.sym[0].mod2();
		int min_idx1=0;
		double min_metric2=(x2-psk8.sym[0]).mod2()+y*psk8.sym[0].mod2();
		int min_idx2=0;
		double tmp_metric=0.0;
		for(int k=1;k<8;k++) {
			tmp_metric=(x1-psk8.sym[k]).mod2()+y*psk8.sym[k].mod2();
			if(tmp_metric<min_metric1) {
				min_idx1=k; min_metric1=tmp_metric;
			}
			tmp_metric=(x2-psk8.sym[k]).mod2()+y*psk8.sym[k].mod2();
			if(tmp_metric<min_metric2) {
				min_idx2=k; min_metric2=tmp_metric;
			}
		}
		decode.vect[t]=psk8.sym[min_idx1];
		decode.vect[t+1]=psk8.sym[min_idx2];
	}
	sum.~ZMatrix();
	noise.~ZMatrix();
	
	int out_size=0;
	WVector binary_out = z8PSK_demod(decode).copy(0, binary_in.taille);
	rtj.~ZMatrix();
	decode.~ZVector();
	WVector diff=binary_out-binary_in;
	cout << endl <<"demodulated and decoded signal : " << binary_out;
	cout << endl <<"error : " << diff;

	
	int sum1 = 0;
	for(int h=0;h<binary_in.taille;h++) {
		sum1 += (binary_in.vect[h]==binary_out.vect[h]) ? 0 : 1;
		if (binary_in.vect[h]!=binary_out.vect[h]) printf("\n @%d in : %d ----> out : %d", h, binary_in.vect[h], binary_out.vect[h]);
	}
	
	printf("\n Sum : %d", sum1);

	

	// print("rtj ...", rtj);




	return 1;
}

#endif
