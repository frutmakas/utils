/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/08/03 13:16:57 $
 * $Revision: 1.1.2.11 $
 * $Id: psk.cpp,v 1.1.2.11 2006/08/03 13:16:57 syed Exp $
 ********************************************************************/
 
// qpskaw gn.cpp : Defines the entry point for the console application.
//

/*#include  <math.h>
#include  <stdio.h>*/
#include <iostream>
using namespace std;

#include "modulation/psk.h"
#include "tools/tools.h"
#ifdef WIN_DOS
#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 
#endif

TPSK8 psk8;
TPSK8_GRAY psk8gray;
TPSK4 psk4;
TPSK4_GRAY psk4gray;
TBPSK bpsk;

ZVector z8PSK_mod(const WVector &data){
	if (!psk8.initialized) {
		fprintf(stderr , "INITIALIZE PSK8 FIRST !!! call psk8init()");
		return ZVector();
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
		return WVector();
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



// gray
ZVector z8PSKGray_mod(const WVector &data){
	if (!psk8gray.initialized) {
		fprintf(stderr , "INITIALIZE psk8gray FIRST !!! call psk8grayinit()");
		return ZVector();
	}
	ZVector data_out((int)ceil(data.taille/3.0));
//	*psk8gray_size = (int)ceil(data.taille/3.0);
//	DCplx *data_out=nspzMalloc(*psk8gray_size);
	for(register int i = 0, j=0; i < data.taille; i+=3, j++) {
#ifndef DEBUG
		data_out.vect[j]=psk8gray.sym[bin2oct(&data.vect[i], data.taille-i)];
#else	
		int oct = bin2oct(&data.vect[i], data.taille-i);
		data_out.vect[j]=psk8gray.sym[oct];
#endif
	}
	return data_out;
}

WVector z8PSKGray_demod(const ZVector &data){
	if (!psk8gray.initialized) {
		fprintf(stderr , "INITIALIZE psk8gray FIRST !!! call psk8grayinit()");
		return WVector();
	}

	WVector data_out(data.taille*3);

//	*binary_size=data.taille*3;
//	short *data_out=nspwMalloc(*binary_size);

	for(int i  =0 ,j=0; i < data.taille; i++, j+=3) {
		int min_pos=0;
		double min_distance = (data.vect[i]-psk8gray.sym[0]).mag();
		double dist = 0.0;
		for(int k = 1; k < 8 ;k++) {
			dist = (data.vect[i]-psk8gray.sym[k]).mag();
			if (dist < min_distance) {
				min_distance=dist; min_pos=k;
			}
		}
		oct2bin(min_pos, &(data_out.vect[j]));
	}

	return data_out;
}

/// PSK4 Gray
/*  Modulation QPSK shifted, codage Gray  */
ZVector zPSK4ShiftedGray_mod(const WVector &data){
	long i;
	double amp;
	ZVector x(data.taille>>1);
	amp=sqrt(2.0)/2.0;
	for (i=0; i<x.taille; i++){
		x.vect[i].re=(data.vect[2*i]==0) ? amp : -amp;
		x.vect[i].im=(data.vect[2*i+1]==0) ? amp : -amp;
	}
	return x;
}
/*-------------------------------------------------------*/
/*  Demodulation QPSK shifted, codage Gray  */
WVector zPSK4ShiftedGray_demod(const ZVector &x)
{
	WVector data_hat(x.taille*2);
	long i;
	for (i=0; i<x.taille; i++){
		data_hat.vect[2*i] = (x.vect[i].re >= 0) ? 0 : 1;
		data_hat.vect[2*i+1] = (x.vect[i].im >= 0) ? 0 : 1;
	}
	return data_hat;
}


ZVector zPSK4Gray_mod(const WVector &data){
	ZVector data_out(data.taille>>1);
	for(register int i = 0, j=0; i < data.taille; i+=2, j++) {
		int quad = bin2quad(&data.vect[i], data.taille-i);
		data_out.vect[j]=psk4gray.sym[quad];
	}
	return data_out;
}

/*-------------------------------------------------------*/
/*  Demodulation QPSK shifted, codage Gray  */
WVector zPSK4Gray_demod(const ZVector &data) {
	WVector data_out(data.taille*2);
	for(int i  =0 ,j=0; i < data.taille; i++, j+=2) {
		int min_pos=0;
		double min_distance = (data.vect[i]-psk4gray.sym[0]).mag();
		double dist = 0.0;
		for(int k = 1; k < 4 ;k++) {
			dist = (data.vect[i]-psk4gray.sym[k]).mag();
			if (dist < min_distance) {
				min_distance=dist; min_pos=k;
			}
		}
		quad2bin(min_pos, &(data_out.vect[j]));
	}
	return data_out;
}

ZVector z4PSK_mod(const WVector &data){
	ZVector data_out(data.taille>>1);
	for(register int i = 0, j=0; i < data.taille; i+=2, j++) {
		int quad = bin2quad(&data.vect[i], data.taille-i);
		data_out.vect[j]=psk4.sym[quad];
	}
	return data_out;
}

WVector z4PSK_demod(const ZVector &data){
	WVector data_out(data.taille*2);
	for(int i  =0 ,j=0; i < data.taille; i++, j+=2) {
		int min_pos=0;
		double min_distance = (data.vect[i]-psk4.sym[0]).mag();
		double dist = 0.0;
		for(int k = 1; k < 4 ;k++) {
			dist = (data.vect[i]-psk4.sym[k]).mag();
			if (dist < min_distance) {
				min_distance=dist; min_pos=k;
			}
		}
		quad2bin(min_pos, &(data_out.vect[j]));
	}
	return data_out;
}

ZVector z4PSKGray_mod(const WVector &data){
	ZVector data_out(data.taille>>1);
	for(register int i = 0, j=0; i < data.taille; i+=2, j++) {
		int quad = bin2quad(&data.vect[i], data.taille-i);
		data_out.vect[j]=psk4gray.sym[quad];
	}
	return data_out;
}

WVector z4PSKGray_demod(const ZVector &data){
	WVector data_out(data.taille*2);
	for(int i  =0 ,j=0; i < data.taille; i++, j+=2) {
		int min_pos=0;
		double min_distance = (data.vect[i]-psk4gray.sym[0]).mag();
		double dist = 0.0;
		for(int k = 1; k < 4 ;k++) {
			dist = (data.vect[i]-psk4gray.sym[k]).mag();
			if (dist < min_distance) {
				min_distance=dist; min_pos=k;
			}
		}
		quad2bin(min_pos, &(data_out.vect[j]));
	}
	return data_out;
}

DVector dBPSK_mod(const WVector &data) {
    DVector data_out(data.taille);
    for(register int i=0;i<data.taille;i++) {
        data_out.vect[i]=data.vect[i] ? 1 : -1;
    }    
    return data_out;
}    

ZVector zBPSK_mod(const WVector &data) {
    ZVector data_out(data.taille);
    for(register int i=0;i<data.taille;i++) {
        data_out.vect[i]=data.vect[i] ? 1 : -1;
    }    
    return data_out;
}    

WVector zBPSK_demod(const ZVector &data) {
    WVector data_out(data.taille);
    for(register int i=0;i<data.taille;i++) {
        data_out.vect[i]=data.vect[i].re > 0 ? 1 : 0;
    }    
    return data_out;
}    

WVector dBPSK_demod(const DVector &data) {
    WVector data_out(data.taille);
    for(register int i=0;i<data.taille;i++) {
        data_out.vect[i]=data.vect[i] > 0 ? 1 : 0;
    }    
    return data_out;
}    

/*void psk8init(){
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
}*/
