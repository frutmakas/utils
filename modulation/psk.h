/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/09/07 10:07:35 $
 * $Revision: 1.1.2.14 $
 * $Id: psk.h,v 1.1.2.14 2006/09/07 10:07:35 syed Exp $
 ********************************************************************/
#ifndef _PSK_H_
#define _PSK_H_

#include "globaldef.h"

#include "tools/utilitis.h"
#include <iostream>
using namespace std;

void psk8init();

class TConstellation {
protected :
    WMatrix bitmap;
public:
	DCplx *sym;
	bool initialized;
	int size, bit_per_sym;
	string id;
	ZVector (*modulate_function)(const WVector &data);
	WVector (*demodulate_function)(const ZVector &data);
    WVector getBinary(const int idx) {
        if(idx > size || idx < 0) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);

        WVector bin(bit_per_sym);
        for(int i =0;i<bit_per_sym;i++) {
            bin.vect[i]=bitmap.mat[idx][i];
        }
        return bin;
    }
    DCplx decode(const DCplx &c) {
        double min=1e50, m;
        int idx = 0;
        for(int x = 0; x < size; x++) {
            if((m=(c-sym[x]).mag())<min) {
                min = m;
                idx = x;
            }
        }
        return sym[idx];
    }
	TConstellation() { initialized=false; sym=NULL; size=0; modulate_function=NULL; demodulate_function=NULL;}
	//~TConstellation();
};

ZVector z8PSK_mod(const WVector &data);
WVector z8PSK_demod(const ZVector &data);
//WVector z8PSK_getBinary(const int idx);

ZVector z8PSKGray_mod(const WVector &data);
WVector z8PSKGray_demod(const ZVector &data);
//WVector z8PSKGray_getBinary(const int idx);

ZVector z4PSKGray_mod(const WVector &data);
WVector z4PSKGray_demod(const ZVector &data);
//WVector z4PSKGray_getBinary(const int idx);

ZVector z4PSK_mod(const WVector &data);
WVector z4PSK_demod(const ZVector &data);
//WVector z4PSK_getBinary(const int idx);

DVector dBPSK_mod(const WVector &data);
WVector dBPSK_demod(const DVector &data);

ZVector zBPSK_mod(const WVector &data);
WVector zBPSK_demod(const ZVector &data);

class TPSK8 : public TConstellation {
public:
	//DCplx sym[8];
	bool initialized;
	TPSK8(){
		sym = new DCplx[8];
		size=8;
		bit_per_sym=3;
		initialized = true;
		modulate_function = &z8PSK_mod;
		demodulate_function = &z8PSK_demod;
        bitmap = WMatrix(size, bit_per_sym);
        int i=0,j=0;
        bitmap.mat[j][i++] = 0; bitmap.mat[j][i++] = 0; bitmap.mat[j++][i++] = 0; i=0;
        bitmap.mat[j][i++] = 0; bitmap.mat[j][i++] = 0; bitmap.mat[j++][i++] = 1; i=0;
        bitmap.mat[j][i++] = 0; bitmap.mat[j][i++] = 1; bitmap.mat[j++][i++] = 0; i=0;
        bitmap.mat[j][i++] = 0; bitmap.mat[j][i++] = 1; bitmap.mat[j++][i++] = 1; i=0;
        bitmap.mat[j][i++] = 1; bitmap.mat[j][i++] = 0; bitmap.mat[j++][i++] = 0; i=0;
        bitmap.mat[j][i++] = 1; bitmap.mat[j][i++] = 0; bitmap.mat[j++][i++] = 1; i=0;
        bitmap.mat[j][i++] = 1; bitmap.mat[j][i++] = 1; bitmap.mat[j++][i++] = 0; i=0;
        bitmap.mat[j][i++] = 1; bitmap.mat[j][i++] = 1; bitmap.mat[j++][i++] = 1; i=0;

		id = "8-PSK";

		sym[0].re = 1.0;
		sym[0].im = 0.0;
		sym[1].re = 0.70711;
		sym[1].im = 0.70711;
		sym[2].re = 0.0;
		sym[2].im = 1.0;
		sym[3].re = -0.70711;
		sym[3].im = 0.70711;

		sym[4].re = -1.0;
		sym[4].im = 0.0;
		sym[5].re = -0.70711;
		sym[5].im = -0.70711;
		sym[6].re = 0.0;
		sym[6].im = -1.0;
		sym[7].re = 0.70711;
		sym[7].im = -0.70711;
	}
	~TPSK8() {
		if(sym) delete[] sym;
	}
};

extern TPSK8 psk8;

class TPSK8_GRAY : public TConstellation {
public:
	//DCplx sym[8];
	bool initialized;
	TPSK8_GRAY(){
		sym = new DCplx[8];
		size=8;
		bit_per_sym=3;
		initialized = true;
		id = "8-PSK (Gray)";

		modulate_function = &z8PSKGray_mod;
		demodulate_function = &z8PSKGray_demod;
        bitmap = WMatrix(size, bit_per_sym);
		sym[0].re = 1.0; // 000
		sym[0].im = 0.0;
		sym[1].re = 0.70711; // 001
		sym[1].im = 0.70711;
		sym[3].re = 0.0; // 011
		sym[3].im = 1.0;
		sym[2].re = -0.70711; //010
		sym[2].im = 0.70711;

		sym[6].re = -1.0; //110
		sym[6].im = 0.0;
		sym[7].re = -0.70711;//111
		sym[7].im = -0.70711;
		sym[5].re = 0.0;//101
		sym[5].im = -1.0;
		sym[4].re = 0.70711; //100
		sym[4].im = -0.70711;

        int i=0,j=0;
        bitmap.mat[j][i++] = 0; bitmap.mat[j][i++] = 0; bitmap.mat[j++][i++] = 0; i=0;
        bitmap.mat[j][i++] = 0; bitmap.mat[j][i++] = 0; bitmap.mat[j++][i++] = 1; i=0;
        bitmap.mat[j][i++] = 0; bitmap.mat[j][i++] = 1; bitmap.mat[j++][i++] = 0; i=0;
        bitmap.mat[j][i++] = 0; bitmap.mat[j][i++] = 1; bitmap.mat[j++][i++] = 1; i=0;
        bitmap.mat[j][i++] = 1; bitmap.mat[j][i++] = 0; bitmap.mat[j++][i++] = 0; i=0;
        bitmap.mat[j][i++] = 1; bitmap.mat[j][i++] = 0; bitmap.mat[j++][i++] = 1; i=0;
        bitmap.mat[j][i++] = 1; bitmap.mat[j][i++] = 1; bitmap.mat[j++][i++] = 0; i=0;
        bitmap.mat[j][i++] = 1; bitmap.mat[j][i++] = 1; bitmap.mat[j++][i++] = 1; i=0;

	}

	~TPSK8_GRAY() {
		if(sym) delete[] sym;
	}

};

extern TPSK8_GRAY psk8gray;

class TPSK4_GRAY : public TConstellation  {
public:
	//DCplx sym[4];
	bool initialized;
	TPSK4_GRAY(){
		initialized = true;
		size=4;
		bit_per_sym=2;
		id = "QPSK (Gray)";
		modulate_function = &z4PSKGray_mod;
		demodulate_function = &z4PSKGray_demod;
        bitmap = WMatrix(size, bit_per_sym);

		sym = new DCplx[4];
		double val = 1.0/sqrt(2.0);
		sym[0].re = val;
		sym[0].im = val;
		sym[1].re = -val;
		sym[1].im = val;

		sym[3].re = -val;
		sym[3].im = -val;
		sym[2].re = val;
		sym[2].im = -val;
        bitmap = WMatrix(size, bit_per_sym);
        int i=0,j=0;
        bitmap.mat[j][i++] = 0; bitmap.mat[j++][i++] = 0; i=0;
        bitmap.mat[j][i++] = 0; bitmap.mat[j++][i++] = 1; i=0;
        bitmap.mat[j][i++] = 1; bitmap.mat[j++][i++] = 0; i=0;
        bitmap.mat[j][i++] = 1; bitmap.mat[j++][i++] = 1; i=0;

	}

	~TPSK4_GRAY() {
		if(sym) delete[] sym;
	}

};

extern TPSK4_GRAY psk4gray;


class TPSK4 : public TConstellation  {
public:
	//DCplx sym[4];
	bool initialized;
	TPSK4(){
		sym = new DCplx[4];
		size=4;
		bit_per_sym=2;
		id = "QPSK";
		modulate_function = &z4PSK_mod;
		demodulate_function = &z4PSK_demod;
        bitmap = WMatrix(size, bit_per_sym);

		initialized = true;
		double val = 1.0/sqrt(2.0);
		sym[0].re = val;
		sym[0].im = val;
		sym[1].re = -val;
		sym[1].im = val;

		sym[2].re = -val;
		sym[2].im = -val;
		sym[3].re = val;
		sym[3].im = -val;
        bitmap = WMatrix(size, bit_per_sym);
        int i=0,j=0;
        bitmap.mat[j][i++] = 0; bitmap.mat[j++][i++] = 0; i=0;
        bitmap.mat[j][i++] = 0; bitmap.mat[j++][i++] = 1; i=0;
        bitmap.mat[j][i++] = 1; bitmap.mat[j++][i++] = 0; i=0;
        bitmap.mat[j][i++] = 1; bitmap.mat[j++][i++] = 1; i=0;
	}
	~TPSK4() {
		if(sym) delete[] sym;
	}
};

extern TPSK4 psk4;


class TBPSK : public TConstellation  {
public:
	//DCplx sym[4];
	bool initialized;
	TBPSK(){
		sym = new DCplx[2];
		size=2;
		bit_per_sym=1;
		id = "BPSK [-1/+1]";
		modulate_function = &zBPSK_mod;
		demodulate_function = &zBPSK_demod;
		initialized = true;
		double val = 1.0;
		sym[0].re = val;
		sym[0].im = 0;
		sym[1].re = -val;
		sym[1].im = 0;
	}
	~TBPSK() {
		if(sym) delete[] sym;
	}
};

extern TBPSK bpsk;

/*
class TBPSK_z : public TConstellation  {
public:
	//DCplx sym[4];
	bool initialized;
	TBPSK_d(){
		sym = new DCplx[2];
		size=2;
		bit_per_sym=1;
		modulate_function = &zBPSK_mod;
		demodulate_function = &zBPSK_demod;
		initialized = true;
		double val = 1.0/sqrt(2);
		sym[0].re = val;
		sym[0].im = val;
		sym[1].re = -val;
		sym[1].im = val;

		sym[2].re = -val;
		sym[2].im = -val;
		sym[3].re = val;
		sym[3].im = -val;
	}
	~TPSK4() {
		if(sym) delete[] sym;
	}

};

extern TPSK4 psk4;
*/

#endif
