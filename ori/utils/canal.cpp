// qpskawgn.cpp : Defines the entry point for the console application.
//


#include  <math.h>


#include "canal.h"
#include "globaldef.h"

#ifdef _NSP_INTEL_

/** coef puissance decroissant */
DCplx *canal_nip(DCplx *data, int data_size, DCplx *coef, int coef_size) {
	DCplx *newdata=nspzMalloc(data_size);
	DCplx *tmpdata=nspzMalloc(coef_size);
	int iteration=data_size+coef_size-1;
	for(int i=0;i<data_size;i++) {
		if(i<coef_size) {
			nspzbZero(tmpdata, coef_size);
			nspzbCopy(data, (tmpdata+coef_size-i-1),i+1); // 0 0 x x x
		} else {
			nspzbCopy((data+i-coef_size+1), tmpdata, coef_size);
		}
		nspzbMpy2(coef, tmpdata, coef_size);
		nspzSum(tmpdata, coef_size, &(newdata[i]));
	}
	nspFree(tmpdata);
	return newdata;
}

/** coef puissance decroissant */
void canal(DCplx *data, int data_size, DCplx *coef, int coef_size) {
	DCplx *newdata=nspzMalloc(data_size);
	DCplx *tmpdata=nspzMalloc(coef_size);
	int iteration=data_size+coef_size-1;
	for(int i=0;i<data_size;i++) {
		if(i<coef_size) {
			nspzbZero(tmpdata, coef_size);
			nspzbCopy(data, (tmpdata+coef_size-i-1),i+1); // 0 0 x x x
		} else {
			nspzbCopy((data+i-coef_size+1), tmpdata, coef_size);
		}
		nspzbMpy2(coef, tmpdata, coef_size);
		nspzSum(tmpdata, coef_size, &(newdata[i]));
	}
	nspzbCopy(newdata, data, data_size);
	nspFree(tmpdata);
	nspFree(newdata);
}

void canalconv(DCplx *data, int data_size, DCplx *coef, int coef_size) {
	DCplx *newdata=nspzMalloc(data_size+coef_size-1);
	nspzConv(data, data_size, coef, coef_size ,newdata);
	nspzbCopy(newdata, data, data_size);
	nspFree(newdata);
}

DCplx *canalconv_nip(DCplx *data, int data_size, DCplx *coef, int coef_size) {
	DCplx *newdata=nspzMalloc(data_size+coef_size-1);
	nspzConv(data, data_size, coef, coef_size ,newdata);
	return newdata;
}

#else //_NSP_INTEL_
#include "utilitis.h"

ZVector canal_nip(const ZVector &data, const ZVector &coef) {
	return data.conv_nip(coef);
}

DVector canal_nip(const DVector &data, const DVector &coef) {
	return data.conv_nip(coef);
}

ZVector canal_nip(const DVector &data, const ZVector &coef) {
	return data.conv_nip(coef);
}

ZVector canal_nip(const ZVector &data, const DVector &coef) {
	return data.conv_nip(coef);
}

/** coef puissance decroissant */
void canal(ZVector &data, const ZVector &coef) {
	data.fixedconv(coef);
}

/** coef puissance decroissant */
void canal(DVector &data, const DVector &coef) {
	data.fixedconv(coef);
}


#endif //_NSP_INTEL_
