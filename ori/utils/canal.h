#ifndef _CANAL_H_
#define _CANAL_H_

#include "globaldef.h"

#ifdef _NSP_INTEL_

#define nsp_UsesSampleGen
#define nsp_UsesVector
#define nsp_UsesConvolution

#include "nsp.h"

DCplx *canal_nip(DCplx *data, int data_size, DCplx *coef, int coef_size) ;
void canal(DCplx *data, int data_size, DCplx *coef, int coef_size) ;
void canalconv(DCplx *data, int data_size, DCplx *coef, int coef_size) ;
DCplx *canalconv_nip(DCplx *data, int data_size, DCplx *coef, int coef_size) ;

#else 
#include "utilitis.h"

ZVector canal_nip(const ZVector &data, const ZVector &coef);
DVector canal_nip(const DVector &data, const DVector &coef);
ZVector canal_nip(const DVector &data, const ZVector &coef);
ZVector canal_nip(const ZVector &data, const DVector &coef);
void canal(ZVector &data, const ZVector &coef);
void canal(DVector &data, const DVector &coef);

#endif
#endif