#ifndef _GOLD_H_
#define _GOLD_H_
#include "globaldef.h"

#ifdef _NSP_INTEL_

#include <nsp.h>

double *dGoldCode(int InitialValue[2], int Polynomial[2], int polydeg, int idxwidth, int idx, int *code_size);
DCplx *zGoldCode(int InitialValue[2], int Polynomial[2], int polydeg, int idxwidth, int idx, int *code_size);

DCplx *zSpreadData(DCplx *data, int data_size, DCplx *spreadcode, int spread_size, int *spreaddata_size);
DCplx *zSpreadData(DCplx *data, int data_size, double *spreadcode, int spread_size, int *spreaddata_size);
double *dSpreadData(double *data, int data_size, double *spreadcode, int spread_size, int *spreaddata_size);
double *dUnspreadData(double *data, int data_size, double *spreadcode, int spread_size, int *unspreaded_size);

#else 

DVector dGoldCode(int InitialValue[2], int Polynomial[2], int polydeg, int idxwidth, int idx);
ZVector zGoldCode(int InitialValue[2], int Polynomial[2], int polydeg, int idxwidth, int idx);

ZVector SpreadData(const ZVector &data, const ZVector &spreadcode);
ZVector SpreadData(const ZVector &data, const DVector &spreadcode);
ZVector SpreadData(const DVector &data, const ZVector &spreadcode);
DVector zSpreadData(const DVector &data, const DVector &spreadcode);

DVector UnspreadData(const DVector &data, const DVector &spreadcode);

#endif
#endif;
