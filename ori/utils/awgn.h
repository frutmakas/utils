#ifndef _AWGN_H_
#define _AWGN_H_

#include "globaldef.h"
#ifdef _NSP_INTEL_

#define nsp_UsesSampleGen
#define nsp_UsesVector

#include "nsp.h"


extern NSPZRandGausState zgausStatePtr;
extern NSPDRandGausState dgausStatePtr;
extern NSPWRandUniState  uniStatePtr;

DCplx zAWGN_nip(DCplx data);
void zAWGN(DCplx *data);
DCplx *zbAWGN_nip(DCplx *data, int size);
void zbAWGN(DCplx *data, int size);


double dAWGN_nip(double data) ;
void dAWGN(double *data);
double *dbAWGN_nip(double *data, int size);
void dbAWGN(double *data, int size);

#else // _NSP_INTEL_

#include "utilitis.h"
#include "randgen.h"
#include "myutils.h"

extern zRandGausStatePtr zgausStatePtr;
extern dRandGausStatePtr  dgausStatePtr;
//extern zRandGausStatePtr   uniStatePtr;

DCplx zAWGN_nip(DCplx data) ;
void zAWGN(DCplx *data);
ZVector zbAWGN_nip(const ZVector &data);
void zbAWGN(const ZVector &data);
ZVector zbAWGN_nip(const DVector &data);
double dAWGN_nip(double data) ;
void dAWGN(double *data);
DVector dbAWGN_nip(const DVector &data);
ZVector dbAWGN_nip(const ZVector &data);
void dbAWGN(const DVector &data);
void dbAWGN(const ZVector &data);

#endif



#endif