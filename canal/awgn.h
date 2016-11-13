#ifndef _AWGN_H_
#define _AWGN_H_

#include "globaldef.h"

#include "tools/utilitis.h"
#include "rand/randgen.h"
#include "tools/myutils.h"

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
