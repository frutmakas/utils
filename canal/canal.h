#ifndef _CANAL_H_
#define _CANAL_H_

#include "globaldef.h"
#include "tools/utilitis.h"

ZVector canal_nip(const ZVector &data, const ZVector &coef);
DVector canal_nip(const DVector &data, const DVector &coef);
ZVector canal_nip(const DVector &data, const ZVector &coef);
ZVector canal_nip(const ZVector &data, const DVector &coef);
void canal(ZVector &data, const ZVector &coef);
void canal(DVector &data, const DVector &coef);

#endif

