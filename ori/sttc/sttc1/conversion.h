#ifndef _CONVERSION_STTC1_H_
#define _CONVERSION_STTC1_H_

#include "globaldef.h"

#ifdef _NSP_INTEL_

#include <nsp.h>
#include "nspmatrice.h"
#include "nspdcplx.h"

ZMatrix G2(DCplx *ech);
ZMatrix G3(DCplx *ech);
ZMatrix G4(DCplx *ech);

ZMatrix G2(DCplx *ech, int nb_ech);
ZMatrix G3(DCplx *ech, int nb_ech);
ZMatrix G4(DCplx *ech, int nb_ech);

#else //_NSP_INTEL
ZMatrix G2(const ZVector &ech);
ZMatrix G2(const ZVector &ech, int nb_ech);
ZMatrix G3(const ZVector &ech);
ZMatrix G3(const ZVector &ech, int nb_ech);
ZMatrix G4(const ZVector &ech);
ZMatrix G4(const ZVector &ech, int nb_ech);

#endif //_NSP_INTEL_

#endif