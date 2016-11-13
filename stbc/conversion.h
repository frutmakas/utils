#ifndef _CONVERSION_STTC1_H_
#define _CONVERSION_STTC1_H_

#include "globaldef.h"
ZMatrix G2(const ZVector &ech);
ZMatrix G2(const ZVector &ech, int nb_ech);
ZMatrix G3(const ZVector &ech);
ZMatrix G3(const ZVector &ech, int nb_ech);
ZMatrix G4(const ZVector &ech);
ZMatrix G4(const ZVector &ech, int nb_ech);
ZMatrix G2_type2(const ZVector &ech, int nb_ech);

#endif