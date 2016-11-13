/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/01/20 13:39:52 $
 * $Revision: 1.1.2.1 $
 * $Id: gold.h,v 1.1.2.1 2004/01/20 13:39:52 syed Exp $
 ********************************************************************/

#ifndef _GOLD_H_
#define _GOLD_H_
#include "globaldef.h"

DVector dGoldCode(int InitialValue[2], int Polynomial[2], int polydeg, int idxwidth, int idx, bool normalize=true);
ZVector zGoldCode(int InitialValue[2], int Polynomial[2], int polydeg, int idxwidth, int idx, bool normalize=true);
#endif
