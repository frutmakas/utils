/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/04/05 18:49:16 $
 * $Revision: 1.1.2.5 $
 * $Id: fft.h,v 1.1.2.5 2004/04/05 18:49:16 syed Exp $
 ********************************************************************/
 
#ifndef _FFT_H_
#define _FFT_H_

#include "tools/utilitis.h"
#define SHARED_NORMALIZATION
void four1(double  data[], unsigned long nn, int isign);

ZVector fft(const ZVector &td);

ZVector ifft(const ZVector &fd);

ZVector zero_center(const ZVector &ffted);

ZVector butterfly(const ZVector &zero_centered);

#endif // _FFT_H_

