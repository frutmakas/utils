/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/02/11 21:48:37 $
 * $Revision: 1.1.2.2 $
 * $Id: spread.h,v 1.1.2.2 2004/02/11 21:48:37 syed Exp $
 ********************************************************************/
#ifndef __CDMA__SPREAD_H_
#define __CDMA__SPREAD_H_

ZVector SpreadData(const ZVector &data, const ZVector &spreadcode);
ZVector SpreadData(const ZVector &data, const DVector &spreadcode);
ZVector SpreadData(const DVector &data, const ZVector &spreadcode);
DVector SpreadData(const DVector &data, const DVector &spreadcode);

DVector UnspreadData(const DVector &data, const DVector &spreadcode);
ZVector UnspreadData(const ZVector &data, const ZVector &spreadcode);
ZVector UnspreadData(const ZVector &data, const DVector &spreadcode);
ZVector UnspreadData(const DVector &data, const ZVector &spreadcode);


#endif

