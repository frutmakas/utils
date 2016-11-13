/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2003/06/11 14:18:45 $
 * $Revision: 1.1.2.3 $
 * $Id: utils2lapack.h,v 1.1.2.3 2003/06/11 14:18:45 syed Exp $
 ********************************************************************/
#ifndef _UTILS_2_LAPCK_H_
#define _UTILS_2_LAPCK_H_
//#include "blaswrap.h"
#include "tools/all.h"
#include "clapack.h"


inline doublecomplex DCplxToDoubleComplex(const DCplx &m)  {
	doublecomplex c; 
	c.r = m.re; c.i=m.im; return c;
}

doublecomplex * ZMatrixToDoubleComplexVector(const ZMatrix &m) ;
doublereal * DMatrixToDoubleRealVector(const DMatrix &m);

ZMatrix & DoubleComplexVectorToZMatrix(doublecomplex *c, int line, int col, ZMatrix &m);
DMatrix & DoubleRealVectorToDMatrix(doublereal *c, int line, int col, DMatrix &m);

void printdoublecomplexasmatrix(ofstream &os, doublecomplex *u, int line, int col);
void printdoublecomplexasvector(ofstream &os, doublecomplex *u, int size);
#endif