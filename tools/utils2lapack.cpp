/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2003/06/11 14:18:45 $
 * $Revision: 1.1.2.3 $
 * $Id: utils2lapack.cpp,v 1.1.2.3 2003/06/11 14:18:45 syed Exp $
 ********************************************************************/

#include "tools/all.h"
#include "tools/utils2lapack.h"
#include <iostream.h>

doublecomplex * ZMatrixToDoubleComplexVector(const ZMatrix &m) {
	// lapack represent the matrix in an array packed by columns
	if(m.line<=0 || m.col<=0 || m.mat==NULL ) throw CUtilitisException(INVALID_MATRIX);
	doublecomplex *c = new doublecomplex[m.line*m.col];
	if(c==NULL) throw CUtilitisException(OUT_OF_MEMORY);
	register int i,j;
	for(i=0;i<m.col; i++)
		for(j=0;j<m.line;j++)  {
			//c[i*m.line+j]=DCplxToDoubleComplex(m.mat[i][j]);
			c[i*m.line+j].r=m.mat[j][i].re; 
			c[i*m.line+j].i=m.mat[j][i].im;
		}
	return c;
}

ZMatrix & DoubleComplexVectorToZMatrix(doublecomplex *c, int line, int col, ZMatrix &m) {
	if(line<=0 || col<=0 || c==NULL) throw CUtilitisException(INVALID_VECTOR);
	m = ZMatrix(line,col);
	register i,j;
	for(i=0;i<col;i++) 
		for(j=0;j<line;j++) {
			m.mat[j][i].re=c[i*line+j].r;
			m.mat[j][i].im=c[i*line+j].i;
		}
	return m;
}

doublereal * DMatrixToDoubleRealVector(const DMatrix &m) {
	// lapack represent the matrix in an array packed by columns
	if(m.line<=0 || m.col<=0 || m.mat==NULL ) throw CUtilitisException(INVALID_MATRIX);
	doublereal *c = new doublereal[m.line*m.col];
	if(c==NULL) throw CUtilitisException(OUT_OF_MEMORY);
	register int i,j;
	for(i=0;i<m.col; i++)
		for(j=0;j<m.line;j++)  {
			c[i*m.line+j]=m.mat[j][i];
//			c[i*m.line+j].r=m.mat[j][i].re; 
//			c[i*m.line+j].i=m.mat[j][i].im;
		}
	return c;
}


DMatrix & DoubleRealVectorToDMatrix(doublereal *c, int line, int col, DMatrix &m) {
	if(line<=0 || col<=0 || c==NULL) throw CUtilitisException(INVALID_VECTOR);
	m = DMatrix(line,col);
	register i,j;
	for(i=0;i<col;i++) 
		for(j=0;j<line;j++) 
			m.mat[j][i]=c[i*line+j];
	return m;
}

void printdoublecomplexasmatrix(ofstream &os, doublecomplex *u, int line, int col){
	os << endl << "["; 
	for(int i=0; i < line; i++) { 
		os << "[" ;
		for(int j=0; j< col; j++) {
			os << u[i*col+j].r << "+" << u[i*col+j].i<< "i " ;
		}
		os << "]"<< endl;
	}
        os << "]" << endl;
}


void printdoublecomplexasvector(ofstream &os, doublecomplex *u, int size){
	os << endl << "["; 
	for(int i=0; i < size; i++) os << u[i].r<<"+"<<u[i].i<<"i ";
	os << "]"<<endl;
}