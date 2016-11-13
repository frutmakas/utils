/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/06/23 09:05:41 $
 * $Revision: 1.1.2.24 $
 * $Id: ext_utilitis.cpp,v 1.1.2.24 2006/06/23 09:05:41 syed Exp $
 ********************************************************************/

//#include <fstream.h>
#include <iostream>
#include <fstream>
using namespace std;
#include "tools/utilitis.h"
#include "tools/superclass.h"
#include "tools/ext_utilitis.h"
#include "globaldef.h"

UVector<ZVector> SplitMatrixByCol(const ZMatrix &m) { // on aura col vecteurs de taille line
	if (m.mat==NULL || m.line<=0 || m.col<=0) throw CExtUtilitisException(INVALID_MATRIX);
	UVector<ZVector> u(m.col);
	for(register int i=0;i<m.col;i++) {
		u.vect[i]=ZVector(m.line);
		if(u.vect[i].vect==NULL) throw CExtUtilitisException(MEMORY_ALLOCATION_ERROR);
		for(register int j=0;j<m.line;j++) {
			u.vect[i].vect[j]=m.mat[j][i];
		}
	}
	return u;
}

UVector<ZVector> SplitMatrixByLine(const ZMatrix &m) {
	if (m.mat==NULL || m.line<=0 || m.col<=0) throw CExtUtilitisException(INVALID_MATRIX);
	UVector<ZVector> u(m.line);
	for(register int i=0;i<m.line;i++) {
		u.vect[i]=ZVector(m.col);
		if(u.vect[i].vect==NULL) throw CExtUtilitisException(MEMORY_ALLOCATION_ERROR);
		for(register int j=0;j<m.col;j++) {
			u.vect[i].vect[j]=m.mat[i][j];
		}
	}
	return u;
}

UVector<DVector> SplitMatrixByCol(const DMatrix &m) {
	if (m.mat==NULL || m.line<=0 || m.col<=0) throw CExtUtilitisException(INVALID_MATRIX);
	UVector<DVector> u(m.col);
	for(register int i=0;i<m.col;i++) {
		u.vect[i]=DVector(m.line);
		if(u.vect[i].vect==NULL) throw CExtUtilitisException(MEMORY_ALLOCATION_ERROR);
		for(register int j=0;j<m.line;j++) {
			u.vect[i].vect[j]=m.mat[j][i];
		}
	}
	return u;
}

UVector<DVector> SplitMatrixByLine(const DMatrix &m) {
	if (m.mat==NULL || m.line<=0 || m.col<=0) throw CExtUtilitisException(INVALID_MATRIX);
	UVector<DVector> u(m.line);
	for(register int i=0;i<m.line;i++) {
		u.vect[i]=DVector(m.col);
		if(u.vect[i].vect==NULL) throw CExtUtilitisException(MEMORY_ALLOCATION_ERROR);
		for(register int j=0;j<m.col;j++) {
			u.vect[i].vect[j]=m.mat[i][j];
		}
	}
	return u;
}

UVector<WVector> SplitMatrixByCol(const WMatrix &m) {
	if (m.mat==NULL || m.line<=0 || m.col<=0) throw CExtUtilitisException(INVALID_MATRIX);
	UVector<WVector> u(m.col);
	for(register int i=0;i<m.col;i++) {
		u.vect[i]=WVector(m.line);
		if(u.vect[i].vect==NULL) throw CExtUtilitisException(MEMORY_ALLOCATION_ERROR);
		for(register int j=0;j<m.line;j++) {
			u.vect[i].vect[j]=m.mat[j][i];
		}
	}
	return u;
}

UVector<WVector> SplitMatrixByLine(const WMatrix &m) {
	if (m.mat==NULL || m.line<=0 || m.col<=0) throw CExtUtilitisException(INVALID_MATRIX);
	UVector<WVector> u(m.line);
	for(register int i=0;i<m.line;i++) {
		u.vect[i]=WVector(m.col);
		if(u.vect[i].vect==NULL) throw CExtUtilitisException(MEMORY_ALLOCATION_ERROR);
		for(register int j=0;j<m.col;j++) {
			u.vect[i].vect[j]=m.mat[i][j];
		}
	}
	return u;
}

ZMatrix JoinUVectorInCol(const UVector<ZVector> &u) { // on aura une matrice de taille u.vect[0].taille x u.taille
	if(u.taille<=0 || u.vect==NULL) throw CExtUtilitisException(INVALID_UVECTOR);
	//if(u.vect[0].taille<=0 || u.vect[0].vect==NULL) CExtUtilitisException(EMPTY_UVECTOR);
	ZMatrix m(u.vect[0].taille, u.taille);
	for(register int i=0;i<u.taille;i++) {
		for(register int j=0;j<u.vect[i].taille;j++) {
			if(j>m.line) throw CExtUtilitisException(INCONSISTENT_UVECTOR);
			m.mat[j][i]=u.vect[i].vect[j];
		}
	}
	return m;
}

ZMatrix JoinUVectorInLine(const UVector<ZVector> &u) { // on aura une matrice de taille u.vect[0].taille x u.taille
	if(u.taille<=0 || u.vect==NULL) throw CExtUtilitisException(INVALID_UVECTOR);
	//if(u.vect[0].taille<=0 || u.vect[0].vect==NULL) CExtUtilitisException(EMPTY_UVECTOR);
	ZMatrix m(u.taille, u.vect[0].taille);
	for(register int i=0;i<u.taille;i++) {
		for(register int j=0;j<u.vect[i].taille;j++) {
			if(j>m.col) throw CExtUtilitisException(INCONSISTENT_UVECTOR);
			m.mat[i][j]=u.vect[i].vect[j];
		}
	}
	return m;
}

DMatrix JoinUVectorInCol(const UVector<DVector> &u) { // on aura une matrice de taille u.vect[0].taille x u.taille
	if(u.taille<=0 || u.vect==NULL) throw CExtUtilitisException(INVALID_UVECTOR);
	//if(u.vect[0].taille<=0 || u.vect[0].vect==NULL) CExtUtilitisException(EMPTY_UVECTOR);
	DMatrix m(u.vect[0].taille, u.taille);
	for(register int i=0;i<u.taille;i++) {
		for(register int j=0;j<u.vect[i].taille;j++) {
			if(j>m.line) throw CExtUtilitisException(INCONSISTENT_UVECTOR);
			m.mat[j][i]=u.vect[i].vect[j];
		}
	}
	return m;
}

DMatrix JoinUVectorInLine(const UVector<DVector> &u) { // on aura une matrice de taille u.vect[0].taille x u.taille
	if(u.taille<=0 || u.vect==NULL) throw CExtUtilitisException(INVALID_UVECTOR);
	//if(u.vect[0].taille<=0 || u.vect[0].vect==NULL) CExtUtilitisException(EMPTY_UVECTOR);
	DMatrix m(u.taille, u.vect[0].taille);
	for(register int i=0;i<u.taille;i++) {
		for(register int j=0;j<u.vect[i].taille;j++) {
			if(j>m.col) throw CExtUtilitisException(INCONSISTENT_UVECTOR);
			m.mat[i][j]=u.vect[i].vect[j];
		}
	}
	return m;
}

WMatrix JoinUVectorInCol(const UVector<WVector> &u) { // on aura une matrice de taille u.vect[0].taille x u.taille
	if(u.taille<=0 || u.vect==NULL) throw CExtUtilitisException(INVALID_UVECTOR);
	//if(u.vect[0].taille<=0 || u.vect[0].vect==NULL) CExtUtilitisException(EMPTY_UVECTOR);
	WMatrix m(u.vect[0].taille, u.taille);
	for(register int i=0;i<u.taille;i++) {
		for(register int j=0;j<u.vect[i].taille;j++) {
			if(j>m.line) throw CExtUtilitisException(INCONSISTENT_UVECTOR);
			m.mat[j][i]=u.vect[i].vect[j];
		}
	}
	return m;
}

WMatrix JoinUVectorInLine(const UVector<WVector> &u) { // on aura une matrice de taille u.vect[0].taille x u.taille
	if(u.taille<=0 || u.vect==NULL) throw CExtUtilitisException(INVALID_UVECTOR);
	//if(u.vect[0].taille<=0 || u.vect[0].vect==NULL) CExtUtilitisException(EMPTY_UVECTOR);
	WMatrix m(u.taille, u.vect[0].taille);
	for(register int i=0;i<u.taille;i++) {
		for(register int j=0;j<u.vect[i].taille;j++) {
			if(j>m.col) throw CExtUtilitisException(INCONSISTENT_UVECTOR);
			m.mat[i][j]=u.vect[i].vect[j];
		}
	}
	return m;
}

ZMatrix VectorInColMatrixRepresentation(const ZVector &zv) { // on aura une matrice de taille zv.taille x 1
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	ZMatrix zm(zv.taille,1);
	for(register int i=0;i<zv.taille;i++) {
		zm.mat[i][0]=zv.vect[i];
	}
	return zm;
}

ZMatrix VectorInLineMatrixRepresentation(const ZVector &zv) { // on aura une matrice de taille 1 x zv.taille 
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	ZMatrix zm(1, zv.taille);
	for(register int i=0;i<zv.taille;i++) {
		zm.mat[0][i]=zv.vect[i];
	}
	return zm;
}

DMatrix VectorInColMatrixRepresentation(const DVector &zv) { // on aura une matrice de taille zv.taille x 1
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	DMatrix zm(zv.taille,1);
	for(register int i=0;i<zv.taille;i++) {
		zm.mat[i][0]=zv.vect[i];
	}
	return zm;
}

DMatrix VectorInLineMatrixRepresentation(const DVector &zv) { // on aura une matrice de taille 1 x zv.taille 
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	DMatrix zm(1, zv.taille);
	for(register int i=0;i<zv.taille;i++) {
		zm.mat[0][i]=zv.vect[i];
	}
	return zm;
}

WMatrix VectorInColMatrixRepresentation(const WVector &zv) { // on aura une matrice de taille zv.taille x 1
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	WMatrix zm(zv.taille,1);
	for(register int i=0;i<zv.taille;i++) {
		zm.mat[i][0]=zv.vect[i];
	}
	return zm;
}

WMatrix VectorInLineMatrixRepresentation(const WVector &zv) { // on aura une matrice de taille 1 x zv.taille 
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	WMatrix zm(1, zv.taille);
	for(register int i=0;i<zv.taille;i++) {
		zm.mat[0][i]=zv.vect[i];
	}
	return zm;
}

ZVector MatrixInVectorRepresentation(const ZMatrix &zm) {
	if(zm.mat==NULL || zm.line<=0 || zm.col<=0) throw CExtUtilitisException(INVALID_MATRIX);
	
	if(zm.line==1) {
		ZVector zv(zm.col);
		for(register int i =0; i <zm.col;i++) zv.vect[i]=zm.mat[0][i];
		return zv;
	} 
	if(zm.col==1){
		ZVector zv(zm.line);
		for(register int i =0; i <zm.line;i++) zv.vect[i]=zm.mat[i][0];
		return zv;
	}
	throw  CExtUtilitisException(INVALID_MATRIX);
}

DVector MatrixInVectorRepresentation(const DMatrix &zm) {
	if(zm.mat==NULL || zm.line<=0 || zm.col<=0) throw CExtUtilitisException(INVALID_MATRIX);
	
	if(zm.line==1) {
		DVector zv(zm.col);
		for(register int i =0; i <zm.col;i++) zv.vect[i]=zm.mat[0][i];
		return zv;
	} 
	if(zm.col==1){
		DVector zv(zm.line);
		for(register int i =0; i <zm.line;i++) zv.vect[i]=zm.mat[i][0];
		return zv;
	}
	throw  CExtUtilitisException(INVALID_MATRIX);
}

WVector MatrixInVectorRepresentation(const WMatrix &zm) {
	if(zm.mat==NULL || zm.line<=0 || zm.col<=0) throw CExtUtilitisException(INVALID_MATRIX);
	
	if(zm.line==1) {
		WVector zv(zm.col);
		for(register int i =0; i <zm.col;i++) zv.vect[i]=zm.mat[0][i];
		return zv;
	} 
	if(zm.col==1){
		WVector zv(zm.line);
		for(register int i =0; i <zm.line;i++) zv.vect[i]=zm.mat[i][0];
		return zv;
	}
	throw  CExtUtilitisException(INVALID_MATRIX);
}


ZVector MatrixAsColumnPackedVector(const ZMatrix &z) {
	if(z.line<=0 || z.col <=0 || z.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	ZVector v(z.line*z.col);
	for(register int i=0;i<z.col;i++){
		for(register int j=0;j<z.line;j++) {
			v.vect[i*z.line+j]=z.mat[j][i];
		}
	}
	return v;
}

ZVector MatrixAsLinePackedVector(const ZMatrix &z) {
	if(z.line<=0 || z.col <=0 || z.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	ZVector v(z.line*z.col);
	for(register int i=0;i<z.line;i++){
		for(register int j=0;j<z.col;j++) {
			v.vect[i*z.col+j]=z.mat[i][j];
		}
	}
	return v;
}

DVector MatrixAsColumnPackedVector(const DMatrix &z) {
	if(z.line<=0 || z.col <=0 || z.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DVector v(z.line*z.col);
	for(register int i=0;i<z.col;i++){
		for(register int j=0;j<z.line;j++) {
			v.vect[i*z.line+j]=z.mat[j][i];
		}
	}
	return v;
}

DVector MatrixAsLinePackedVector(const DMatrix &z) {
	if(z.line<=0 || z.col <=0 || z.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DVector v(z.line*z.col);
	for(register int i=0;i<z.line;i++){
		for(register int j=0;j<z.col;j++) {
			v.vect[i*z.col+j]=z.mat[i][j];
		}
	}
	return v;
}

WVector MatrixAsColumnPackedVector(const WMatrix &z) {
	if(z.line<=0 || z.col <=0 || z.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	WVector v(z.line*z.col);
	for(register int i=0;i<z.col;i++){
		for(register int j=0;j<z.line;j++) {
			v.vect[i*z.line+j]=z.mat[j][i];
		}
	}
	return v;
}

WVector MatrixAsLinePackedVector(const WMatrix &z) {
	if(z.line<=0 || z.col <=0 || z.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	WVector v(z.line*z.col);
	for(register int i=0;i<z.line;i++){
		for(register int j=0;j<z.col;j++) {
			v.vect[i*z.col+j]=z.mat[i][j];
		}
	}
	return v;
}

ZMatrix ColumnPackedVectorAsMatrix(const ZVector &z, int line, int col) {
	if(z.taille<=0 || z.vect ==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(line*col!=z.taille) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	ZMatrix m(line,col);
	register int i,j,k=0;
	for(i=0;i<col;i++) {
		for(j=0;j<line;j++) {
			m.mat[j][i]=z.vect[k++];
		}
	}
	return m;
}

ZMatrix LinePackedVectorAsMatrix(const ZVector &z, int line, int col) {
	if(z.taille<=0 || z.vect ==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(line*col!=z.taille) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	ZMatrix m(line,col);
	register int i,j,k=0;
	for(i=0;i<line;i++) {
		for(j=0;j<col;j++) {
			m.mat[i][j]=z.vect[k++];
		}
	}
	return m;
}

DMatrix ColumnPackedVectorAsMatrix(const DVector &z, int line, int col) {
	if(z.taille<=0 || z.vect ==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(line*col!=z.taille) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	DMatrix m(line,col);
	register int i,j,k=0;
	for(i=0;i<col;i++) {
		for(j=0;j<line;j++) {
			m.mat[j][i]=z.vect[k++];
		}
	}
	return m;
}

DMatrix LinePackedVectorAsMatrix(const DVector &z, int line, int col) {
	if(z.taille<=0 || z.vect ==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(line*col!=z.taille) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	DMatrix m(line,col);
	register int i,j,k=0;
	for(i=0;i<line;i++) {
		for(j=0;j<col;j++) {
			m.mat[i][j]=z.vect[k++];
		}
	}
	return m;
}

WMatrix ColumnPackedVectorAsMatrix(const WVector &z, int line, int col) {
	if(z.taille<=0 || z.vect ==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(line*col!=z.taille) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	WMatrix m(line,col);
	register int i,j,k=0;
	for(i=0;i<col;i++) {
		for(j=0;j<line;j++) {
			m.mat[j][i]=z.vect[k++];
		}
	}
	return m;
}

WMatrix LinePackedVectorAsMatrix(const WVector &z, int line, int col) {
	if(z.taille<=0 || z.vect ==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(line*col!=z.taille) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	WMatrix m(line,col);
	register int i,j,k=0;
	for(i=0;i<line;i++) {
		for(j=0;j<col;j++) {
			m.mat[i][j]=z.vect[k++];
		}
	}
	return m;
}


WVector ColXInWMatrixToVector(const WMatrix &zm, int col) { // col = [0;zm.col[
	if(zm.mat==NULL || zm.line<=0 || zm.col<=0) throw CExtUtilitisException(INVALID_MATRIX);
	if(zm.col<=col) throw CExtUtilitisException(MATRIX_INDEX_OUT_OF_RANGE);
	WVector zv(zm.line);
	for(register int i =0; i <zm.line;i++) zv.vect[i]=zm.mat[i][col];
	return zv;
}

WVector LineXInWMatrixToVector(const WMatrix &zm, int line) { // line = [0;zm.line[
	if(zm.mat==NULL || zm.line<=0 || zm.col<=0) throw CExtUtilitisException(INVALID_MATRIX);
	if(zm.line<=line) throw CExtUtilitisException(MATRIX_INDEX_OUT_OF_RANGE);
	WVector zv(zm.col);
	for(register int i =0; i <zm.col;i++) zv.vect[i]=zm.mat[line][i];
	return zv;
}

DVector ColXInDMatrixToVector(const DMatrix &zm, int col) { // col = [0;zm.col[
	if(zm.mat==NULL || zm.line<=0 || zm.col<=0) throw CExtUtilitisException(INVALID_MATRIX);
	if(zm.col<=col) throw CExtUtilitisException(MATRIX_INDEX_OUT_OF_RANGE);
	DVector zv(zm.line);
	for(register int i =0; i <zm.line;i++) zv.vect[i]=zm.mat[i][col];
	return zv;
}

DVector LineXInDMatrixToVector(const DMatrix &zm, int line) { // line = [0;zm.line[
	if(zm.mat==NULL || zm.line<=0 || zm.col<=0) throw CExtUtilitisException(INVALID_MATRIX);
	if(zm.line<=line) throw CExtUtilitisException(MATRIX_INDEX_OUT_OF_RANGE);
	DVector zv(zm.col);
	for(register int i =0; i <zm.col;i++) zv.vect[i]=zm.mat[line][i];
	return zv;
}

ZVector ColXInZMatrixToVector(const ZMatrix &zm, int col) { // col = [0;zm.col[
	if(zm.mat==NULL || zm.line<=0 || zm.col<=0) throw CExtUtilitisException(INVALID_MATRIX);
	if(zm.col<=col) throw CExtUtilitisException(MATRIX_INDEX_OUT_OF_RANGE);
	ZVector zv(zm.line);
	for(register int i =0; i <zm.line;i++) zv.vect[i]=zm.mat[i][col];
	return zv;
}

ZVector LineXInZMatrixToVector(const ZMatrix &zm, int line) { // line = [0;zm.line[
	if(zm.mat==NULL || zm.line<=0 || zm.col<=0) throw CExtUtilitisException(INVALID_MATRIX);
	if(zm.line<=line) throw CExtUtilitisException(MATRIX_INDEX_OUT_OF_RANGE);
	ZVector zv(zm.col);
	for(register int i =0; i <zm.col;i++) zv.vect[i]=zm.mat[line][i];
	return zv;
}

ZMatrix ZVectorToDiagonalMatrix(const ZVector &zv) {
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	ZMatrix zm(zv.taille,zv.taille);
	for(register int i=0;i<zv.taille;i++) zm.mat[i][i]=zv.vect[i];
	return zm;
}

DMatrix DVectorToDiagonalMatrix(const DVector &zv) {
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	DMatrix zm(zv.taille,zv.taille);
	for(register int i=0;i<zv.taille;i++) zm.mat[i][i]=zv.vect[i];
	return zm;
}

WMatrix WVectorToDiagonalMatrix(const WVector &zv) {
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	WMatrix zm(zv.taille,zv.taille);
	for(register int i=0;i<zv.taille;i++) zm.mat[i][i]=zv.vect[i];
	return zm;
}

ZMatrix LineXinMatrixAsDiagonal(const ZMatrix &m, int x) { 
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CExtUtilitisException(INVALID_MATRIX);
	if(x>=m.col) throw CExtUtilitisException(INDEX_OUT_OF_RANGE);
	ZMatrix n(m.col, m.col);
	for(int i=0;i<m.col;i++) {
		n.mat[i][i]=m.mat[x][i];
	}
	return n;
}

DMatrix LineXinMatrixAsDiagonal(const DMatrix &m, int x) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CExtUtilitisException(INVALID_MATRIX);
	if(x>=m.col) throw CExtUtilitisException(INDEX_OUT_OF_RANGE);
	DMatrix n(m.col, m.col);
	for(int i=0;i<m.col;i++) {
		n.mat[i][i]=m.mat[x][i];
	}
	return n;
}

WMatrix LineXinMatrixAsDiagonal(const WMatrix &m, int x) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CExtUtilitisException(INVALID_MATRIX);
	if(x>=m.col) throw CExtUtilitisException(INDEX_OUT_OF_RANGE);
	WMatrix n(m.col, m.col);
	for(int i=0;i<m.col;i++) {
		n.mat[i][i]=m.mat[x][i];
	}
	return n;
}

ZMatrix ColXinMatrixAsDiagonal(const ZMatrix &m, int x) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CExtUtilitisException(INVALID_MATRIX);
	if(x>=m.line) throw CExtUtilitisException(INDEX_OUT_OF_RANGE);
	ZMatrix n(m.line, m.line);
	for(int i=0;i<m.line;i++) {
		n.mat[i][i]=m.mat[i][x];
	}
	return n;
}

DMatrix ColXinMatrixAsDiagonal(const DMatrix &m, int x) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CExtUtilitisException(INVALID_MATRIX);
	if(x>=m.line) throw CExtUtilitisException(INDEX_OUT_OF_RANGE);
	DMatrix n(m.line, m.line);
	for(int i=0;i<m.col;i++) {
		n.mat[i][i]=m.mat[i][x];
	}
	return n;
}

WMatrix ColXinMatrixAsDiagonal(const WMatrix &m, int x) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CExtUtilitisException(INVALID_MATRIX);
	if(x>=m.line) throw CExtUtilitisException(INDEX_OUT_OF_RANGE);
	WMatrix n(m.line, m.line);
	for(int i=0;i<m.col;i++) {
		n.mat[i][i]=m.mat[i][x];
	}
	return n;
}

ZVector extractdiag(const ZMatrix &m) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CExtUtilitisException(INVALID_MATRIX);
	if(m.line != m.col) throw CExtUtilitisException(MATRIX_NOT_SQUARE);
	ZVector v(m.line);
	for(int i=0;i<m.line;i++) {
		v.vect[i]=m.mat[i][i];
	}
	return v;
}

WVector extractdiag(const WMatrix &m) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CExtUtilitisException(INVALID_MATRIX);
	if(m.line != m.col) throw CExtUtilitisException(MATRIX_NOT_SQUARE);
	WVector v(m.line);
	for(int i=0;i<m.line;i++) {
		v.vect[i]=m.mat[i][i];
	}
	return v;
}

DVector extractdiag(const DMatrix &m) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CExtUtilitisException(INVALID_MATRIX);
	if(m.line != m.col) throw CExtUtilitisException(MATRIX_NOT_SQUARE);
	DVector v(m.line);
	for(int i=0;i<m.line;i++) {
		v.vect[i]=m.mat[i][i];
	}
	return v;
}

ZMatrix diag(const ZVector &zv) {
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	ZMatrix z(zv.taille,zv.taille);
	for(int d=0;d<zv.taille;d++) {
		z.mat[d][d]=zv.vect[d];
	}
	return z;
}

DMatrix diag(const DVector &zv) {
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	DMatrix z(zv.taille,zv.taille);
	for(int d=0;d<zv.taille;d++) {
		z.mat[d][d]=zv.vect[d];
	}
	return z;
}

WMatrix diag(const WVector &zv) {
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	WMatrix z(zv.taille,zv.taille);
	for(int d=0;d<zv.taille;d++) {
		z.mat[d][d]=zv.vect[d];
	}
	return z;
}

ZMatrix diag(const UVector<ZMatrix> &zv) {
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	if(zv.vect[0].mat==NULL || zv.vect[0].line<=0 || zv.vect[0].col<=0) 
		throw CExtUtilitisException(INVALID_MATRIX);
	ZMatrix z(zv.taille*zv.vect[0].line,zv.taille*zv.vect[0].col);
	for(int d=0;d<zv.taille;d++) {
		if(zv.vect[d].mat==NULL || zv.vect[d].line<=0|| zv.vect[d].col<=0)
			throw CExtUtilitisException(INVALID_MATRIX);
		if(zv.vect[d].line!=zv.vect[0].line|| zv.vect[d].col!=zv.vect[0].line)
			throw CExtUtilitisException(INCONSISTENT_MATRIX_SIZE_IN_UVECTOR);
		for(int l=0;l<zv.vect[d].line;l++) {
			for(int c=0;c<zv.vect[d].col;c++) {
				z.mat[d*zv.vect[d].line+l][d*zv.vect[d].col+c]=zv.vect[d].mat[l][c];
			}
		}
	}
	return z;
}

DMatrix diag(const UVector<DMatrix> &zv) {
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	if(zv.vect[0].mat==NULL || zv.vect[0].line<=0 || zv.vect[0].col<=0) 
		throw CExtUtilitisException(INVALID_MATRIX);
	DMatrix z(zv.taille*zv.vect[0].line,zv.taille*zv.vect[0].col);
	for(int d=0;d<zv.taille;d++) {
		if(zv.vect[d].mat==NULL || zv.vect[d].line<=0|| zv.vect[d].col<=0)
			throw CExtUtilitisException(INVALID_MATRIX);
		if(zv.vect[d].line!=zv.vect[0].line|| zv.vect[d].col!=zv.vect[0].line)
			throw CExtUtilitisException(INCONSISTENT_MATRIX_SIZE_IN_UVECTOR);
		for(int l=0;l<zv.vect[d].line;l++) {
			for(int c=0;c<zv.vect[d].col;c++) {
				z.mat[d*zv.vect[d].line+l][d*zv.vect[d].col+c]=zv.vect[d].mat[l][c];
			}
		}
	}
	return z;
}

WMatrix diag(const UVector<WMatrix> &zv) {
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	if(zv.vect[0].mat==NULL || zv.vect[0].line<=0 || zv.vect[0].col<=0) 
		throw CExtUtilitisException(INVALID_MATRIX);
	WMatrix z(zv.taille*zv.vect[0].line,zv.taille*zv.vect[0].col);
	for(int d=0;d<zv.taille;d++) {
		if(zv.vect[d].mat==NULL || zv.vect[d].line<=0|| zv.vect[d].col<=0)
			throw CExtUtilitisException(INVALID_MATRIX);
		if(zv.vect[d].line!=zv.vect[0].line|| zv.vect[d].col!=zv.vect[0].line)
			throw CExtUtilitisException(INCONSISTENT_MATRIX_SIZE_IN_UVECTOR);
		for(int l=0;l<zv.vect[d].line;l++) {
			for(int c=0;c<zv.vect[d].col;c++) {
				z.mat[d*zv.vect[d].line+l][d*zv.vect[d].col+c]=zv.vect[d].mat[l][c];
			}
		}
	}
	return z;
}

ZMatrix col(const ZMatrix &z1, const ZMatrix &z2) {
	if(z1.line<=0 || z1.col<=0 || z1.mat==NULL || 
		z2.line<=0|| z2.col<=0 || z2.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(z2.line!=z1.line) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix z(z1.line, z1.col+z2.col);
	InsertSubMatrixIntoMatrix(z, z1,0,0);
	InsertSubMatrixIntoMatrix(z, z2,0,z1.col);
	return z;

}

DMatrix col(const DMatrix &z1, const DMatrix &z2) {
	if(z1.line<=0 || z1.col<=0 || z1.mat==NULL || 
		z2.line<=0|| z2.col<=0 || z2.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(z2.line!=z1.line) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	DMatrix z(z1.line, z1.col+z2.col);
	InsertSubMatrixIntoMatrix(z, z1,0,0);
	InsertSubMatrixIntoMatrix(z, z2,0,z1.col);
	return z;

}

WMatrix col(const WMatrix &z1, const WMatrix &z2) {
	if(z1.line<=0 || z1.col<=0 || z1.mat==NULL || 
		z2.line<=0|| z2.col<=0 || z2.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(z2.line!=z1.line) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	WMatrix z(z1.line, z1.col+z2.col);
	InsertSubMatrixIntoMatrix(z, z1,0,0);
	InsertSubMatrixIntoMatrix(z, z2,0,z1.col);
	return z;
}

ZMatrix line(const ZMatrix &z1, const ZMatrix &z2) {
	if(z1.line<=0 || z1.col<=0 || z1.mat==NULL || 
		z2.line<=0|| z2.col<=0 || z2.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(z2.col!=z1.col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix z(z1.line+z2.line, z1.col);
	InsertSubMatrixIntoMatrix(z, z1,0,0);
	InsertSubMatrixIntoMatrix(z, z2,z1.line, 0);
	return z;

}

DMatrix line(const DMatrix &z1, const DMatrix &z2) {
	if(z1.line<=0 || z1.col<=0 || z1.mat==NULL || 
		z2.line<=0|| z2.col<=0 || z2.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(z2.col!=z1.col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	DMatrix z(z1.line+z2.line, z1.col);
	InsertSubMatrixIntoMatrix(z, z1,0,0);
	InsertSubMatrixIntoMatrix(z, z2,z1.line, 0);
	return z;

}

WMatrix line(const WMatrix &z1, const WMatrix &z2) {
	if(z1.line<=0 || z1.col<=0 || z1.mat==NULL || 
		z2.line<=0|| z2.col<=0 || z2.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(z2.col!=z1.col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	WMatrix z(z1.line+z2.line, z1.col);
	InsertSubMatrixIntoMatrix(z, z1,0,0);
	InsertSubMatrixIntoMatrix(z, z2,z1.line, 0);
	return z;
}

ZMatrix col(const UVector<ZMatrix> &zv) {
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	if(zv.vect[0].mat==NULL || zv.vect[0].line<=0 || zv.vect[0].col<=0) 
		throw CExtUtilitisException(INVALID_MATRIX);
	ZMatrix z(zv.vect[0].line,zv.taille*zv.vect[0].col);
	for(int d=0;d<zv.taille;d++) {
		if(zv.vect[d].mat==NULL || zv.vect[d].line<=0|| zv.vect[d].col<=0)
			throw CExtUtilitisException(INVALID_MATRIX);
		if(zv.vect[d].line!=zv.vect[0].line|| zv.vect[d].col!=zv.vect[0].line)
			throw CExtUtilitisException(INCONSISTENT_MATRIX_SIZE_IN_UVECTOR);
		for(int l=0;l<zv.vect[d].line;l++) {
			for(int c=0;c<zv.vect[d].col;c++) {
				z.mat[l][d*zv.vect[d].col+c]=zv.vect[d].mat[l][c];
			}
		}
	}
	return z;
}

DMatrix col(const UVector<DMatrix> &zv) {
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	if(zv.vect[0].mat==NULL || zv.vect[0].line<=0 || zv.vect[0].col<=0) 
		throw CExtUtilitisException(INVALID_MATRIX);
	DMatrix z(zv.vect[0].line,zv.taille*zv.vect[0].col);
	for(int d=0;d<zv.taille;d++) {
		if(zv.vect[d].mat==NULL || zv.vect[d].line<=0|| zv.vect[d].col<=0)
			throw CExtUtilitisException(INVALID_MATRIX);
		if(zv.vect[d].line!=zv.vect[0].line|| zv.vect[d].col!=zv.vect[0].line)
			throw CExtUtilitisException(INCONSISTENT_MATRIX_SIZE_IN_UVECTOR);
		for(int l=0;l<zv.vect[d].line;l++) {
			for(int c=0;c<zv.vect[d].col;c++) {
				z.mat[l][d*zv.vect[d].col+c]=zv.vect[d].mat[l][c];
			}
		}
	}
	return z;
}

WMatrix col(const UVector<WMatrix> &zv) {
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	if(zv.vect[0].mat==NULL || zv.vect[0].line<=0 || zv.vect[0].col<=0) 
		throw CExtUtilitisException(INVALID_MATRIX);
	WMatrix z(zv.vect[0].line,zv.taille*zv.vect[0].col);
	for(int d=0;d<zv.taille;d++) {
		if(zv.vect[d].mat==NULL || zv.vect[d].line<=0|| zv.vect[d].col<=0)
			throw CExtUtilitisException(INVALID_MATRIX);
		if(zv.vect[d].line!=zv.vect[0].line|| zv.vect[d].col!=zv.vect[0].line)
			throw CExtUtilitisException(INCONSISTENT_MATRIX_SIZE_IN_UVECTOR);
		for(int l=0;l<zv.vect[d].line;l++) {
			for(int c=0;c<zv.vect[d].col;c++) {
				z.mat[l][d*zv.vect[d].col+c]=zv.vect[d].mat[l][c];
			}
		}
	}
	return z;
}

ZMatrix line(const UVector<ZMatrix> &zv) {
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	if(zv.vect[0].mat==NULL || zv.vect[0].line<=0 || zv.vect[0].col<=0) 
		throw CExtUtilitisException(INVALID_MATRIX);
	ZMatrix z(zv.taille*zv.vect[0].line, zv.vect[0].col);
	for(int d=0;d<zv.taille;d++) {
		if(zv.vect[d].mat==NULL || zv.vect[d].line<=0|| zv.vect[d].col<=0)
			throw CExtUtilitisException(INVALID_MATRIX);
		if(zv.vect[d].line!=zv.vect[0].line|| zv.vect[d].col!=zv.vect[0].line)
			throw CExtUtilitisException(INCONSISTENT_MATRIX_SIZE_IN_UVECTOR);
		for(int l=0;l<zv.vect[d].line;l++) {
			for(int c=0;c<zv.vect[d].col;c++) {
				z.mat[d*zv.vect[d].line+l][c]=zv.vect[d].mat[l][c];
			}
		}
	}
	return z;
}

DMatrix line(const UVector<DMatrix> &zv) {
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	if(zv.vect[0].mat==NULL || zv.vect[0].line<=0 || zv.vect[0].col<=0) 
		throw CExtUtilitisException(INVALID_MATRIX);
	DMatrix z(zv.taille*zv.vect[0].line, zv.vect[0].col);
	for(int d=0;d<zv.taille;d++) {
		if(zv.vect[d].mat==NULL || zv.vect[d].line<=0|| zv.vect[d].col<=0)
			throw CExtUtilitisException(INVALID_MATRIX);
		if(zv.vect[d].line!=zv.vect[0].line|| zv.vect[d].col!=zv.vect[0].line)
			throw CExtUtilitisException(INCONSISTENT_MATRIX_SIZE_IN_UVECTOR);
		for(int l=0;l<zv.vect[d].line;l++) {
			for(int c=0;c<zv.vect[d].col;c++) {
				z.mat[d*zv.vect[d].line+l][c]=zv.vect[d].mat[l][c];
			}
		}
	}
	return z;
}

WMatrix line(const UVector<WMatrix> &zv) {
	if(zv.vect==NULL || zv.taille<=0) throw CExtUtilitisException(INVALID_VECTOR);
	if(zv.vect[0].mat==NULL || zv.vect[0].line<=0 || zv.vect[0].col<=0) 
		throw CExtUtilitisException(INVALID_MATRIX);
	WMatrix z(zv.taille*zv.vect[0].line, zv.vect[0].col);
	for(int d=0;d<zv.taille;d++) {
		if(zv.vect[d].mat==NULL || zv.vect[d].line<=0|| zv.vect[d].col<=0)
			throw CExtUtilitisException(INVALID_MATRIX);
		if(zv.vect[d].line!=zv.vect[0].line|| zv.vect[d].col!=zv.vect[0].line)
			throw CExtUtilitisException(INCONSISTENT_MATRIX_SIZE_IN_UVECTOR);
		for(int l=0;l<zv.vect[d].line;l++) {
			for(int c=0;c<zv.vect[d].col;c++) {
				z.mat[d*zv.vect[d].line+l][c]=zv.vect[d].mat[l][c];
			}
		}
	}
	return z;
}



ZMatrix SubMatrix(const ZMatrix & m, int line, int col, int nb_line, int nb_col) {
	if(m.line<=0 || m.col <=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.line<line || m.col < col) throw CUtilitisException(__FILE__, __LINE__, MATRIX_INDEX_OUT_OF_RANGE);
	if(m.line<nb_line || m.col < nb_col) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	if(m.line<line+nb_line || m.col < col+nb_col) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	ZMatrix z(nb_line, nb_col);
	for(int i=line, ii=0; ii <nb_line; ii++, i++) {
		for(int j=col, jj=0; jj <nb_col; jj++, j++) {
			z.mat[ii][jj]=m.mat[i][j];
		}
	}
	return z;
}

DMatrix SubMatrix(const DMatrix & m, int line, int col, int nb_line, int nb_col) {
	if(m.line<=0 || m.col <=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.line<line || m.col < col) throw CUtilitisException(__FILE__, __LINE__, MATRIX_INDEX_OUT_OF_RANGE);
	if(m.line<nb_line || m.col < nb_col) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	if(m.line<line+nb_line || m.col < col+nb_col) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	DMatrix z(nb_line, nb_col);
	for(register int i=line, ii=0; ii <nb_line; ii++, i++) {
		for(register int j=col, jj=0; jj <nb_col; jj++, j++) {
			z.mat[ii][jj]=m.mat[i][j];
		}
	}
	return z;
}

WMatrix SubMatrix(const WMatrix & m, int line, int col, int nb_line, int nb_col) {
	if(m.line<=0 || m.col <=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.line<line || m.col < col) throw CUtilitisException(__FILE__, __LINE__, MATRIX_INDEX_OUT_OF_RANGE);
	if(m.line<nb_line || m.col < nb_col) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	if(m.line<line+nb_line || m.col < col+nb_col) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	WMatrix z(nb_line, nb_col);
	for(register int i=line, ii=0; ii <nb_line; ii++, i++) {
		for(register int j=col, jj=0; jj <nb_col; jj++, j++) {
			z.mat[ii][jj]=m.mat[i][j];
		}
	}
	return z;
}

void InsertSubMatrixIntoMatrix(const ZMatrix &dest, const ZMatrix &submatrix, int line_pos, int col_pos) {
	if(submatrix.line<=0 || submatrix.col <=0 || submatrix.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(dest.line<=0 || dest.col <=0 || dest.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);

	if(dest.line < submatrix.line || dest.col < submatrix.col) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	if(dest.line < line_pos || dest.col < col_pos) throw CUtilitisException(__FILE__, __LINE__, MATRIX_INDEX_OUT_OF_RANGE);
	if(dest.line < submatrix.line+line_pos || dest.col < submatrix.col+col_pos) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	
	for(register int i=line_pos, ii=0; ii <submatrix.line; ii++, i++) {
		for(register int j=col_pos, jj=0; jj <submatrix.col; jj++, j++) {
			dest.mat[i][j]=submatrix.mat[ii][jj];
		}
	}
}

void InsertSubMatrixIntoMatrix(const DMatrix &dest, const DMatrix &submatrix, int line_pos, int col_pos) {
	if(submatrix.line<=0 || submatrix.col <=0 || submatrix.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(dest.line<=0 || dest.col <=0 || dest.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);

	if(dest.line < submatrix.line || dest.col < submatrix.col) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	if(dest.line < line_pos || dest.col < col_pos) throw CUtilitisException(__FILE__, __LINE__, MATRIX_INDEX_OUT_OF_RANGE);
	if(dest.line < submatrix.line+line_pos || dest.col < submatrix.col+col_pos) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	
	for(register int i=line_pos, ii=0; ii <submatrix.line; ii++, i++) {
		for(register int j=col_pos, jj=0; jj <submatrix.col; jj++, j++) {
			dest.mat[i][j]=submatrix.mat[ii][jj];
		}
	}
}

void InsertSubMatrixIntoMatrix(const WMatrix &dest, const WMatrix &submatrix, int line_pos, int col_pos) {
	if(submatrix.line<=0 || submatrix.col <=0 || submatrix.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(dest.line<=0 || dest.col <=0 || dest.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);

	if(dest.line < submatrix.line || dest.col < submatrix.col) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	if(dest.line < line_pos || dest.col < col_pos) throw CUtilitisException(__FILE__, __LINE__, MATRIX_INDEX_OUT_OF_RANGE);
	if(dest.line < submatrix.line+line_pos || dest.col < submatrix.col+col_pos) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	
	for(register int i=line_pos, ii=0; ii <submatrix.line; ii++, i++) {
		for(register int j=col_pos, jj=0; jj <submatrix.col; jj++, j++) {
			dest.mat[i][j]=submatrix.mat[ii][jj];
		}
	}
}

ZMatrix ExtractSubMatrixFromMatrix(const ZMatrix &z, int line_pos, int col_pos, int nb_line, int nb_col) {
	if(z.line<=0 || z.col <=0 || z.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(z.line < nb_line || z.col < nb_col) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	if(z.line < line_pos || z.col < col_pos) throw CUtilitisException(__FILE__, __LINE__, MATRIX_INDEX_OUT_OF_RANGE);
	if(z.line < nb_line+line_pos || z.col < nb_col+col_pos) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	ZMatrix m(nb_line, nb_col);
	for(register int i=line_pos, ii=0;ii<nb_line;ii++,i++) {
		for(register int j=col_pos,jj=0;jj<nb_col;jj++,j++) {
			m.mat[ii][jj]=z.mat[i][j];
		}
	}
	return m;
}


DMatrix ExtractSubMatrixFromMatrix(const DMatrix &z, int line_pos, int col_pos, int nb_line, int nb_col) {
	if(z.line<=0 || z.col <=0 || z.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(z.line < nb_line || z.col < nb_col) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	if(z.line < line_pos || z.col < col_pos) throw CUtilitisException(__FILE__, __LINE__, MATRIX_INDEX_OUT_OF_RANGE);
	if(z.line < nb_line+line_pos || z.col < nb_col+col_pos) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	DMatrix m(nb_line, nb_col);
	for(register int i=line_pos, ii=0;ii<nb_line;ii++,i++) {
		for(register int j=col_pos,jj=0;jj<nb_col;jj++,j++) {
			m.mat[ii][jj]=z.mat[i][j];
		}
	}
	return m;
}

WMatrix ExtractSubMatrixFromMatrix(const WMatrix &z, int line_pos, int col_pos, int nb_line, int nb_col) {
	if(z.line<=0 || z.col <=0 || z.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(z.line < nb_line || z.col < nb_col) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	if(z.line < line_pos || z.col < col_pos) throw CUtilitisException(__FILE__, __LINE__, MATRIX_INDEX_OUT_OF_RANGE);
	if(z.line < nb_line+line_pos || z.col < nb_col+col_pos) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
	WMatrix m(nb_line, nb_col);
	for(register int i=line_pos, ii=0;ii<nb_line;ii++,i++) {
		for(register int j=col_pos,jj=0;jj<nb_col;jj++,j++) {
			m.mat[ii][jj]=z.mat[i][j];
		}
	}
	return m;
}

WMatrix ExtractLineXFromMatrix(const WMatrix &z, int line_number){
	if(z.line<=0 || z.col<=0 ||z.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(line_number >= z.line) throw CExtUtilitisException(MATRIX_INDEX_OUT_OF_RANGE);
	WMatrix n(1, z.col);
	for(register int i=0;i<z.col;i++) {
		n.mat[0][i]=z.mat[line_number][i];
	}
	return n;
}

DMatrix ExtractLineXFromMatrix(const DMatrix &z, int line_number){
	if(z.line<=0 || z.col<=0 ||z.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(line_number >= z.line) throw CExtUtilitisException(MATRIX_INDEX_OUT_OF_RANGE);
	DMatrix n(1, z.col);
	for(register int i=0;i<z.col;i++) {
		n.mat[0][i]=z.mat[line_number][i];
	}
	return n;
}

ZMatrix ExtractLineXFromMatrix(const ZMatrix &z, int line_number){
	if(z.line<=0 || z.col<=0 ||z.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(line_number >= z.line) throw CExtUtilitisException(MATRIX_INDEX_OUT_OF_RANGE);
	ZMatrix n(1, z.col);
	for(register int i=0;i<z.col;i++) {
		n.mat[0][i]=z.mat[line_number][i];
	}
	return n;
}

WMatrix ExtractColXFromMatrix(const WMatrix &z, int col_number){
	if(z.line<=0 || z.col<=0 ||z.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(col_number >= z.col) throw CExtUtilitisException(MATRIX_INDEX_OUT_OF_RANGE);
	WMatrix n(z.line, 1);
	for(register int i=0;i<z.line;i++) {
		n.mat[i][0]=z.mat[i][col_number];
	}
	return n;
}

DMatrix ExtractColXFromMatrix(const DMatrix &z, int col_number){
	if(z.line<=0 || z.col<=0 ||z.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(col_number >= z.col) throw CExtUtilitisException(MATRIX_INDEX_OUT_OF_RANGE);
	DMatrix n(z.line, 1);
	for(register int i=0;i<z.line;i++) {
		n.mat[i][0]=z.mat[i][col_number];
	}
	return n;
}

ZMatrix ExtractColXFromMatrix(const ZMatrix &z, int col_number){
	if(z.line<=0 || z.col<=0 ||z.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(col_number >= z.col) throw CExtUtilitisException(MATRIX_INDEX_OUT_OF_RANGE);
	ZMatrix n(z.line, 1);
	for(register int i=0;i<z.line;i++) {
		n.mat[i][0]=z.mat[i][col_number];
	}
	return n;
}

ofstream &matlaboutput(ofstream &os, const char* varname, const ZMatrix &m, int precision) {
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	int prevprec = os.precision();
	os.precision(precision);
	os << endl;
	os << varname << " = [" << endl;
	for(int i=0;i<m.line;i++) {
		os << "[";
		for(int j=0;j<m.col;j++) {
			os << m.mat[i][j] << " ";
		}
		os << "]" << endl;
	}
	os << "];" << endl;
	os.precision(prevprec);
	return os;
}


ofstream &matlaboutput(ofstream &os, const char* varname, const DMatrix &m, int precision) {
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	int prevprec = os.precision();
	os.precision(precision);
	os << endl;
	os << varname << " = [" << endl;
	for(int i=0;i<m.line;i++) {
		os << "[";
		for(int j=0;j<m.col;j++) {
			os << m.mat[i][j] << " ";
		}
		os << "]" << endl;
	}
	os << "];" << endl;
	os.precision(prevprec);
	return os;
}

ofstream &matlaboutput(ofstream &os, const char* varname, const WMatrix &m, int precision) {
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	int prevprec = os.precision();
	os.precision(precision);
	os << endl;
	os << varname << " = [" << endl;
	for(int i=0;i<m.line;i++) {
		os << "[";
		for(int j=0;j<m.col;j++) {
			os << m.mat[i][j] << " ";
		}
		os << "]" << endl;
	}
	os << "];" << endl;
	os.precision(prevprec);
	return os;
}

ofstream &matlaboutput(ofstream &os, const char* varname, const ZVector &m, int precision) {
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	int prevprec = os.precision();
	os.precision(precision);
	os << endl;
	os << varname << " = [";
	for(int i=0;i<m.taille;i++) {
			os << m.vect[i] << " ";
	}
	os << "];" << endl;
	os.precision(prevprec);
	return os;
}

ofstream &matlaboutput(ofstream &os, const char* varname, const DVector &m, int precision) {
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	int prevprec = os.precision();
	os.precision(precision);
	os << endl;
	os << varname << " = [";
	for(int i=0;i<m.taille;i++) {
			os << m.vect[i] << " ";
	}
	os << "];" << endl;
	os.precision(prevprec);
	return os;
}

ofstream &matlaboutput(ofstream &os, const char* varname, const WVector &m, int precision) {
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	int prevprec = os.precision();
	os.precision(precision);
	os << endl;
	os << varname << " = [";
	for(int i=0;i<m.taille;i++) {
			os << m.vect[i] << " ";
	}
	os << "];" << endl;
	os.precision(prevprec);
	return os;
}

ofstream &matlaboutput(ofstream &os, const char* varname, const DCplx &m, int precision) {
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	int prevprec = os.precision();
	os.precision(precision);
	os << endl;
	os << varname << " = " <<  m << endl;
	os.precision(prevprec);
	return os;
}

ofstream &matlaboutput(ofstream &os, const char* varname, double m, int precision) {
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	int prevprec = os.precision();
	os.precision(precision);
	os << endl;
	os << varname << " = " << m << ";" << endl;
	os.precision(prevprec);
	return os;
}


ZMatrix ConcatUVectorAsColumnMatrix(const UVector<ZVector> &u) {
	if(u.vect==NULL || u.taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_UVECTOR);
	if(u.vect[0].taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_UVECTOR);
	ZMatrix z(u.vect[0].taille*u.taille,1);
	for(register int i=0, j=0;i<u.taille;i++) {
		if(u.vect[i].taille!=u.vect[0].taille) throw CUtilitisException(__FILE__, __LINE__, INCONSISTENT_UVECTOR);
		for(register int k=0;k<u.vect[i].taille;k++) {
			z.mat[j++][0]=u.vect[i].vect[k];
		}
	}
	return z;
}

DMatrix ConcatUVectorAsColumnMatrix(const UVector<DVector> &u) {
	if(u.vect==NULL || u.taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_UVECTOR);
	if(u.vect[0].taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_UVECTOR);
	DMatrix z(u.vect[0].taille*u.taille,1);
	for(register int i=0, j=0;i<u.taille;i++) {
		if(u.vect[i].taille!=u.vect[0].taille) throw CUtilitisException(__FILE__, __LINE__, INCONSISTENT_UVECTOR);
		for(register int k=0;k<u.vect[i].taille;k++) {
			z.mat[j++][0]=u.vect[i].vect[k];
		}
	}
	return z;
}

WMatrix ConcatUVectorAsColumnMatrix(const UVector<WVector> &u) {
	if(u.vect==NULL || u.taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_UVECTOR);
	if(u.vect[0].taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_UVECTOR);
	WMatrix z(u.vect[0].taille*u.taille,1);
	for(register int i=0, j=0;i<u.taille;i++) {
		if(u.vect[i].taille!=u.vect[0].taille) throw CUtilitisException(__FILE__, __LINE__, INCONSISTENT_UVECTOR);
		for(register int k=0;k<u.vect[i].taille;k++) {
			z.mat[j++][0]=u.vect[i].vect[k];
		}
	}
	return z;
}

ZMatrix ConcatUVectorAsLineMatrix(const UVector<ZVector> &u) {
	if(u.vect==NULL || u.taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_UVECTOR);
	if(u.vect[0].taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_UVECTOR);
	ZMatrix z(1,u.vect[0].taille*u.taille);
	for(register int i=0, j=0;i<u.taille;i++) {
		if(u.vect[i].taille!=u.vect[0].taille) throw CUtilitisException(__FILE__, __LINE__, INCONSISTENT_UVECTOR);
		for(register int k=0;k<u.vect[i].taille;k++) {
			z.mat[0][j++]=u.vect[i].vect[k];
		}
	}
	return z;
}

DMatrix ConcatUVectorAsLineMatrix(const UVector<DVector> &u) {
	if(u.vect==NULL || u.taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_UVECTOR);
	if(u.vect[0].taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_UVECTOR);
	DMatrix z(1,u.vect[0].taille*u.taille);
	for(register int i=0, j=0;i<u.taille;i++) {
		if(u.vect[i].taille!=u.vect[0].taille) throw CUtilitisException(__FILE__, __LINE__, INCONSISTENT_UVECTOR);
		for(register int k=0;k<u.vect[i].taille;k++) {
			z.mat[0][j++]=u.vect[i].vect[k];
		}
	}
	return z;
}

WMatrix ConcatUVectorAsLineMatrix(const UVector<WVector> &u) {
	if(u.vect==NULL || u.taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_UVECTOR);
	if(u.vect[0].taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_UVECTOR);
	WMatrix z(1,u.vect[0].taille*u.taille);
	for(register int i=0, j=0;i<u.taille;i++) {
		if(u.vect[i].taille!=u.vect[0].taille) throw CUtilitisException(__FILE__, __LINE__, INCONSISTENT_UVECTOR);
		for(register int k=0;k<u.vect[i].taille;k++) {
			z.mat[0][j++]=u.vect[i].vect[k];
		}
	}
	return z;
}

DMatrix ReshapeColumnwise(const DVector &m, int new_line, int new_col) {
	if(m.taille<=0 || m.vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.taille!=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	DMatrix n(new_line, new_col);
	for(register int nc=0;nc<new_col;nc++) {
		for(register int nl=0;nl<new_line;nl++) {
			n.mat[nl][nc]=m.vect[nc*new_line+nl];
		}
	}
	return n;
}

WMatrix ReshapeColumnwise(const WVector &m, int new_line, int new_col) {
	if(m.taille<=0 || m.vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.taille!=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	WMatrix n(new_line, new_col);
	for(register int nc=0;nc<new_col;nc++) {
		for(register int nl=0;nl<new_line;nl++) {
			n.mat[nl][nc]=m.vect[nc*new_line+nl];
		}
	}
	return n;
}

ZMatrix ReshapeColumnwise(const ZVector &m, int new_line, int new_col) {
	if(m.taille<=0 || m.vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.taille!=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(new_line, new_col);
	for(register int nc=0;nc<new_col;nc++) {
		for(register int nl=0;nl<new_line;nl++) {
			n.mat[nl][nc]=m.vect[nc*new_line+nl];
		}
	}
	return n;
}

DMatrix ReshapeLinewise(const DVector &m, int new_line, int new_col) {
	if(m.taille<=0 || m.vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.taille!=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	DMatrix n(new_line, new_col);
	for(register int nl=0;nl<new_line;nl++) {
		for(register int nc=0;nc<new_col;nc++) {
			n.mat[nl][nc]=m.vect[nl*new_col+nc];
		}
	}
	return n;
}

WMatrix ReshapeLinewise(const WVector &m, int new_line, int new_col) {
	if(m.taille<=0 || m.vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.taille!=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	WMatrix n(new_line, new_col);
	for(register int nl=0;nl<new_line;nl++) {
		for(register int nc=0;nc<new_col;nc++) {
			n.mat[nl][nc]=m.vect[nl*new_col+nc];
		}
	}
	return n;
}

ZMatrix ReshapeLinewise(const ZVector &m, int new_line, int new_col) {
	if(m.taille<=0 || m.vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.taille!=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(new_line, new_col);
	for(register int nl=0;nl<new_line;nl++) {
		for(register int nc=0;nc<new_col;nc++) {
			n.mat[nl][nc]=m.vect[nl*new_col+nc];
		}
	}
	return n;
}

DMatrix ReshapeColumnwise(const DMatrix &m, int new_line, int new_col) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.line*m.col !=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	DMatrix n(new_line, new_col);
	register int nc=0,nl=0;
	for(register int mc=0;mc<m.col;mc++) {
		for(register int ml=0;ml<m.line;ml++) {
			n.mat[nl++][nc] =m.mat[ml][mc];
			if(nl>=new_line) {
				nl=0;nc++;
			}
		}
	}
	return n;
}

ZMatrix ReshapeColumnwise(const ZMatrix &m, int new_line, int new_col) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.line*m.col !=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(new_line, new_col);
	register int nc=0,nl=0;
	for(register int mc=0;mc<m.col;mc++) {
		for(register int ml=0;ml<m.line;ml++) {
			n.mat[nl++][nc] =m.mat[ml][mc];
			if(nl>=new_line) {
				nl=0;nc++;
			}
		}
	}
	return n;
}

WMatrix ReshapeColumnwise(const WMatrix &m, int new_line, int new_col) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.line*m.col !=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	WMatrix n(new_line, new_col);
	register int nc=0,nl=0;
	for(register int mc=0;mc<m.col;mc++) {
		for(register int ml=0;ml<m.line;ml++) {
			n.mat[nl++][nc] =m.mat[ml][mc];
			if(nl>=new_line) {
				nl=0;nc++;
			}
		}
	}
	return n;
}

DMatrix ReshapeLinewise(const DMatrix &m, int new_line, int new_col) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.line*m.col !=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	DMatrix n(new_line, new_col);
	register int nc=0,nl=0;
	for(register int ml=0;ml<m.line;ml++) {
		for(register int mc=0;mc<m.col;mc++) {
			n.mat[nl][nc++] =m.mat[ml][mc];
			if(nc>=new_col) {
				nc=0;nl++;
			}
		}
	}
	return n;
}

ZMatrix ReshapeLinewise(const ZMatrix &m, int new_line, int new_col) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.line*m.col !=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(new_line, new_col);
	register int nc=0,nl=0;
	for(register int ml=0;ml<m.line;ml++) {
		for(register int mc=0;mc<m.col;mc++) {
			n.mat[nl][nc++] =m.mat[ml][mc];
			if(nc>=new_col) {
				nc=0;nl++;
			}
		}
	}
	return n;
}

WMatrix ReshapeLinewise(const WMatrix &m, int new_line, int new_col) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.line*m.col !=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	WMatrix n(new_line, new_col);
	register int nc=0,nl=0;
	for(register int ml=0;ml<m.line;ml++) {
		for(register int mc=0;mc<m.col;mc++) {
			n.mat[nl][nc++] =m.mat[ml][mc];
			if(nc>=new_col) {
				nc=0;nl++;
			}
		}
	}
	return n;
}

DMatrix ReshapeLinewiseTranspose(const DMatrix &m, int new_line, int new_col) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.line*m.col !=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	DMatrix n(new_line, new_col);
	register int nc=0,nl=0;
	for(register int ml=0;ml<m.line;ml++) {
		for(register int mc=0;ml<m.col;mc++) {
			n.mat[nl][nc++] =m.mat[ml][mc];
			if(nc>=new_col) {
				nc=0;nl++;
			}
		}
	}
	return n;
}

ZMatrix ReshapeLinewiseTranspose(const ZMatrix &m, int new_line, int new_col) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.line*m.col !=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(new_line, new_col);
	register int nc=0,nl=0;
	for(register int ml=0;ml<m.line;ml++) {
		for(register int mc=0;ml<m.col;mc++) {
			n.mat[nl][nc++] =m.mat[ml][mc];
			if(nc>=new_col) {
				nc=0;nl++;
			}
		}
	}
	return n;
}


WMatrix ReshapeLinewiseTranspose(const WMatrix &m, int new_line, int new_col) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.line*m.col !=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	WMatrix n(new_line, new_col);
	register int nc=0,nl=0;
	for(register int ml=0;ml<m.line;ml++) {
		for(register int mc=0;ml<m.col;mc++) {
			n.mat[nl][nc++] =m.mat[ml][mc];
			if(nc>=new_col) {
				nc=0;nl++;
			}
		}
	}
	return n;
}

DMatrix ReshapeColumnwiseTranspose(const DMatrix &m, int new_line, int new_col) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.line*m.col !=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	DMatrix n(new_line, new_col);
	register int nc=0,nl=0;
	for(register int mc=0;mc<m.col;mc++) {
		for(register int ml=0;ml<m.line;ml++) {
			n.mat[nl++][nc] =m.mat[ml][mc];
			if(nl>=new_line) {
				nl=0;nc++;
			}
		}
	}
	return n;
}

ZMatrix ReshapeColumnwiseTranspose(const ZMatrix &m, int new_line, int new_col) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.line*m.col !=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(new_line, new_col);
	register int nc=0,nl=0;
	for(register int mc=0;mc<m.col;mc++) {
		for(register int ml=0;ml<m.line;ml++) {
			n.mat[nl++][nc] =m.mat[ml][mc];
			if(nl>=new_line) {
				nl=0;nc++;
			}
		}
	}
	return n;
}

WMatrix ReshapeColumnwiseTranspose(const WMatrix &m, int new_line, int new_col) {
	if(m.line<=0 || m.col<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(m.line*m.col !=new_line*new_col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	WMatrix n(new_line, new_col);
	register int nc=0,nl=0;
	for(register int mc=0;mc<m.col;mc++) {
		for(register int ml=0;ml<m.line;ml++) {
			n.mat[nl++][nc] =m.mat[ml][mc];
			if(nl>=new_line) {
				nl=0;nc++;
			}
		}
	}
	return n;
}

void RepeatMatrixDiagonally(const ZMatrix &src, ZMatrix &dest, int repeat_number, DCplx fill_with) {
	if(src.line<=0 || src.col<=0 || src.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(repeat_number<=0) throw CUtilitisException(__FILE__, __LINE__, INVALID_ARGUMENT);
	int new_line=repeat_number*src.line, new_col=repeat_number*src.col;
	dest = ZMatrix(new_line, new_col, fill_with);
	for(int c=0;c<repeat_number;c++) {
		register int line_offset=c*src.line, col_offset=c*src.col;
		for(register int i=0;i<src.line;i++) {
			for(register int j=0;j<src.col;j++) {
				dest.mat[i+line_offset][j+col_offset]=src.mat[i][j];
			}
		}
	}
}

void RepeatMatrixDiagonally(const DMatrix &src, DMatrix &dest, int repeat_number, double fill_with) {
	if(src.line<=0 || src.col<=0 || src.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(repeat_number<=0) throw CUtilitisException(__FILE__, __LINE__, INVALID_ARGUMENT);
	int new_line=repeat_number*src.line, new_col=repeat_number*src.col;
	dest = DMatrix(new_line, new_col, fill_with);
	for(int c=0;c<repeat_number;c++) {
		register int line_offset=c*src.line, col_offset=c*src.col;
		for(register int i=0;i<src.line;i++) {
			for(register int j=0;j<src.col;j++) {
				dest.mat[i+line_offset][j+col_offset]=src.mat[i][j];
			}
		}
	}
}

void RepeatMatrixDiagonally(const WMatrix &src, WMatrix &dest, int repeat_number, int  fill_with) {
	if(src.line<=0 || src.col<=0 || src.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(repeat_number<=0) throw CUtilitisException(__FILE__, __LINE__, INVALID_ARGUMENT);
	int new_line=repeat_number*src.line, new_col=repeat_number*src.col;
	dest = WMatrix(new_line, new_col, fill_with);
	for(int c=0;c<repeat_number;c++) {
		register int line_offset=c*src.line, col_offset=c*src.col;
		for(register int i=0;i<src.line;i++) {
			for(register int j=0;j<src.col;j++) {
				dest.mat[i+line_offset][j+col_offset]=src.mat[i][j];
			}
		}
	}
}

int StreamVector(const WVector &src, UVector<WVector> &dest, int nbstream) {
	if(src.taille<=0 || src.vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(nbstream>src.taille) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	int ispadded= src.taille%nbstream;
	int streamsize=(int)ceil((double)src.taille/(double)nbstream);
	dest = UVector<WVector>(nbstream);
	for(int t=0;t<nbstream;t++) {
		dest.vect[t]=WVector(streamsize);
	}
	for(int i=0,j=0,k=0;i<src.taille;) {
		dest.vect[j++].vect[k] = src.vect[i++];
		if(j==nbstream) {
			j=0;k++;
		}
	}
	return ispadded;
}

int StreamVector(const DVector &src, UVector<DVector> &dest, int nbstream) {
	if(src.taille<=0 || src.vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(nbstream>src.taille) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	int ispadded= src.taille%nbstream;
	int streamsize=(int)ceil((double)src.taille/(double)nbstream);
	dest = UVector<DVector>(nbstream);
	for(int t=0;t<nbstream;t++) {
		dest.vect[t]=DVector(streamsize);
	}
	for(int i=0,j=0,k=0;i<src.taille;) {
		dest.vect[j++].vect[k] = src.vect[i++];
		if(j==nbstream) {
			j=0;k++;
		}
	}
	return ispadded;
}

int StreamVector(const ZVector &src, UVector<ZVector> &dest, int nbstream) {
	if(src.taille<=0 || src.vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(nbstream>src.taille) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	int ispadded= src.taille%nbstream;
	int streamsize=(int)ceil((double)src.taille/(double)nbstream);
	dest = UVector<ZVector>(nbstream);
	for(int t=0;t<nbstream;t++) {
		dest.vect[t]=ZVector(streamsize);
	}
	for(int i=0,j=0,k=0;i<src.taille;) {
		dest.vect[j++].vect[k] = src.vect[i++];
		if(j==nbstream) {
			j=0;k++;
		}
	}
	return ispadded;
}

void AddColXinMatrixWithY(const ZMatrix &m, int x, DCplx y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(x>m.col) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	for(int i=0;i<m.line;i++) {
		m.mat[i][x] += y;
	}
}

void AddColXinMatrixWithY(const DMatrix &m, int x, double y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(x>m.col) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	for(int i=0;i<m.line;i++) {
		m.mat[i][x] += y;
	}
}

void AddColXinMatrixWithY(const DMatrix &m, int x, int y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(x>m.col) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	if (y==0) return;
	for(int i=0;i<m.line;i++) {
		m.mat[i][x] += y;
	}
}

void MultColXinMatrixWithY(const ZMatrix &m, int x, DCplx y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(x>m.col) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	for(int i=0;i<m.line;i++) {
		m.mat[i][x] *= y;
	}
}

void MultColXinMatrixWithY(const DMatrix &m, int x, double y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(x>m.col) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	for(int i=0;i<m.line;i++) {
		m.mat[i][x] *= y;
	}
}

void MultColXinMatrixWithY(const DMatrix &m, int x, int y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(x>m.col) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	if (y==1) return;
	for(int i=0;i<m.line;i++) {
		m.mat[i][x] *= y;
	}
}


void AddLineXinMatrixWithY(const ZMatrix &m, int x, DCplx y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(x>m.line) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	for(int i=0;i<m.col;i++) {
		m.mat[x][i] += y;
	}
}

void AddLineXinMatrixWithY(const DMatrix &m, int x, double y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(x>m.line) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	for(int i=0;i<m.col;i++) {
		m.mat[x][i] += y;
	}
}

void AddLineXinMatrixWithY(const DMatrix &m, int x, int y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(x>m.line) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	if (y==0) return;
	for(int i=0;i<m.col;i++) {
		m.mat[x][i] += y;
	}
}

void MultLineXinMatrixWithY(const ZMatrix &m, int x, DCplx y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(x>m.line) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	for(int i=0;i<m.col;i++) {
		m.mat[i][x] *= y;
	}
}

void MultLineXinMatrixWithY(const DMatrix &m, int x, double y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(x>m.line) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	for(int i=0;i<m.col;i++) {
		m.mat[i][x] *= y;
	}
}

void MultLineXinMatrixWithY(const DMatrix &m, int x, int y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(x>m.line) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	if (y==1) return;
	for(int i=0;i<m.col;i++) {
		m.mat[i][x] *= y;
	}
}

void AddEachColinMatrixWithY(const ZMatrix &m, const ZVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.col) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.col; i++) {
		for(register int j=0;j<m.line;j++) {
			m.mat[j][i]+=y.vect[i];
		}
	}
}

void AddEachColinMatrixWithY(const ZMatrix &m, const WVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.col) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.col; i++) {
		if (y.vect[i]==0) continue;
		for(register int j=0;j<m.line;j++) {
			m.mat[j][i]+=y.vect[i];
		}
	}
}

void AddEachColinMatrixWithY(const ZMatrix &m, const DVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.col) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.col; i++) {
		for(register int j=0;j<m.line;j++) {
			m.mat[j][i]+=y.vect[i];
		}
	}
}

void AddEachColinMatrixWithY(const DMatrix &m, const DVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.col) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.col; i++) {
		for(register int j=0;j<m.line;j++) {
			m.mat[j][i]+=y.vect[i];
		}
	}
}

void AddEachColinMatrixWithY(const DMatrix &m, const WVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.col) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.col; i++) {
		if (y.vect[i]==0) continue;
		for(register int j=0;j<m.line;j++) {
			m.mat[j][i]+=y.vect[i];
		}
	}
}

void AddEachColinMatrixWithY(const WMatrix &m, const WVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.col) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.col; i++) {
		if (y.vect[i]==0) continue;
		for(register int j=0;j<m.line;j++) {
			m.mat[j][i]+=y.vect[i];
		}
	}
}

void AddEachLineinMatrixWithY(const ZMatrix &m, const ZVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.line) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.line; i++) {
		for(register int j=0;j<m.col;j++) {
			m.mat[i][j]+=y.vect[i];
		}
	}
}

void AddEachLineinMatrixWithY(const ZMatrix &m, const WVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.line) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.line; i++) {
		for(register int j=0;j<m.col;j++) {
			m.mat[i][j]+=y.vect[i];
		}
	}
}

void AddEachLineinMatrixWithY(const ZMatrix &m, const DVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.line) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.line; i++) {
		for(register int j=0;j<m.col;j++) {
			m.mat[i][j]+=y.vect[i];
		}
	}
}

void AddEachLineinMatrixWithY(const DMatrix &m, const DVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.line) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.line; i++) {
		for(register int j=0;j<m.col;j++) {
			m.mat[i][j]+=y.vect[i];
		}
	}
}

void AddEachLineinMatrixWithY(const DMatrix &m, const WVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.line) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.line; i++) {
		if (y.vect[i]==0) continue;
		for(register int j=0;j<m.col;j++) {
			m.mat[i][j]+=y.vect[i];
		}
	}
}

void AddEachLineinMatrixWithY(const WMatrix &m, const WVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.line) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.line; i++) {
		if (y.vect[i]==0) continue;
		for(register int j=0;j<m.col;j++) {
			m.mat[i][j]+=y.vect[i];
		}
	}
}

void MultEachColinMatrixWithY(const ZMatrix &m, const ZVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.col) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.col; i++) {
		for(register int j=0;j<m.line;j++) {
			m.mat[j][i]*=y.vect[i];
		}
	}
}

void MultEachColinMatrixWithY(const ZMatrix &m, const WVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.col) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.col; i++) {
		if (y.vect[i]==1) continue;
		for(register int j=0;j<m.line;j++) {
			m.mat[j][i]*=y.vect[i];
		}
	}
}

void MultEachColinMatrixWithY(const ZMatrix &m, const DVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.col) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.col; i++) {
		for(register int j=0;j<m.line;j++) {
			m.mat[j][i]*=y.vect[i];
		}
	}
}

void MultEachColinMatrixWithY(const DMatrix &m, const DVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.col) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.col; i++) {
		for(register int j=0;j<m.line;j++) {
			m.mat[j][i]*=y.vect[i];
		}
	}
}

void MultEachColinMatrixWithY(const DMatrix &m, const WVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.col) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.col; i++) {
		if (y.vect[i]==1) continue;
		for(register int j=0;j<m.line;j++) {
			m.mat[j][i]*=y.vect[i];
		}
	}
}

void MultEachColinMatrixWithY(const WMatrix &m, const WVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.col) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.col; i++) {
		if (y.vect[i]==1) continue;
		for(register int j=0;j<m.line;j++) {
			m.mat[j][i]*=y.vect[i];
		}
	}
}

void MultEachLineinMatrixWithY(const ZMatrix &m, const ZVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.line) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.line; i++) {
		for(register int j=0;j<m.col;j++) {
			m.mat[i][j]*=y.vect[i];
		}
	}
}

void MultEachLineinMatrixWithY(const ZMatrix &m, const WVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.line) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.line; i++) {
		if (y.vect[i]==1) continue;
		for(register int j=0;j<m.col;j++) {
			m.mat[i][j]*=y.vect[i];
		}
	}
}

void MultEachLineinMatrixWithY(const ZMatrix &m, const DVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.line) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.line; i++) {
		for(register int j=0;j<m.col;j++) {
			m.mat[i][j]*=y.vect[i];
		}
	}
}

void MultEachLineinMatrixWithY(const DMatrix &m, const DVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.line) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.line; i++) {
		for(register int j=0;j<m.col;j++) {
			m.mat[i][j]*=y.vect[i];
		}
	}
}

void MultEachLineinMatrixWithY(const DMatrix &m, const WVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.line) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.line; i++) {
		if (y.vect[i]==1) continue;
		for(register int j=0;j<m.col;j++) {
			m.mat[i][j]*=y.vect[i];
		}
	}
}

void MultEachLineinMatrixWithY(const WMatrix &m, const WVector y) {
	if(m.col<=0 || m.line<=0 || m.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(y.taille!=m.line) throw CUtilitisException(__FILE__, __LINE__, INVALID_DIMENSION);
	for(register int i=0; i<m.line; i++) {
		if (y.vect[i]==1) continue;
		for(register int j=0;j<m.col;j++) {
			m.mat[i][j]*=y.vect[i];
		}
	}
}

ZMatrix InterleaveTwoMatrixColumnWise(const ZMatrix &z1, const ZMatrix &z2) {
	if(z1.col<=0 || z1.line<=0 || z1.mat==NULL || z2.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(z1.col!=z2.col || z1.line!=z2.line ) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(z1.line, 2*z1.col);
	for(int i=0, j=0,k=0;i<z1.col;i++) {
		for(j=0;j<z1.line;j++) {
			n.mat[j][k]=z1.mat[j][i];
		}
		k++;
		for(j=0;j<z2.line;j++) {
			n.mat[j][k]=z2.mat[j][i];
		}
		k++;
	}
	return n;
}

DMatrix InterleaveTwoMatrixColumnWise(const DMatrix &z1, const DMatrix &z2) {
	if(z1.col<=0 || z1.line<=0 || z1.mat==NULL || z2.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(z1.col!=z2.col || z1.line!=z2.line ) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	DMatrix n(z1.line, 2*z1.col);
	for(int i=0,j=0,k=0;i<z1.col;i++) {
		for(j=0;j<z1.line;j++) {
			n.mat[j][k]=z1.mat[j][i];
		}
		k++;
		for(j=0;j<z2.line;j++) {
			n.mat[j][k]=z2.mat[j][i];
		}
		k++;
	}
	return n;
}

WMatrix InterleaveTwoMatrixColumnWise(const WMatrix &z1, const WMatrix &z2) {
	if(z1.col<=0 || z1.line<=0 || z1.mat==NULL || z2.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(z1.col!=z2.col || z1.line!=z2.line ) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	WMatrix n(z1.line, 2*z1.col);
	for(int i=0,j=0,k=0;i<z1.col;i++) {
		for(j=0;j<z1.line;j++) {
			n.mat[j][k]=z1.mat[j][i];
		}
		k++;
		for(j=0;j<z2.line;j++) {
			n.mat[j][k]=z2.mat[j][i];
		}
		k++;
	}
	return n;
}

ZMatrix InterleaveTwoMatrixLineWise(const ZMatrix &z1, const ZMatrix &z2) {
	if(z1.col<=0 || z1.line<=0 || z1.mat==NULL || z2.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(z1.col!=z2.col || z1.line!=z2.line ) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(2*z1.line, z1.col);
	for(int i=0,j=0,k=0;i<z1.line;i++) {
		for(j=0;j<z1.col;j++) {
			n.mat[k][j]=z1.mat[i][j];
		}
		k++;
		for(j=0;j<z2.col;j++) {
			n.mat[k][j]=z2.mat[i][j];
		}
		k++;
	}
	return n;
}

DMatrix InterleaveTwoMatrixLineWise(const DMatrix &z1, const DMatrix &z2) {
	if(z1.col<=0 || z1.line<=0 || z1.mat==NULL || z2.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(z1.col!=z2.col || z1.line!=z2.line ) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	DMatrix n(2*z1.line, z1.col);
	for(int i=0,j=0,k=0;i<z1.line;i++) {
		for(j=0;j<z1.col;j++) {
			n.mat[k][j]=z1.mat[i][j];
		}
		k++;
		for(j=0;j<z2.col;j++) {
			n.mat[k][j]=z2.mat[i][j];
		}
		k++;
	}
	return n;
}

WMatrix InterleaveTwoMatrixLineWise(const WMatrix &z1, const WMatrix &z2) {
	if(z1.col<=0 || z1.line<=0 || z1.mat==NULL || z2.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(z1.col!=z2.col || z1.line!=z2.line ) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	WMatrix n(2*z1.line, z1.col);
	for(int i=0,j=0,k=0;i<z1.line;i++) {
		for(j=0;j<z1.col;j++) {
			n.mat[k][j]=z1.mat[i][j];
		}
		k++;
		for(j=0;j<z2.col;j++) {
			n.mat[k][j]=z2.mat[i][j];
		}
		k++;
	}
	return n;
}

DMatrix exp(const DMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(d.line, d.col); 
	for (register int i=0;i<d.line;i++) {
		for(register int j=0;j<d.col;j++){
			n.mat[i][j]=exp(d.mat[i][j]);
		}
	}
	return n;
}

void LowerBound(const DMatrix &d, double lowerlimit) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	for (register int i=0;i<d.line;i++) {
		for(register int j=0;j<d.col;j++){
			d.mat[i][j]=d.mat[i][j]<lowerlimit ? lowerlimit : d.mat[i][j];
		}
	}

}

void HigherBound(const DMatrix &d, double higherlimit) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	for (register int i=0;i<d.line;i++) {
		for(register int j=0;j<d.col;j++){
			d.mat[i][j]=d.mat[i][j]>higherlimit ? higherlimit : d.mat[i][j];
		}
	}

}

void Bound(const DMatrix &d, double lowerlimit, double higherlimit) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	for (register int i=0;i<d.line;i++) {
		for(register int j=0;j<d.col;j++){
			d.mat[i][j]=d.mat[i][j]>higherlimit ? higherlimit : (d.mat[i][j]<lowerlimit ? lowerlimit : d.mat[i][j]);
		}
	}
}

DMatrix ColMax(const DMatrix& d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(2,d.col);
	double max;
	int maxpos;
	for(int j=0;j<d.col;j++) {
		max=d.mat[0][j];maxpos=0;
		for(int i=0;i<d.line;i++) {
			if(d.mat[i][j]>max) {
				max=d.mat[i][j]; maxpos=i;
			}
		}
		n.mat[0][j]=max;
		n.mat[1][j]=maxpos;
	}
	return n;
}

DMatrix ColMax(const DMatrix& d, WMatrix &pos) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(1,d.col);
	pos=WMatrix(1,d.col);
	double max;
	int maxpos;
	for(int j=0;j<d.col;j++) {
		max=d.mat[0][j];maxpos=0;
		for(int i=0;i<d.line;i++) {
			if(d.mat[i][j]>max) {
				max=d.mat[i][j]; maxpos=i;
			}
		}
		n.mat[0][j]=max;
		pos.mat[0][j]=maxpos;
	}
	return n;
}

DMatrix ColMin(const DMatrix& d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(2,d.col);
	double min;
	int minpos;
	for(int j=0;j<d.col;j++) {
		min=d.mat[0][j];minpos=0;
		for(int i=0;i<d.line;i++) {
			if(d.mat[i][j]<min) {
				min=d.mat[i][j]; minpos=i;
			}
		}
		n.mat[0][j]=min;
		n.mat[1][j]=minpos;
	}
	return n;
}

DMatrix ColMin(const DMatrix& d, WMatrix &pos) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(1,d.col);
	pos=WMatrix(1,d.col);
	double min;
	int minpos;
	for(int j=0;j<d.col;j++) {
		min=d.mat[0][j];minpos=0;
		for(int i=0;i<d.line;i++) {
			if(d.mat[i][j]<min) {
				min=d.mat[i][j]; minpos=i;
			}
		}
		n.mat[0][j]=min;
		pos.mat[0][j]=minpos;
	}
	return n;
}

DMatrix RowMax(const DMatrix& d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(d.line,2);
	double max;
	int maxpos;
	for(int i=0;i<d.line;i++) {
		max=d.mat[i][0];maxpos=0;
		for(int j=0;j<d.col;j++) {
			if(d.mat[i][j]>max) {
				max=d.mat[i][j]; maxpos=i;
			}
		}
		n.mat[i][0]=max;
		n.mat[i][1]=maxpos;
	}
	return n;
}

DMatrix RowMax(const DMatrix& d, WMatrix &pos) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(d.line,2);
	pos=WMatrix(d.line,2);
	double max;
	int maxpos;
	for(int i=0;i<d.line;i++) {
		max=d.mat[i][0];maxpos=0;
		for(int j=0;j<d.col;j++) {
			if(d.mat[i][j]>max) {
				max=d.mat[i][j]; maxpos=i;
			}
		}
		n.mat[i][0]=max;
		pos.mat[i][0]=maxpos;
	}
	return n;
}

DMatrix RowMin(const DMatrix& d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(d.line,2);
	double min;
	int minpos;
	for(int i=0;i<d.line;i++) {
		min=d.mat[i][0];minpos=0;
		for(int j=0;j<d.col;j++) {
			if(d.mat[i][j]>min) {
				min=d.mat[i][j]; minpos=i;
			}
		}
		n.mat[i][0]=min;
		n.mat[i][1]=minpos;
	}
	return n;
}

DMatrix RowMin(const DMatrix& d, WMatrix &pos) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(d.line,2);
	pos=WMatrix(d.line,2);
	double min;
	int minpos;
	for(int i=0;i<d.line;i++) {
		min=d.mat[i][0];minpos=0;
		for(int j=0;j<d.col;j++) {
			if(d.mat[i][j]>min) {
				min=d.mat[i][j]; minpos=i;
			}
		}
		n.mat[i][0]=min;
		pos.mat[i][0]=minpos;
	}
	return n;
}

WMatrix ColMax(const WMatrix& d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	WMatrix n(2,d.col);
	int max;
	int maxpos;
	for(int j=0;j<d.col;j++) {
		max=d.mat[0][j];maxpos=0;
		for(int i=0;i<d.line;i++) {
			if(d.mat[i][j]>max) {
				max=d.mat[i][j]; maxpos=i;
			}
		}
		n.mat[0][j]=max;
		n.mat[1][j]=maxpos;
	}
	return n;
}

WMatrix ColMax(const WMatrix& d, WMatrix &pos) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	WMatrix n(1,d.col);
	pos=WMatrix(1,d.col);
	int max;
	int maxpos;
	for(int j=0;j<d.col;j++) {
		max=d.mat[0][j];maxpos=0;
		for(int i=0;i<d.line;i++) {
			if(d.mat[i][j]>max) {
				max=d.mat[i][j]; maxpos=i;
			}
		}
		n.mat[0][j]=max;
		pos.mat[0][j]=maxpos;
	}
	return n;
}

WMatrix ColMin(const WMatrix& d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	WMatrix n(2,d.col);
	int min;
	int minpos;
	for(int j=0;j<d.col;j++) {
		min=d.mat[0][j];minpos=0;
		for(int i=0;i<d.line;i++) {
			if(d.mat[i][j]<min) {
				min=d.mat[i][j]; minpos=i;
			}
		}
		n.mat[0][j]=min;
		n.mat[1][j]=minpos;
	}
	return n;
}

WMatrix ColMin(const WMatrix& d, WMatrix &pos) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	WMatrix n(1,d.col);
	pos=WMatrix(1,d.col);
	int min;
	int minpos;
	for(int j=0;j<d.col;j++) {
		min=d.mat[0][j];minpos=0;
		for(int i=0;i<d.line;i++) {
			if(d.mat[i][j]<min) {
				min=d.mat[i][j]; minpos=i;
			}
		}
		n.mat[0][j]=min;
		pos.mat[0][j]=minpos;
	}
	return n;
}

WMatrix RowMax(const WMatrix& d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	WMatrix n(d.line,2);
	int max;
	int maxpos;
	for(int i=0;i<d.line;i++) {
		max=d.mat[i][0];maxpos=0;
		for(int j=0;j<d.col;j++) {
			if(d.mat[i][j]>max) {
				max=d.mat[i][j]; maxpos=i;
			}
		}
		n.mat[i][0]=max;
		n.mat[i][1]=maxpos;
	}
	return n;
}

WMatrix RowMax(const WMatrix& d, WMatrix &pos) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	WMatrix n(d.line,2);
	pos=WMatrix(d.line,2);
	int max;
	int maxpos;
	for(int i=0;i<d.line;i++) {
		max=d.mat[i][0];maxpos=0;
		for(int j=0;j<d.col;j++) {
			if(d.mat[i][j]>max) {
				max=d.mat[i][j]; maxpos=i;
			}
		}
		n.mat[i][0]=max;
		pos.mat[i][0]=maxpos;
	}
	return n;
}

WMatrix RowMin(const WMatrix& d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	WMatrix n(d.line,2);
	int min;
	int minpos;
	for(int i=0;i<d.line;i++) {
		min=d.mat[i][0];minpos=0;
		for(int j=0;j<d.col;j++) {
			if(d.mat[i][j]>min) {
				min=d.mat[i][j]; minpos=i;
			}
		}
		n.mat[i][0]=min;
		n.mat[i][1]=minpos;
	}
	return n;
}

WMatrix RowMin(const WMatrix& d, WMatrix &pos) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	WMatrix n(d.line,2);
	pos=WMatrix(d.line,2);
	int min;
	int minpos;
	for(int i=0;i<d.line;i++) {
		min=d.mat[i][0];minpos=0;
		for(int j=0;j<d.col;j++) {
			if(d.mat[i][j]>min) {
				min=d.mat[i][j]; minpos=i;
			}
		}
		n.mat[i][0]=min;
		pos.mat[i][0]=minpos;
	}
	return n;
}

DMatrix ColSum(const DMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(1,d.col,0.0);
	for(register int j=0;j<d.col;j++) {
		for(register int i=0;i<d.line;i++) {
			n.mat[0][j]+=d.mat[i][j];
		}
	}
	return n;
}

DMatrix RowSum(const DMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(d.line,1,0.0);
	for(register int i=0;i<d.line;i++) {
		for(register int j=0;j<d.col;j++) {
			n.mat[i][0]+=d.mat[i][j];
		}
	}
	return n;
}

DMatrix ColAvg(const DMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(1,d.col,0.0);
	double coef=1.0/d.line;
	for(register int j=0;j<d.col;j++) {
		for(register int i=0;i<d.line;i++) {
			n.mat[0][j]+=d.mat[i][j];
		}
	}
	return n*coef;
}

DMatrix RowAvg(const DMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(d.line,1,0.0);
	double coef=1.0/d.col;
	for(register int i=0;i<d.line;i++) {
		for(register int j=0;j<d.col;j++) {
			n.mat[i][0]+=d.mat[i][j];
		}
	}
	return n*coef;
}


WMatrix ColSum(const WMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	WMatrix n(1,d.col,0);
	for(register int j=0;j<d.col;j++) {
		for(register int i=0;i<d.line;i++) {
			n.mat[0][j]+=d.mat[i][j];
		}
	}
	return n;
}

WMatrix RowSum(const WMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	WMatrix n(d.line,1,0);
	for(register int i=0;i<d.line;i++) {
		for(register int j=0;j<d.col;j++) {
			n.mat[i][0]+=d.mat[i][j];
		}
	}
	return n;
}

DMatrix ColAvg(const WMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(1,d.col,0.0);
	double coef=1.0/d.line;
	for(register int j=0;j<d.col;j++) {
		for(register int i=0;i<d.line;i++) {
			n.mat[0][j]+=d.mat[i][j];
		}
	}
	return n*coef;
}

DMatrix RowAvg(const WMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(d.line,1,0.0);
	double coef=1.0/d.col;
	for(register int i=0;i<d.line;i++) {
		for(register int j=0;j<d.col;j++) {
			n.mat[i][0]+=d.mat[i][j];
		}
	}
	return n*coef;
}


WVector ExtractDiagonal(const WMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(d.col!=d.line) throw CUtilitisException(__FILE__, __LINE__, MATRIX_NOT_SQUARE);
	WVector diag(d.line);
	for(register int i=0;i<d.line;i++) {
		diag.vect[i]=d.mat[i][i];
	}
	return diag;
}

DVector ExtractDiagonal(const DMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(d.col!=d.line) throw CUtilitisException(__FILE__, __LINE__, MATRIX_NOT_SQUARE);
	DVector diag(d.line);
	for(register int i=0;i<d.line;i++) {
		diag.vect[i]=d.mat[i][i];
	}
	return diag;
}

ZVector ExtractDiagonal(const ZMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(d.col!=d.line) throw CUtilitisException(__FILE__, __LINE__, MATRIX_NOT_SQUARE);
	ZVector diag(d.line);
	for(register int i=0;i<d.line;i++) {
		diag.vect[i]=d.mat[i][i];
	}
	return diag;
}


WMatrix ExtractDiagonalAsColumnMatrix(const WMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(d.col!=d.line) throw CUtilitisException(__FILE__, __LINE__, MATRIX_NOT_SQUARE);
	WMatrix diag(d.line, 1);
	for(register int i=0;i<d.line;i++) {
		diag.mat[i][0]=d.mat[i][i];
	}
	return diag;
}

DMatrix ExtractDiagonalAsColumnMatrix(const DMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(d.col!=d.line) throw CUtilitisException(__FILE__, __LINE__, MATRIX_NOT_SQUARE);
	DMatrix diag(d.line, 1);
	for(register int i=0;i<d.line;i++) {
		diag.mat[i][0]=d.mat[i][i];
	}
	return diag;
}

ZMatrix ExtractDiagonalAsColumnMatrix(const ZMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(d.col!=d.line) throw CUtilitisException(__FILE__, __LINE__, MATRIX_NOT_SQUARE);
	ZMatrix diag(d.line, 1);
	for(register int i=0;i<d.line;i++) {
		diag.mat[i][0]=d.mat[i][i];
	}
	return diag;
}

WMatrix ExtractDiagonalAsLineMatrix(const WMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(d.col!=d.line) throw CUtilitisException(__FILE__, __LINE__, MATRIX_NOT_SQUARE);
	WMatrix diag(1, d.line);
	for(register int i=0;i<d.line;i++) {
		diag.mat[0][i]=d.mat[i][i];
	}
	return diag;
}

DMatrix ExtractDiagonalAsLineMatrix(const DMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(d.col!=d.line) throw CUtilitisException(__FILE__, __LINE__, MATRIX_NOT_SQUARE);
	DMatrix diag(1, d.line);
	for(register int i=0;i<d.line;i++) {
		diag.mat[0][i]=d.mat[i][i];
	}
	return diag;
}

ZMatrix ExtractDiagonalAsLineMatrix(const ZMatrix &d) {
	if(d.col<=0 || d.line<=0 || d.mat==NULL || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(d.col!=d.line) throw CUtilitisException(__FILE__, __LINE__, MATRIX_NOT_SQUARE);
	ZMatrix diag(1, d.line);
	for(register int i=0;i<d.line;i++) {
		diag.mat[0][i]=d.mat[i][i];
	}
	return diag;
}


ZMatrix ColumnCirculant(const ZVector& col) {
	if(col.vect==NULL || col.taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	ZMatrix n(col.taille, col.taille);
	for(register int i=0; i < col.taille; i++) { // col
		for(register int j=0; j < col.taille; j++) { // line
			n.mat[(j+i)%col.taille][i]=col.vect[j];
		}
	}
	return n;
}

DMatrix ColumnCirculant(const DVector& col) {
	if(col.vect==NULL || col.taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	DMatrix n(col.taille, col.taille);
	for(register int i=0; i < col.taille; i++) { // col
		for(register int j=0; j < col.taille; j++) { // line
			n.mat[(j+i)%col.taille][i]=col.vect[j];
		}
	}
	return n;
}

WMatrix ColumnCirculant(const WVector& col) {
	if(col.vect==NULL || col.taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	WMatrix n(col.taille, col.taille);
	for(register int i=0; i < col.taille; i++) { // col
		for(register int j=0; j < col.taille; j++) { // line
			n.mat[(j+i)%col.taille][i]=col.vect[j];
		}
	}
	return n;
}

ZMatrix LineCirculant(const ZVector& col) {
	if(col.vect==NULL || col.taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	ZMatrix n(col.taille, col.taille);
	for(register int i=0; i < col.taille; i++) { // col
		for(register int j=0; j < col.taille; j++) { // line
			n.mat[i][(j+i)%col.taille]=col.vect[j];
		}
	}
	return n;
}

DMatrix LineCirculant(const DVector& col) {
	if(col.vect==NULL || col.taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	DMatrix n(col.taille, col.taille);
	for(register int i=0; i < col.taille; i++) { // col
		for(register int j=0; j < col.taille; j++) { // line
			n.mat[i][(j+i)%col.taille]=col.vect[j];
		}
	}
	return n;
}

WMatrix LineCirculant(const WVector& col) {
	if(col.vect==NULL || col.taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	WMatrix n(col.taille, col.taille);
	for(register int i=0; i < col.taille; i++) { // col
		for(register int j=0; j < col.taille; j++) { // line
			n.mat[i][(j+i)%col.taille]=col.vect[j];
		}
	}
	return n;
}

WMatrix DiagPackedVectorAsMatrix(const WVector& v, int matline) {
	if(v.vect==NULL || v.taille <= 0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(v.taille%matline) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	WMatrix m(matline, v.taille);
	register int j=0;
	for(register int i=0;i<v.taille;i++) {
		m.mat[j++][i]=v.vect[i];
		if(j++==matline) j=0;
	}
	return m;
}

DMatrix DiagPackedVectorAsMatrix(const DVector& v, int matline) {
	if(v.vect==NULL || v.taille <= 0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(v.taille%matline) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	DMatrix m(matline, v.taille);
	register int j=0;
	for(register int i=0;i<v.taille;i++) {
		m.mat[j++][i]=v.vect[i];
		if(j++==matline) j=0;
	}
	return m;
}

ZMatrix DiagPackedVectorAsMatrix(const ZVector& v, int matline) {
	if(v.vect==NULL || v.taille <= 0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(v.taille%matline) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	ZMatrix m(matline, v.taille);
	register int j=0;
	for(register int i=0;i<v.taille;i++) {
		m.mat[j++][i]=v.vect[i];
		if(j++==matline) j=0;
	}
	return m;
}

