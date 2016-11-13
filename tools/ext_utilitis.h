/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/06/06 09:19:09 $
 * $Revision: 1.1.2.19 $
 * $Id: ext_utilitis.h,v 1.1.2.19 2006/06/06 09:19:09 syed Exp $
 ********************************************************************/
#ifndef _TOOLS_EXT_UTILITIS_H_
#define _TOOLS_EXT_UTILITIS_H_
//#include <fstream.h>
#include <iostream>
#include <fstream>
using namespace std;
#include "tools/utilitis.h"
#include "tools/superclass.h"
#include "tools/ext_utilitis.h"
#include "globaldef.h"

class CExtUtilitisException{
public : 
	int err_mesg;
	CExtUtilitisException(int code) : err_mesg(code){}
};


UVector<ZVector> SplitMatrixByCol(const ZMatrix &m); // on aura col vecteurs de taille line

UVector<ZVector> SplitMatrixByLine(const ZMatrix &m);

UVector<DVector> SplitMatrixByCol(const DMatrix &m);

UVector<DVector> SplitMatrixByLine(const DMatrix &m);

UVector<WVector> SplitMatrixByCol(const WMatrix &m);

UVector<WVector> SplitMatrixByLine(const WMatrix &m);

ZMatrix JoinUVectorInCol(const UVector<ZVector> &u); // on aura une matrice de taille u.vect[0].taille x u.taille

ZMatrix JoinUVectorInLine(const UVector<ZVector> &u); // on aura une matrice de taille u.vect[0].taille x u.taille

DMatrix JoinUVectorInCol(const UVector<DVector> &u); // on aura une matrice de taille u.vect[0].taille x u.taille

DMatrix JoinUVectorInLine(const UVector<DVector> &u); // on aura une matrice de taille u.vect[0].taille x u.taille

WMatrix JoinUVectorInCol(const UVector<WVector> &u); // on aura une matrice de taille u.vect[0].taille x u.taille

WMatrix JoinUVectorInLine(const UVector<WVector> &u); // on aura une matrice de taille u.vect[0].taille x u.taille

ZMatrix VectorInColMatrixRepresentation(const ZVector &zv); // on aura une matrice de taille zv.taille x 1

ZMatrix VectorInLineMatrixRepresentation(const ZVector &zv); // on aura une matrice de taille 1 x zv.taille 

DMatrix VectorInColMatrixRepresentation(const DVector &zv); // on aura une matrice de taille zv.taille x 1

DMatrix VectorInLineMatrixRepresentation(const DVector &zv); // on aura une matrice de taille 1 x zv.taille 

WMatrix VectorInColMatrixRepresentation(const WVector &zv); // on aura une matrice de taille zv.taille x 1

WMatrix VectorInLineMatrixRepresentation(const WVector &zv); // on aura une matrice de taille 1 x zv.taille 

ZVector MatrixInVectorRepresentation(const ZMatrix &zm);

DVector MatrixInVectorRepresentation(const DMatrix &zm);

WVector MatrixInVectorRepresentation(const WMatrix &zm);

ZVector MatrixAsColumnPackedVector(const ZMatrix &z);

ZVector MatrixAsLinePackedVector(const ZMatrix &z);

DVector MatrixAsColumnPackedVector(const DMatrix &z);

DVector MatrixAsLinePackedVector(const DMatrix &z);

WVector MatrixAsColumnPackedVector(const WMatrix &z);

WVector MatrixAsLinePackedVector(const WMatrix &z);

WVector ColXInWMatrixToVector(const WMatrix &zm, int col); // col = [0;zm.col[

WVector LineXInWMatrixToVector(const WMatrix &zm, int line); // line = [0;zm.line[

DVector ColXInDMatrixToVector(const DMatrix &zm, int col); // col = [0;zm.col[

DVector LineXInDMatrixToVector(const DMatrix &zm, int line); // line = [0;zm.line[

ZVector ColXInZMatrixToVector(const ZMatrix &zm, int col); // col = [0;zm.col[

ZVector LineXInZMatrixToVector(const ZMatrix &zm, int line); // line = [0;zm.line[

ZMatrix ZVectorToDiagonalMatrix(const ZVector &zv);

DMatrix DVectorToDiagonalMatrix(const DVector &zv);

WMatrix WVectorToDiagonalMatrix(const WVector &zv);

ZMatrix LineXinMatrixAsDiagonal(const ZMatrix &m, int x); // x starts with 0

DMatrix LineXinMatrixAsDiagonal(const DMatrix &m, int x); // x starts with 0

WMatrix LineXinMatrixAsDiagonal(const WMatrix &m, int x); // x starts with 0

ZMatrix ColXinMatrixAsDiagonal(const ZMatrix &m, int x); // x starts with 0

DMatrix ColXinMatrixAsDiagonal(const DMatrix &m, int x); // x starts with 0

WMatrix ColXinMatrixAsDiagonal(const WMatrix &m, int x); // x starts with 0

ZMatrix diag(const ZVector &zv); // le vector est le composant d'une matrice diagonale

DMatrix diag(const DVector &zv);// le vector est le composant d'une matrice diagonale

WMatrix diag(const WVector &zv);// le vector est le composant d'une matrice diagonale

ZMatrix diag(const UVector<ZMatrix> &zv); // les matrices sont rangé diagonalement dans une plus grosse matrice

DMatrix diag(const UVector<DMatrix> &zv);// les matrices sont rangé diagonalement dans une plus grosse matrice

WMatrix diag(const UVector<WMatrix> &zv);// les matrices sont rangé diagonalement dans une plus grosse matrice

ZMatrix col(const ZMatrix &z1, const ZMatrix &z2); // concatener deux matrices horizontalement

DMatrix col(const DMatrix &z1, const DMatrix &z2); // concatener deux matrices horizontalement

WMatrix col(const WMatrix &z1, const WMatrix &z2); // concatener deux matrices horizontalement

ZMatrix line(const ZMatrix &z1, const ZMatrix &z2); // concatener deux matrices verticalement

DMatrix line(const DMatrix &z1, const DMatrix &z2);// concatener deux matrices verticalement

WMatrix line(const WMatrix &z1, const WMatrix &z2); // concatener deux matrices verticalement

ZMatrix col(const UVector<ZMatrix> &zv); // concatener plusieurs matrices horizontalement

DMatrix col(const UVector<DMatrix> &zv);// concatener plusieurs matrices horizontalement

WMatrix col(const UVector<WMatrix> &zv);// concatener plusieurs matrices horizontalement

ZMatrix line(const UVector<ZMatrix> &zv); // concatener plusieurs matrices verticalement

DMatrix line(const UVector<DMatrix> &zv);// concatener plusieurs matrices verticalement

WMatrix line(const UVector<WMatrix> &zv);// concatener plusieurs matrices verticalement

ofstream &matlaboutput(ofstream &os, const char* varname, double m, int precision=5);
ofstream &matlaboutput(ofstream &os, const char* varname, const ZMatrix &m, int precision=5);

ofstream &matlaboutput(ofstream &os, const char* varname, const DMatrix &m, int precision=5);

ofstream &matlaboutput(ofstream &os, const char* varname, const WMatrix &m, int precision=5);

ofstream &matlaboutput(ofstream &os, const char* varname, const ZVector &m, int precision=5);

ofstream &matlaboutput(ofstream &os, const char* varname, const DVector &m, int precision=5);

ofstream &matlaboutput(ofstream &os, const char* varname, const WVector &m, int precision=5);

ofstream &matlaboutput(ofstream &os, const char* varname, const DCplx &m, int precision=5);

ZMatrix SubMatrix(const ZMatrix & m, int line, int col, int nb_line, int nb_col); // extraire une sous matrice

DMatrix SubMatrix(const DMatrix & m, int line, int col, int nb_line, int nb_col);// extraire une sous matrice

WMatrix SubMatrix(const WMatrix & m, int line, int col, int nb_line, int nb_col);// extraire une sous matrice

void InsertSubMatrixIntoMatrix(const ZMatrix &dest, const ZMatrix &submatrix, int line_pos, int col_pos); // inserer une sous matrice

void InsertSubMatrixIntoMatrix(const DMatrix &dest, const DMatrix &submatrix, int line_pos, int col_pos); // inserer une sous matrice

void InsertSubMatrixIntoMatrix(const WMatrix &dest, const WMatrix &submatrix, int line_pos, int col_pos); // inserer une sous matrice

ZMatrix ExtractSubMatrixFromMatrix(const ZMatrix &z, int line_pos, int col_pos, int nb_line, int nb_col);// extraire une sous matrice

DMatrix ExtractSubMatrixFromMatrix(const DMatrix &z, int line_pos, int col_pos, int nb_line, int nb_col);// extraire une sous matrice

WMatrix ExtractSubMatrixFromMatrix(const WMatrix &z, int line_pos, int col_pos, int nb_line, int nb_col);// extraire une sous matrice

WMatrix ExtractLineXFromMatrix(const WMatrix &z, int line_number); // extraire ligne X de la matrice

DMatrix ExtractLineXFromMatrix(const DMatrix &z, int line_number);// extraire ligne X de la matrice

ZMatrix ExtractLineXFromMatrix(const ZMatrix &z, int line_number);// extraire ligne X de la matrice

WMatrix ExtractColXFromMatrix(const WMatrix &z, int col_number); // extraire colonne X de la matrice

DMatrix ExtractColXFromMatrix(const DMatrix &z, int col_number);// extraire colonne X de la matrice

ZMatrix ExtractColXFromMatrix(const ZMatrix &z, int col_number);// extraire colonne X de la matrice

ZMatrix ConcatUVectorAsColumnMatrix(const UVector<ZVector> &u);

DMatrix ConcatUVectorAsColumnMatrix(const UVector<DVector> &u);

WMatrix ConcatUVectorAsColumnMatrix(const UVector<WVector> &u);

ZMatrix ConcatUVectorAsLineMatrix(const UVector<ZVector> &u);

DMatrix ConcatUVectorAsLineMatrix(const UVector<DVector> &u);

WMatrix ConcatUVectorAsLineMatrix(const UVector<WVector> &u);

ZMatrix ColumnPackedVectorAsMatrix(const ZVector &z, int line, int col);

ZMatrix LinePackedVectorAsMatrix(const ZVector &z, int line, int col);

DMatrix ColumnPackedVectorAsMatrix(const DVector &z, int line, int col);

DMatrix LinePackedVectorAsMatrix(const DVector &z, int line, int col);

WMatrix ColumnPackedVectorAsMatrix(const WVector &z, int line, int col);

WMatrix LinePackedVectorAsMatrix(const WVector &z, int line, int col);

DMatrix ReshapeColumnwise(const DVector &m, int new_line, int new_col);

WMatrix ReshapeColumnwise(const WVector &m, int new_line, int new_col);

ZMatrix ReshapeColumnwise(const ZVector &m, int new_line, int new_col);

DMatrix ReshapeLinewise(const DVector &m, int new_line, int new_col);

WMatrix ReshapeLinewise(const WVector &m, int new_line, int new_col);

ZMatrix ReshapeLinewise(const ZVector &m, int new_line, int new_col);

DMatrix ReshapeColumnwise(const DMatrix &m, int new_line, int new_col);

ZMatrix ReshapeColumnwise(const ZMatrix &m, int new_line, int new_col);

WMatrix ReshapeColumnwise(const WMatrix &m, int new_line, int new_col);

DMatrix ReshapeLinewise(const DMatrix &m, int new_line, int new_col);

ZMatrix ReshapeLinewise(const ZMatrix &m, int new_line, int new_col);

WMatrix ReshapeLinewise(const WMatrix &m, int new_line, int new_col);

DMatrix ReshapeLinewiseTranspose(const DMatrix &m, int new_line, int new_col);

ZMatrix ReshapeLinewiseTranspose(const ZMatrix &m, int new_line, int new_col);

WMatrix ReshapeLinewiseTranspose(const WMatrix &m, int new_line, int new_col);

DMatrix ReshapeColumnwiseTranspose(const DMatrix &m, int new_line, int new_col);

ZMatrix ReshapeColumnwiseTranspose(const ZMatrix &m, int new_line, int new_col);

WMatrix ReshapeColumnwiseTranspose(const WMatrix &m, int new_line, int new_col);

void RepeatMatrixDiagonally(const ZMatrix &src, ZMatrix &dest, int repeat_number, DCplx fill_with);

void RepeatMatrixDiagonally(const DMatrix &src, DMatrix &dest, int repeat_number, double fill_with);

void RepeatMatrixDiagonally(const WMatrix &src, WMatrix &dest, int repeat_number, int  fill_with);

int StreamVector(const WVector &src, UVector<WVector> &dest, int nbstream); 

int StreamVector(const DVector &src, UVector<DVector> &dest, int nbstream); 

int StreamVector(const ZVector &src, UVector<ZVector> &dest, int nbstream); 

void AddColXinMatrixWithY(const ZMatrix &m, int x, DCplx y);

void AddColXinMatrixWithY(const DMatrix &m, int x, double y);

void AddColXinMatrixWithY(const DMatrix &m, int x, int y);

void MultColXinMatrixWithY(const ZMatrix &m, int x, DCplx y);

void MultColXinMatrixWithY(const DMatrix &m, int x, double y);

void MultColXinMatrixWithY(const DMatrix &m, int x, int y);

void AddLineXinMatrixWithY(const ZMatrix &m, int x, DCplx y);

void AddLineXinMatrixWithY(const DMatrix &m, int x, double y);

void AddLineXinMatrixWithY(const DMatrix &m, int x, int y);

void MultLineXinMatrixWithY(const ZMatrix &m, int x, DCplx y);

void MultLineXinMatrixWithY(const DMatrix &m, int x, double y);
	
void MultLineXinMatrixWithY(const DMatrix &m, int x, int y);

void AddEachColinMatrixWithY(const ZMatrix &m, const ZVector y);

void AddEachColinMatrixWithY(const ZMatrix &m, const WVector y);

void AddEachColinMatrixWithY(const ZMatrix &m, const DVector y);

void AddEachColinMatrixWithY(const DMatrix &m, const DVector y);

void AddEachColinMatrixWithY(const DMatrix &m, const WVector y);

void AddEachColinMatrixWithY(const WMatrix &m, const WVector y);

void AddEachLineinMatrixWithY(const ZMatrix &m, const ZVector y);

void AddEachLineinMatrixWithY(const ZMatrix &m, const WVector y);

void AddEachLineinMatrixWithY(const ZMatrix &m, const DVector y);

void AddEachLineinMatrixWithY(const DMatrix &m, const DVector y);

void AddEachLineinMatrixWithY(const DMatrix &m, const WVector y);

void AddEachLineinMatrixWithY(const WMatrix &m, const WVector y);
	
void MultEachColinMatrixWithY(const ZMatrix &m, const ZVector y);
	
void MultEachColinMatrixWithY(const ZMatrix &m, const WVector y);
	
void MultEachColinMatrixWithY(const ZMatrix &m, const DVector y);

void MultEachColinMatrixWithY(const DMatrix &m, const DVector y);

void MultEachColinMatrixWithY(const DMatrix &m, const WVector y);

void MultEachColinMatrixWithY(const WMatrix &m, const WVector y);

void MultEachLineinMatrixWithY(const ZMatrix &m, const ZVector y);

void MultEachLineinMatrixWithY(const ZMatrix &m, const WVector y);

void MultEachLineinMatrixWithY(const ZMatrix &m, const DVector y);

void MultEachLineinMatrixWithY(const DMatrix &m, const DVector y);

void MultEachLineinMatrixWithY(const DMatrix &m, const WVector y);

void MultEachLineinMatrixWithY(const WMatrix &m, const WVector y);



ZMatrix InterleaveTwoMatrixLineWise(const ZMatrix &z1, const ZMatrix &z2);
	
DMatrix InterleaveTwoMatrixLineWise(const DMatrix &z1, const DMatrix &z2);
	
WMatrix InterleaveTwoMatrixLineWise(const WMatrix &z1, const WMatrix &z2);

DMatrix exp(const DMatrix &d);

void LowerBound(const DMatrix &d, double lowerlimit);
	
void HigherBound(const DMatrix &d, double higherlimit);

void Bound(const DMatrix &d, double lowerlimit, double higherlimit);

DMatrix ColMax(const DMatrix& d);

DMatrix ColMax(const DMatrix& d, WMatrix &pos);

DMatrix ColMin(const DMatrix& d);

DMatrix ColMin(const DMatrix& d, WMatrix &pos);

DMatrix RowMax(const DMatrix& d);

DMatrix RowMax(const DMatrix& d, WMatrix &pos);

DMatrix RowMin(const DMatrix& d);

DMatrix RowMin(const DMatrix& d, WMatrix &pos);

WMatrix ColMax(const WMatrix& d);

WMatrix ColMax(const WMatrix& d, WMatrix &pos);

WMatrix ColMin(const WMatrix& d);

WMatrix ColMin(const WMatrix& d, WMatrix &pos);

WMatrix RowMax(const WMatrix& d);

WMatrix RowMax(const WMatrix& d, WMatrix &pos);

WMatrix RowMin(const WMatrix& d);

WMatrix RowMin(const WMatrix& d, WMatrix &pos);

DMatrix ColSum(const DMatrix &d);

DMatrix RowSum(const DMatrix &d);

DMatrix ColAvg(const DMatrix &d);

DMatrix RowAvg(const DMatrix &d);

WMatrix ColSum(const WMatrix &d);

WMatrix RowSum(const WMatrix &d);

DMatrix ColAvg(const WMatrix &d);

DMatrix RowAvg(const WMatrix &d);

ZVector extractdiag(const ZMatrix &m);
	
WVector extractdiag(const WMatrix &m);
	
DVector extractdiag(const DMatrix &m);

WVector ExtractDiagonal(const WMatrix &d);

DVector ExtractDiagonal(const DMatrix &d);

ZVector ExtractDiagonal(const ZMatrix &d);

WMatrix ExtractDiagonalAsColumnMatrix(const WMatrix &d);

DMatrix ExtractDiagonalAsColumnMatrix(const DMatrix &d);

ZMatrix ExtractDiagonalAsColumnMatrix(const ZMatrix &d);

WMatrix ExtractDiagonalAsLineMatrix(const WMatrix &d);

DMatrix ExtractDiagonalAsLineMatrix(const DMatrix &d);

ZMatrix ExtractDiagonalAsLineMatrix(const ZMatrix &d);

ZMatrix ColumnCirculant(const ZVector& col);

DMatrix ColumnCirculant(const DVector& col);

WMatrix ColumnCirculant(const WVector& col);

ZMatrix LineCirculant(const ZVector& col);

DMatrix LineCirculant(const DVector& col);

WMatrix LineCirculant(const WVector& col);

WMatrix DiagPackedVectorAsMatrix(const WVector& v, int matline);

DMatrix DiagPackedVectorAsMatrix(const DVector& v, int matline);

ZMatrix DiagPackedVectorAsMatrix(const ZVector& v, int matline);

#endif
