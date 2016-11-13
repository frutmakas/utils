/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/07/21 10:58:28 $
 * $Revision: 1.1.2.11 $
 * $Id: binmat.cpp,v 1.1.2.11 2006/07/21 10:58:28 syed Exp $
 ********************************************************************/
#include "tools/tools.h"
#include "tools/binmat.h"
#include "globaldef.h"

#ifdef xWIN_DOS
#include <iostream.h>
#include <fstream.h>
#else
#include <iostream>
#include <fstream>
using namespace std;
#endif

const int word_size=sizeof(int)*8;
const int word_mask=-1;
const int word_shift=(int)log_2(word_size);

#ifdef WIN_DOS
//#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ )
#endif

BVector::BVector(int nb_element){
	real_size=nb_element;
	word_count=(int)ceil((double)(nb_element/word_size))+1;
	vect = new int[word_count];
	for(register int i=0;i<word_count;i++) vect[i]=0;
}
BVector::BVector(){
	real_size=word_count=0;
	vect = NULL;
}

BVector::BVector(const BVector &bv){
	if(bv.vect==NULL) throw BVectorException(INVALID_VECTOR);
	real_size=bv.real_size;
	word_count=bv.word_count;
	vect=new int[word_count];
	for(register int i=0;i<word_count;i++) vect[i]=bv.vect[i];
}

BVector::~BVector(){
	real_size=word_count=0;
	if(vect!=NULL) delete[] vect;
	vect=NULL;
}

BVector& BVector::operator=(const BVector &bv){
	if(vect!=NULL) delete[] vect;
	if(bv.vect==NULL) {
		vect=NULL;
		real_size=word_count=0;
		return *this;
	}
	real_size=bv.real_size;
	word_count=bv.word_count;
	vect=new int[word_count];
	for(register int i=0;i<word_count;i++) vect[i]=bv.vect[i];
	return *this;
}

BVector& BVector::operator^=(const BVector &bv){
	if(vect==NULL || bv.vect==NULL) throw BVectorException(INVALID_VECTOR);
	if(real_size!=bv.real_size) throw BVectorException(INVALID_VECTOR_DIMENSION);

	for(register int i=0;i<word_count;i++) vect[i]^=bv.vect[i];
	return *this;
}


int BVector::operator[](int position) const{
	if (position >= real_size) throw BVectorException(INDEX_OUT_OF_RANGE);
	int idx=position>>word_shift;
	return (vect[idx]&(1 << (position&word_mask)))==0?0:1;
}

bool BVector::operator==(const BVector &bv) const{
	if(bv.vect==NULL) throw BVectorException(INVALID_VECTOR);
	for(register int i=0; i<word_count;i++) {
		if (bv.vect[i]!=vect[i]) return false;
	}
	return true;
}

bool BVector::operator!=(const BVector &bv) const{
	if(bv.vect==NULL) throw BVectorException(INVALID_VECTOR);
	for(register int i=0; i<word_count;i++) {
		if (bv.vect[i]!=vect[i]) return true;
	}
	return false;
}

BVector BVector::operator+(const BVector &bv) const { // modulo 2 addition
	if(bv.vect==NULL) throw BVectorException(INVALID_VECTOR);
	BVector nv(real_size);
	for(register int i=0; i<word_count;i++) {
		nv.vect[i]=bv.vect[i]^vect[i];
	}
	return nv;

}

BVector BVector::operator^(int boolean) const{
	int op2=boolean==0?0:word_size-1;;
	BVector nv(real_size);
	for(register int i=0; i<word_count;i++) {
		nv.vect[i]=vect[i]^op2;
	}
	return nv;
}

BVector BVector::operator^(const BVector &bv) const{
	BVector nv(real_size);
	for(register int i=0; i<word_count;i++) {
		nv.vect[i]=vect[i]^bv.vect[i];
	}
	return nv;
}

BVector BVector::operator~() const{
	BVector nv(real_size);
	for(register int i=0; i<word_count;i++) {
		nv.vect[i]=~vect[i];
	}
	return nv;
}

BVector BVector::operator!() const{
	BVector nv(real_size);
	for(register int i=0; i<word_count;i++) {
		nv.vect[i]=vect[i]^word_mask;
	}
	return nv;
}

BVector BVector::operator&(const BVector &bv)const{
	BVector nv(real_size);
	for(register int i=0; i<word_count;i++) {
		nv.vect[i]=vect[i]&bv.vect[i];
	}
	return nv;
}

BVector BVector::operator|(const BVector &bv) const{
	BVector nv(real_size);
	for(register int i=0; i<word_count;i++) {
		nv.vect[i]=vect[i]|bv.vect[i];
	}
	return nv;
}

void BVector::Set(int position) const{
	int idx=position>>word_shift;
	vect[idx]|=(1 << (position&word_mask));
}

void BVector::Reset(int position) const{
	int idx=position>>word_shift;
	int reset_mask=(1 << (position&word_mask))^word_mask;
	vect[idx]&=reset_mask;
}

void BVector::Toggle(int position) const{
	if(this->operator[](position)) Reset(position); else Set(position);
}

void BVector::Clear() const{
	for(register int i=0;i<word_count;i++) vect[i]=0;
}

ostream& operator<<(ostream &os, const BVector& bv){
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	os << endl << bv.real_size << " : [";
	for(register int i=0; i < bv.real_size;i++) {
		os << bv[i] << " ";
	}
	os << "]";
	return os;
}

istream& operator>>(istream &is, BVector &bv) {
	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	int size, bit_value;
	is >> size;
	char junk;
	for(register int j=0;j<2;j++) {
		is>>junk;
	}
	bv=BVector(size);
	for(register int i=0; i <size; i++) {
		is >> bit_value;
		if (bit_value==1) bv.Set(i); else bv.Reset(i);
	}
	is >> junk;
	return is;
}

BMatrix::BMatrix(int line, int col){
	this->line=line;this->col=col;
	mat = new BVector[line];
	for(register int i = 0; i < line;i++) {
		mat[i]=BVector(col);
	}
}

BMatrix::BMatrix(){
	line=col=0; mat=NULL;
}

BMatrix::BMatrix(const BMatrix &bv){
	line=bv.line;col=bv.col;
	mat=new BVector[line];
	for(register int i=0;i<line;i++) {
		mat[i]=bv.mat[i];
	}
}

BMatrix::~BMatrix(){
	if(mat!=NULL) delete[] mat;
	line=col=0;mat=NULL;
}

bool BMatrix::operator==(const BMatrix &bv) const{
	if(line!=bv.line||col!=bv.col) throw BMatrixException(INCOMPATIBLE_MATRIX_DIMENSION);
	if(mat==NULL || bv.mat==NULL) throw BMatrixException(INVALID_MATRIX);
	for(register int i =0; i< line;i++){
		if(mat[i]!=bv.mat[i]) return false;
	}
	return true;
}

bool BMatrix::operator!=(const BMatrix &bv) const{
	if(line!=bv.line||col!=bv.col) throw BMatrixException(INCOMPATIBLE_MATRIX_DIMENSION);
	if(mat==NULL || bv.mat==NULL) throw BMatrixException(INVALID_MATRIX);
	for(register int i =0; i< line;i++){
		if(mat[i]!=bv.mat[i]) return true;
	}
	return false;
}

BMatrix& BMatrix::operator=(const BMatrix &bv){
	if(bv.mat==NULL) {
		delete[] mat;
		line=col=0;mat=NULL;
		return *this;
	}
	line=bv.line;col=bv.col;
	mat=new BVector[line];
	for(register int i=0;i<line;i++) {
		mat[i]=bv.mat[i];
	}
	return *this;
}

BMatrix& BMatrix::operator=(const BIMatrix &bv){
    if(mat!=NULL) delete[] mat;
	if(bv.mat==NULL) {
//		delete[] mat;
		line=col=0;mat=NULL;
		return *this;
	}
	line=bv.line;col=bv.col;
	mat=new BVector[line];
	for(register int i=0;i<line;i++) {
		mat[i]=BVector(col);
		for(register int j=0;j<col;j++) {
			if(bv.mat[j][i]) Set(i,j);
		}
	}
	return *this;
}


BMatrix BMatrix::operator+(const BMatrix &bv) const{
	if(line!=bv.line||col!=bv.col) throw BMatrixException(INCOMPATIBLE_MATRIX_DIMENSION);
	if(mat==NULL || bv.mat==NULL) throw BMatrixException(INVALID_MATRIX);
	BMatrix nm(line,col);
	for(register int i =0; i< line;i++){
		nm.mat[i]=mat[i]+bv.mat[i];
	}
	return nm;
}

BMatrix BMatrix::operator^(const BMatrix &bv) const{
	if(line!=bv.line||col!=bv.col) throw BMatrixException(INCOMPATIBLE_MATRIX_DIMENSION);
	if(mat==NULL || bv.mat==NULL) throw BMatrixException(INVALID_MATRIX);
	BMatrix nm(line,col);
	for(register int i =0; i< line;i++){
		nm.mat[i]=mat[i]^bv.mat[i];
	}
	return nm;
}

BMatrix BMatrix::operator~() const{
	if(mat==NULL) return BMatrix();
	BMatrix nm(line,col);
	for(register int i =0; i< line;i++){
		nm.mat[i]=~mat[i];
	}
	return nm;
}

BMatrix BMatrix::operator&(const BMatrix &bv) const{
	if(line!=bv.line||col!=bv.col) throw BMatrixException(INCOMPATIBLE_MATRIX_DIMENSION);
	if(mat==NULL || bv.mat==NULL) throw BMatrixException(INVALID_MATRIX);
	BMatrix nm(line,col);
	for(register int i =0; i< line;i++){
		nm.mat[i]=mat[i]&bv.mat[i];
	}
	return nm;
}

BMatrix BMatrix::operator|(const BMatrix &bv) const{
	if(line!=bv.line||col!=bv.col) throw BMatrixException(INCOMPATIBLE_MATRIX_DIMENSION);
	if(mat==NULL || bv.mat==NULL) throw BMatrixException(INVALID_MATRIX);
	BMatrix nm(line,col);
	for(register int i =0; i< line;i++){
		nm.mat[i]=mat[i]|bv.mat[i];
	}
	return nm;
}

BMatrix BMatrix::operator*(const BMatrix &bv) const{
	if(bv.mat==NULL || mat==NULL) throw BMatrixException(INVALID_MATRIX);
	if(col!=bv.line) throw BMatrixException(INCOMPATIBLE_MATRIX_DIMENSION);
	BMatrix n(line,bv.col);
	for(register int j=0;j<n.line;j++) { //this.line
		for(register int i=0;i<n.col;i++) { // rhs.col
			int sum=0;
			for(register int k=0;k<col;k++)  // this.col
				sum+=mat[j][k]*bv.mat[k][i];
			if(sum%2) n.Set(j,i); 
		}
	}
	return n;
}


void BMatrix::Set(int line, int col) {
	mat[line].Set(col);
}

void BMatrix::Reset(int line, int col) {
	mat[line].Reset(col);
}

void BMatrix::Toggle(int line, int col) {
	mat[line].Toggle(col);
}


ostream& operator<<(ostream &os, const BMatrix& bv){
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	if(bv.mat==NULL) {
		cerr << endl << "Empty matrix!" << endl;
	}
	os << bv.line << "," << bv.col << " : [" ; 
	for(register int m=0;m<bv.line;m++) {
		os << bv.mat[m];
	}
	os << endl << "]" << endl;
	return os;
}

istream& operator>>(istream &is, BMatrix &bv){
	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	int l,c;
	char junk; 
	is >> l >> junk >> c >> junk >> junk;
	bv.~BMatrix();
	bv=BMatrix(l,c);
	for(register int m=0;m<l;m++){
		is >> bv.mat[m];
		//is >> junk;
	}
	is >> junk;
	return is;
}


/*************************************************/

BIMatrix::BIMatrix(int line, int col){
	this->line=line;this->col=col;
	mat = new BVector[col];
	for(register int i = 0; i < col;i++) {
		mat[i]=BVector(line);
	}
}

BIMatrix::BIMatrix(){
	line=col=0; mat=NULL;
}

BIMatrix::BIMatrix(const BIMatrix &bv){
	line=bv.line;col=bv.col;
	mat=new BVector[col];
	for(register int i=0;i<col;i++) {
		mat[i]=bv.mat[i];
	}
}

BIMatrix::~BIMatrix(){
	if(mat!=NULL) delete[] mat;
	line=col=0;mat=NULL;
}

bool BIMatrix::operator==(const BIMatrix &bv) const{
	if(line!=bv.line||col!=bv.col) throw BIMatrixException(INCOMPATIBLE_MATRIX_DIMENSION);
	if(mat==NULL || bv.mat==NULL) throw BIMatrixException(INVALID_MATRIX);
	for(register int i =0; i< col;i++){
		if(mat[i]!=bv.mat[i]) return false;
	}
	return true;
}

bool BIMatrix::operator!=(const BIMatrix &bv) const{
	if(line!=bv.line||col!=bv.col) throw BIMatrixException(INCOMPATIBLE_MATRIX_DIMENSION);
	if(mat==NULL || bv.mat==NULL) throw BIMatrixException(INVALID_MATRIX);
	for(register int i =0; i< col;i++){
		if(mat[i]!=bv.mat[i]) return true;
	}
	return false;
}

BIMatrix& BIMatrix::operator=(const BIMatrix &bv){
    if(mat!=NULL) delete[] mat;
	if(bv.mat==NULL) {
//		delete[] mat;
		line=col=0;mat=NULL;
		return *this;
	}
	line=bv.line;col=bv.col;
	mat=new BVector[col];
	for(register int i=0;i<col;i++) {
		mat[i]=bv.mat[i];
	}
	return *this;
}

BIMatrix BIMatrix::operator+(const BIMatrix &bv) const{
	if(line!=bv.line||col!=bv.col) throw BIMatrixException(INCOMPATIBLE_MATRIX_DIMENSION);
	if(mat==NULL || bv.mat==NULL) throw BIMatrixException(INVALID_MATRIX);
	BIMatrix nm(line,col);
	for(register int i =0; i< col;i++){
		nm.mat[i]=mat[i]+bv.mat[i];
	}
	return nm;
}

BIMatrix BIMatrix::operator^(const BIMatrix &bv) const{
	if(line!=bv.line||col!=bv.col) throw BIMatrixException(INCOMPATIBLE_MATRIX_DIMENSION);
	if(mat==NULL || bv.mat==NULL) throw BIMatrixException(INVALID_MATRIX);
	BIMatrix nm(line,col);
	for(register int i =0; i< col;i++){
		nm.mat[i]=mat[i]^bv.mat[i];
	}
	return nm;
}

BIMatrix BIMatrix::operator~() const{
	if(mat==NULL) return BIMatrix();
	BIMatrix nm(line,col);
	for(register int i =0; i< col;i++){
		nm.mat[i]=~mat[i];
	}
	return nm;
}

BIMatrix BIMatrix::operator&(const BIMatrix &bv) const{
	if(line!=bv.line||col!=bv.col) throw BIMatrixException(INCOMPATIBLE_MATRIX_DIMENSION);
	if(mat==NULL || bv.mat==NULL) throw BIMatrixException(INVALID_MATRIX);
	BIMatrix nm(line,col);
	for(register int i =0; i< col;i++){
		nm.mat[i]=mat[i]&bv.mat[i];
	}
	return nm;
}

BIMatrix BIMatrix::operator|(const BIMatrix &bv) const{
	if(line!=bv.line||col!=bv.col) throw BIMatrixException(INCOMPATIBLE_MATRIX_DIMENSION);
	if(mat==NULL || bv.mat==NULL) throw BIMatrixException(INVALID_MATRIX);
	BIMatrix nm(line,col);
	for(register int i =0; i< col;i++){
		nm.mat[i]=mat[i]|bv.mat[i];
	}
	return nm;
}

BIMatrix BIMatrix::operator*(const BIMatrix &bv) const{
	if(bv.mat==NULL || mat==NULL) throw BIMatrixException(INVALID_MATRIX);
	if(col!=bv.line) throw BIMatrixException(INCOMPATIBLE_MATRIX_DIMENSION);
	BIMatrix n(line,bv.col);
	for(register int j=0;j<n.line;j++) { //this.line
		for(register int i=0;i<n.col;i++) { // rhs.col
			int sum=0;
			for(register int k=0;k<col;k++)  // this.col
				sum+=mat[k][j]*bv.mat[i][k];
			if(sum%2) n.Set(j,i); 
		}
	}
	return n;
}


void BIMatrix::Set(int line, int col) {
	mat[col].Set(line);
}

void BIMatrix::Reset(int line, int col) {
	mat[col].Reset(line);
}

void BIMatrix::Toggle(int line, int col) {
	mat[col].Toggle(line);
}


ostream& operator<<(ostream &os, const BIMatrix& bv){
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	if(bv.mat==NULL) {
		cerr << endl << "Empty matrix!" << endl;
	}
	os << bv.line << "," << bv.col << " : [" ; 
	for(register int m=0;m<bv.col;m++) {
		os << bv.mat[m];
	}
	os << endl << "]" << endl;
	return os;
}

istream& operator>>(istream &is, BIMatrix &bv){
	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	int l,c;
	char junk; 
	is >> l >> junk >> c >> junk >> junk;
	bv.~BIMatrix();
	bv=BIMatrix(l,c);
	for(register int m=0;m<c;m++){
		is >> bv.mat[m];
		//is >> junk;
	}
	is >> junk;
	return is;
}

void BVector::Force(int position, bool value)
{
	if(position>=real_size) return; // 
	if(value) Set(position);
	else Reset(position);

}
