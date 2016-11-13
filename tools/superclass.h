/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/07/21 10:59:16 $
 * $Revision: 1.1.2.12 $
 * $Id: superclass.h,v 1.1.2.12 2006/07/21 10:59:16 syed Exp $
 ********************************************************************/

#ifndef _SUPERCLASS_C_
#define _SUPERCLASS_C_

#include <iostream>
using namespace std;

#define SUPERCLASS_EMPTY_VECTOR  1
#define SUPERCLASS_INCOMPATIBLE_VECTOR_SIZE  2

class CSuperclassException{
public :
	int err_mesg;
	CSuperclassException(int code) : err_mesg(code){}
};


template<class T> class UVector {
	public:
		int taille;
		T *vect;
		UVector<T>(int size=0);
		UVector<T>(const UVector<T> &u);
		~UVector<T>();
		UVector<T> operator+(const UVector<T> &r);
		UVector<T> operator-(const UVector<T> &r);
		UVector<T> operator*(const UVector<T> &r);
		UVector<T> &operator=(const UVector<T> &r);
		UVector<T> &operator+=(const UVector<T> &r);
		UVector<T> &operator-=(const UVector<T> &r);
		bool operator==(const UVector<T> &r);
		bool operator!=(const UVector<T> &r);
};

template<class T> class UMatrix{
	public:
		int line, col;
		T ** mat;

		UMatrix<T>(int line, int col);
		UMatrix<T>();
		UMatrix<T>(const UMatrix<T> &u);
		~UMatrix<T>();
		UMatrix<T> operator+(const UMatrix<T> &r);
		UMatrix<T> operator-(const UMatrix<T> &r);
		UMatrix<T> operator*(const UMatrix<T> &r);
		UMatrix<T> &operator=(const UMatrix<T> &r);
		bool operator==(const UMatrix<T> &r);
		bool operator!=(const UMatrix<T> &r);
		void transpose();
		UMatrix<T> transpose_nip();
};

template<class T> UVector<T>::UVector(int size){
	if (size<=0) {
		vect=NULL;taille=0;
		return;
	}
	vect = new T[size];
	taille=size;
}

template<class T> UVector<T>::~UVector() {
	if (taille>0 && vect != NULL) {
		delete[] vect;
		vect=NULL;
		taille=0;
	}
}

template<class T> UVector<T> UVector<T>::operator+(const UVector<T> &r){
	if(taille!=r.taille) return;
	UVector<T> nv(taille);
	for(register int i = 0; i < taille; i++) {
		nv.vect[i]=vect[i]+r.vect[i];
	}
	return nv;
}

template<class T> UVector<T> UVector<T>::operator-(const UVector<T> &r){
	if(taille!=r.taille) return;
	UVector<T> nv(taille);
	for(register int i = 0; i < taille; i++) {
		nv.vect[i]=vect[i]-r.vect[i];
	}
	return nv;
}

template<class T> UVector<T> UVector<T>::operator*(const UVector<T> &r){
	if(taille!=r.taille) return;
	UVector<T> nv(taille);
	for(register int i = 0; i < taille; i++) {
		nv.vect[i]=vect[i]*r.vect[i];
	}
	return nv;
}

template<class T> UVector<T> &UVector<T>::operator=(const UVector<T> &r){
	this->~UVector<T>();
	taille=r.taille;
	if(r.vect!=NULL && taille>0){
		vect = new T[taille];
		for(register int i=0;i<taille;i++) {
			vect[i]=r.vect[i];
		}
	}
	return *this;
}

template <class T> UVector<T> &UVector<T>::operator+=(const UVector<T> &r){
         if (r.vect==NULL || taille<=0) throw CSuperclassException(SUPERCLASS_EMPTY_VECTOR);
         if (vect==NULL || taille<=0) {
            *this = r;
            return *this;
         }
         if(r.taille!=taille) throw CSuperclassException(SUPERCLASS_INCOMPATIBLE_VECTOR_SIZE);
         for (register int i=0;i<taille;i++) {
             vect[i]+=r.vect[i];
         }
         return *this;
}

template<class T> UVector<T> &UVector<T>::operator-=(const UVector<T> &r){
         if (r.vect==NULL || taille<=0) throw CSuperclassException(SUPERCLASS_EMPTY_VECTOR);
         if (vect==NULL || taille<=0) {
            *this = r;
            return *this;
         }
         if(r.taille!=taille) throw CSuperclassException(SUPERCLASS_INCOMPATIBLE_VECTOR_SIZE);
         for (register int i=0;i<taille;i++) {
             vect[i]-=r.vect[i];
         }
         return *this;
}

template<class T> bool UVector<T>::operator==(const UVector<T> &r){
	if (taille!=r.taille) return false;
	if(r.vect!=NULL && taille>0){
		for(register int i=0;i<taille;i++) {
			if (vect[i]!=r.vect[i]) return false;
		}
	}
	return true;
}

template<class T> bool UVector<T>::operator!=(const UVector<T> &r){
	if (taille!=r.taille) return true;
	if(r.vect!=NULL && taille>0){
		for(register int i=0;i<taille;i++) {
			if (vect[i]!=r.vect[i]) return true;
		}
	}
	return false;
}

template<class T> UVector<T>::UVector(const UVector<T> &u){
	if (u.taille<=0 || u.vect==NULL) {
		vect=NULL;taille=0;
		return;
	}
	vect = new T[u.taille];
	taille=u.taille;
	for(register int i=0;i<u.taille;i++) {
		vect[i]=u.vect[i];
	}
}

template<class T> UMatrix<T>::UMatrix(int line, int col){
	if (col==0 || line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=col;this->line=line;
	mat=new T*[line];
	for(register int i =0;i<line;i++) {
		mat[i]=new T[col];
	}
}

template<class T> UMatrix<T>::UMatrix(){
	col=line=0;mat=NULL;
}

template<class T> UMatrix<T>::~UMatrix(){
	for(register int i=0;i<line;i++)
		delete[] mat[i];
	delete[] mat;
	line=col=0;
	mat=NULL;
}

template<class T> void UMatrix<T>::transpose() {
	if (mat==NULL ||col==0 ||line==0) return;
	T **tmp = new T*[col];
	register int i=0;
	for(i=0; i<col;i++) {
		tmp[i]=new T[line];
		for(register int j=0;j<line;j++)
			tmp[i][j]=mat[j][i];
	}
	for(i=0;i<line;i++) delete[] mat[i];
	delete[] mat;

	int x=col;col=line;line=x;
	mat=tmp;
}

template<class T> UMatrix<T> UMatrix<T>::transpose_nip(){
	if (mat==NULL ||col==0 ||line==0) return UMatrix();
	UMatrix t; t.col=line;t.line=col;

	t.mat = new T*[col];
	for(register int i=0; i<col;i++) {
		t.mat[i]=new T[line];
		for(register int j=0;j<line;j++)
			t.mat[i][j]=mat[j][i];
	}
	return t;
}

template <class T> UMatrix<T> UMatrix<T>::operator+(const UMatrix<T> &rhs){
	if (col==0 || line==0 || mat==NULL) return UMatrix();
	if (rhs.col!=col || rhs.line!=line) return UMatrix();
	UMatrix<T> n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]+rhs.mat[i][j];
	}
	return n;
}

template <class T> UMatrix<T> UMatrix<T>::operator-(const UMatrix<T> &rhs){
	if (col==0 || line==0 || mat==NULL) return UMatrix();
	if (rhs.col!=col || rhs.line!=line) return UMatrix();
	UMatrix<T> n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]-rhs.mat[i][j];
	}
	return n;
}

template <class T> UMatrix<T> UMatrix<T>::operator*(const UMatrix<T> &rhs){
	if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return UMatrix();
	if (rhs.line!=col) return UMatrix();
	UMatrix<T> n(line, rhs.col);

	for(register int j=0;j<n.line;j++) {
		for(register int i=0;i<n.col;i++) {
			register int sum=0;
			for(register int k=0;k<col;k++)
				sum+=mat[j][k]*rhs.mat[k][i];
			n.mat[j][i]=sum;
		}
	}
	return n;
}

template <class T> UMatrix<T> &UMatrix<T>::operator=(const UMatrix<T> &rhs){
	if(mat!=NULL) {
		for(register int i=0;i<line;i++) delete[] mat[i];
		delete[] mat;
	}

	line=rhs.line;col=rhs.col;

	if(line==0|| col==0 || rhs.mat==NULL) {
		this->mat=NULL; line=col=0; return *this;
	}

	mat=new T*[line];
	for(register int i = 0;i<line;i++) {
		mat[i]=new T[col];
		for(register int j=0;j<col;j++)
			mat[i][j]=rhs.mat[i][j];
	}
	return *this;

}

template <class T> bool UMatrix<T>::operator==(const UMatrix<T> &rhs){
	if(rhs.col!=col || rhs.line!=line) return false;
	for(int i=0;i<line;i++){
		for(int j=0;j<col;j++)
			if (rhs.mat[i][j]!=mat[i][j]) return false;
	}
	return true;
}

template <class T> bool UMatrix<T>::operator!=(const UMatrix<T> &rhs){
	if(rhs.col!=col || rhs.line!=line) return true;
	for(int i=0;i<line;i++){
		for(int j=0;j<col;j++)
			if (rhs.mat[i][j]!=mat[i][j]) return true;
	}
	return false;
}

template<class T> UMatrix<T>::UMatrix(const UMatrix<T> &u){
	if (u.line<=0 || u.col<=0 || u.mat==NULL) {
		mat=NULL;line=col=0;
		return;
	}
	mat = new T*[u.line];
	line=u.line;col=u.col;
	for(register int i=0;i<line;i++) {
		mat[i]=new T[u.col];
		for(register int j=0;j<col;j++) {
			mat[i][j]=u.mat[i][j];
		}
	}
}

template<class T> istream& operator>>(istream &is, UVector<T> &u) {
	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	int size ;
	char junk;
	is >> size >> junk >> junk;
	if(size<=0) { 
		u = UVector<T>();
		return is;
	}
	u = UVector<T>(size);
	for(int i=0;i<size;i++) {
		is >> u.vect[i];
	}
	is >> junk;
	return is;
}

template<class T> ostream& operator<<(ostream &os, UVector<T> &u) {
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	os << u.taille << " : {" << endl;
	for(int i=0;i<u.taille;i++) {
		os << u.vect[i];
	}
	os << "}"<<endl;
	return os;
}

template<class T> istream& operator>>(istream &is, UMatrix<T> &u) {
	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	int line, col ;
	char junk;
	is >> line >> junk >> col >> junk >> junk;
	if(line<=0 || col<=0) { 
		u = UMatrix<T>();
		return is;
	}
	u = UMatrix<T>(line,col);
	for(int i=0;i<line;i++) {
		for(int j=0;j<col;j++) {
			is >> u.mat[i][j];
		}
	}
	is >> junk;
	return is;
}

template<class T> ostream& operator<<(ostream &os, UMatrix<T> &u) {
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	os << u.line << "," << u.col << " : {" << endl;
	for(int i=0;i<u.line;i++) {
		for(int j=0;j<u.col;j++) {
			os << u.mat[i][j];
		}
	}
	os << "}"<<endl;
	return os;
}

#endif
