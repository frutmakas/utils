/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/07/21 10:59:52 $
 * $Revision: 1.1.2.55 $
 * $Id: utilitis.h,v 1.1.2.55 2006/07/21 10:59:52 syed Exp $
 ********************************************************************/

#ifndef _UTILITIS_H_
#define _UTILITIS_H_

/*#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>*/
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

#include "globaldef.h"


#ifndef UTILITIS_M_2_PI 
#define UTILITIS_M_2_PI 6.28318530717959
#endif

#ifndef DOUBLE_EPS
#define DOUBLE_EPS DBL_EPSILON
#endif

#ifndef M_HUGE_NUM
#define M_HUGE_NUM  1e50
#endif  

#ifndef M_PI
#define M_PI   3.14159265358979
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679490
#endif

class WVector;
class ZVector;
class DVector;

class WMatrix;
class ZMatrix;
class DMatrix;


class CUtilitisException{
public : 
	int err_mesg;
	CUtilitisException(int mesg) : err_mesg(mesg) {
		cerr << "CUtilitisException : " << mesg << endl;
		printerror();

	}
	CUtilitisException(const char * filename, int linenumber, int mesg) : err_mesg(mesg) {
		cerr << "CUtilitisException : " << mesg << " at "<< filename << "@" << linenumber << endl;
		printerror();
	} 
protected:
	void printerror() {
		switch (err_mesg) {
			case BAD_FILE : 
				cerr << "CUtilitisException : File descriptor is bad" << endl;
				break;
			case INDEX_OUT_OF_RANGE	: 
				cerr << "CUtilitisException : Index is not within permitted range" << endl;
				break;
 			case INVALID_RANGE	: 
				cerr << "CUtilitisException : Range given is invalid" << endl;
				break;
 			case DECOMPOSE_ERROR	: 
				cerr << "CUtilitisException : Error while decomposing structure " << endl;
				break;
 			case EMPTY_MATRIX	: 
				cerr << "CUtilitisException : The matrix is empty" << endl;
				break;
 			case DIVISION_BY_ZERO	: 
				cerr << "CUtilitisException : Division by zero" << endl;
				break;
 			case INVALID_DIMENSION	: 
				cerr << "CUtilitisException : Dimension is invalid" << endl;
				break;
 			case INVALID_ARGUMENT	: 
				cerr << "CUtilitisException : Argument given is invalid" << endl;
				break;
 			case EMPTY_VECTOR	: 
				cerr << "CUtilitisException : The vector is empty" << endl;
				break;
 			case INVALID_VECTOR	: 
				cerr << "CUtilitisException : The vector is invalid " << endl;
				break;
 			case INVALID_VECTOR_DIMENSION	: 
				cerr << "CUtilitisException : The vector size is invalid" << endl;
				break;
 			case INCOMPATIBLE_VECTOR_DIMENSION	: 
				cerr << "CUtilitisException : The vector dimension is incompatible" << endl;
				break;
 			case  INCOMPATIBLE_MATRIX_DIMENSION	: 
				cerr << "CUtilitisException : The matrix dimension is incompatible" << endl;
				break;
 			case INVALID_MATRIX	: 
				cerr << "CUtilitisException : The matrix is invalid" << endl;
				break;
 			case INVALID_MATRIX_DIMENSION	: 
				cerr << "CUtilitisException : The matrix dimension is invalid" << endl;
				break;
 			case MATRIX_INDEX_OUT_OF_RANGE	: 
				cerr << "CUtilitisException : Matrix index is out of permitted range" << endl;
				break;
 			case MATRIX_NOT_SQUARE	: 
				cerr << "CUtilitisException : The matrix is not square" << endl;
				break;
 			case  INVALID_OFDM_DATA_SIZE_TO_DECODE	: 
				cerr << "CUtilitisException : Decode fails as the OFDM data size invalid" << endl;
				break;
 			case INVALID_FFT_SIZE	: 
				cerr << "CUtilitisException : FFT size is not power of 2" << endl;
				break;
 			case  MEMORY_ALLOCATION_ERROR	: 
				cerr << "CUtilitisException : Error while allocating memory" << endl;
				break;
 			case OUT_OF_MEMORY	: 
				cerr << "CUtilitisException : Not enough memory to continue" << endl;
				break;
 			case  INVALID_UVECTOR	: 
				cerr << "CUtilitisException : UVector is invalid" << endl;
				break;
 			case INCONSISTENT_UVECTOR	: 
				cerr << "CUtilitisException : Element in UVector is not consistent to continue operation" << endl;
				break;
 			case INCONSISTENT_MATRIX_SIZE_IN_UVECTOR	: 
				cerr << "CUtilitisException : Matrix dimension in UVector is not consistent" << endl;
				break;
 			case EMPTY_UVECTOR	: 
				cerr << "CUtilitisException : UVector is empty" << endl;
				break;
			case INVALID_OPERATION : 
				cerr << "CUtilitisException : Assigning a number to a vector is not allowed. Use Set or Reset" << endl;
				break;
			default:
				cerr << "CUtilitisException : Unknown error. Refer to error definition" << endl;
				break;
		}
	}

};

class CMatrixInversionException{
	public:

		int ErrorCode;
		CMatrixInversionException(const char * filename, int linenumber, int mesg) {
			ErrorCode = mesg;
			cerr << "CMatrixInversionException : " << mesg << " at "<< filename << "@" << linenumber << endl;
			printerror();
		} 
		CMatrixInversionException(int code):ErrorCode(code){
			cerr << "CMatrixInversionException : " << code << endl; //<< " at "<< filename << "@" << linenumber << endl;
			printerror();
		}
protected:
	void printerror() {
		switch (ErrorCode) {
 			case INVALID_MATRIX_DIMENSION	: 
				cerr << "CMatrixInversionException : The matrix dimension is invalid" << endl;
				break;
 			case  SINGULAR_MATRIX_ERROR	: 
				cerr << "CMatrixInversionException : Matrix is singular" << endl;
				break;
 			case SINGULAR_MATRIX_1_ERROR	: 
				cerr << "CMatrixInversionException : Matrix is singular (Error type 1)" << endl;
				break;
 			case SINGULAR_MATRIX_2_ERROR	: 
				cerr << "CMatrixInversionException : Matrix is singular (Error type 2)" << endl;
				break;
			default:
				cerr << "CMatrixInversionException : Unknown error. Refer to error definition" << endl;
				break;
		}
	}

};


class DCplx {
	public:
		double re, im;
		DCplx(const double re, const double im);
		DCplx();
		DCplx(const DCplx& c);
		inline DCplx conj() const {return DCplx(re,-im);}

		inline double mag() const { return sqrt(re*re+im*im); }
		inline double mod2() const { return re*re+im*im; }

		inline double phase() const {
			if (re==0 && im==0) {
				return 0.0;
			}
			if (re==0 && im!=0) { 
				return (im>0 ) ?  M_PI_2 : -M_PI_2; 
			}
			return atan(im/re);

		}

		inline DCplx operator+(const DCplx &c) const{ return DCplx(re+c.re, im+c.im); }
		inline DCplx operator+(const double c)const { return DCplx(re+c, im); }

		inline DCplx operator-(const DCplx &c) const{ return DCplx(re-c.re, im-c.im); }
		inline DCplx operator-(const double c) const{ return DCplx(re-c, im); }

		inline DCplx operator*(const DCplx &c) const{ return DCplx(re*c.re-c.im*im, im*c.re+re*c.im); }
		inline DCplx operator*(const double c) const{ return DCplx(re*c, c*im); }

		inline DCplx operator/(const double c) const{ if (c!=0.) return DCplx(re/c, im/c); throw CUtilitisException(__FILE__, __LINE__, DIVISION_BY_ZERO);}

		DCplx log(DCplx c);

		inline DCplx operator/(const DCplx &c) const {
			double t=c.mod2();
			if (t==0.0) throw CUtilitisException(__FILE__, __LINE__, DIVISION_BY_ZERO) ;// return DCplx();
			return DCplx((*this)*c.conj()/t);
		}



		static DCplx exp2RI(const double mag ,const double phase) {

			if (mag<=0.0) return DCplx();

			return DCplx(mag*cos(phase), mag*sin(phase));

		}
		
		inline DCplx operator-(void) const { return DCplx(-re, -im); }

		DCplx &operator=(const DCplx &c) ;
		DCplx &operator=(const double c) ;

		bool operator==(const DCplx &c) const;
		bool operator==(const double c) const;

		bool operator!=(const DCplx &c) const;
		bool operator!=(const double c) const;

		DCplx &operator+=(const DCplx &c) ;

		DCplx &operator+=(const double c) ;

		DCplx &operator-=(const DCplx &c) ;
		DCplx &operator-=(const double c) ;

		DCplx &operator*=(const DCplx &c) ;
		DCplx &operator*=(const double  &c);

		void print() const { printf(" (%.3f, %.3f) ", re,im); }
};

inline DCplx operator+(const double c1, const DCplx &c2) {
	return DCplx(c2.re+c1, c2.im);
}

inline DCplx operator-(const double c1, const DCplx &c2) {
	return DCplx(c1-c2.re, c2.im);
}

inline DCplx operator*(const double c1, const DCplx &c2) {
	return DCplx(c2.re*c1, c2.im*c1);
}

inline DCplx operator/(const double c1, const DCplx &c2) {
	return c2.conj()*c1/c2.mod2();
}

inline bool operator==(const double c1, const DCplx &c2) {
	return c1==c2.re && c2.im==0.0;
}

inline bool operator!=(const double c1, const DCplx &c2) {
	return !(c1==c2.re && c2.im==0.0);
}

class WVector{
	public:
        int * vect;
		int taille;

		WVector(const int &size);
		WVector(const int &size, const int &initial) ;
		WVector(const int &size, const int &initial, const char &mode) ;
		WVector(const int *v, const int size) ;
		WVector(const WVector &c);
		WVector() ;
		~WVector() ;
		WVector operator+(const WVector &c) ; 
		DVector operator+(const DVector &c) ; 
		ZVector operator+(const ZVector &c) ; 
		WVector operator+(const int &c) ; 
		DVector operator+(const double c) ; 
		ZVector operator+(const DCplx &c) ; 

		WVector operator-() ; 
		WVector operator-(const WVector &c) ; 
		DVector operator-(const DVector &c) ; 
		ZVector operator-(const ZVector &c) ; 
		WVector operator-(const int &c) ; 
		DVector operator-(const double c) ; 
		ZVector operator-(const DCplx &c) ; 
		WVector operator<<(int n) const;
		WVector operator>>(int n) const;

		WVector operator*(const WVector &c) ; 
		DVector operator*(const DVector &c) ; 
		ZVector operator*(const ZVector &c) ; 
		WVector operator*(const int &c) ; 
		DVector operator*(const double c) ; 
		ZVector operator*(const DCplx &c) ; 
		WVector operator*(const WMatrix &d) const;
		DVector operator*(const DMatrix &d) const;
		ZVector operator*(const ZMatrix &d) const;

		DVector operator/(const double c) ; 
		ZVector operator/(const DCplx &c) ;

		WVector operator~();
		WVector operator^(const WVector &v);


		WVector &operator=(const WVector &c) ; 
		WVector &operator=(int d) ; 

		bool operator==(const DVector &c); 
		bool operator==(const WVector &c) const;
		bool operator==(const ZVector &c) ; 
		bool operator!=(const WVector &c) const ; 
		bool operator!=(const DVector &c) ; 
		bool operator!=(const ZVector &c) ;

		void print() const {
			if(vect==NULL) return;
			printf("\n[");
			for(register int i=0;i<taille;i++) {
					printf(" %d ", vect[i]);
			}
			printf(" ]\n");
		}

		double DotProd(const DVector &rhs);
		int DotProd(const WVector &rhs);
		DCplx DotProd(const ZVector &rhs);
		WVector copy(int start, int stop) const;
		WVector rev_copy(int start, int stop) const;
		void reset(int value=0);
		void insert(const WVector &rhs, int pos) const;
		WVector RZtoNRZ_nip(void);
		void RZtoNRZ(void);
		int sum();
		int diff_count(const WVector& v);
        double variance();

};

class DVector { 
	public:
		double * vect;
		int taille;

		DVector(const int &size,  double initial=0.0) ;
		DVector(const double *v, const int size) ;
		DVector(const WVector &c);
		DVector(const DVector &c);
		DVector() ;
		~DVector() ;
		DVector operator+(const DVector &c) const; 
		DVector operator+(const WVector &c) const ; 
		ZVector operator+(const ZVector &c)  const; 
		ZVector operator+(const DCplx &c)  const; 
		DVector operator+(const double c)  const; 
		DVector operator-() ; 
		DVector operator-(const DVector &c)  const; 
		DVector operator-(const WVector &c)  const; 
		ZVector operator-(const ZVector &c)  const; 
		ZVector operator-(const DCplx &c)  const; 
		DVector operator-(const double c)  const; 
		ZVector operator*(const ZVector &c)  const; 
		DVector operator*(const DVector &c)  const; 
		DVector operator*(const WVector &c)  const; 
		DVector operator*(const DMatrix &d) const;
		DVector operator*(const WMatrix &d) const;
		ZVector operator*(const ZMatrix &d) const;

		ZVector operator*(const DCplx &c)  const; 
		DVector operator*(const double c)  const; 
		DVector operator/(const double c)  const; 
		ZVector operator/(const DCplx &c)  const;
		DVector &operator=(const DVector &c) ; 
		DVector &operator=(const WVector &c) ; 
		DVector &operator=(int d) ; 

		DVector& operator +=(const DVector &dv);
		DVector& operator +=(const WVector &wv);
		DVector& operator +=(double c);

		DVector& operator *=(double c);

		DVector& operator -=(const DVector &dv);
		DVector& operator -=(const WVector &wv);
		DVector& operator -=(double c);
		DVector operator<<(int n) const;
		DVector operator>>(int n) const;

		bool operator==(const DVector &c)  const; 
		bool operator==(const WVector &c)  const; 
		bool operator==(const ZVector &c)  const; 
		bool operator!=(const DVector &c)  const; 
		bool operator!=(const WVector &c)  const; 
		bool operator!=(const ZVector &c)  const; 

		void fixedconv(const DVector &rhs) const;
		void conv(const DVector &rhs);
		DVector conv_nip(const DVector &rhs) const;
		ZVector conv_nip(const ZVector &rhs) const;
		DVector copy(int start, int stop) const ;
		DVector rev_copy(int start, int stop) const ;
		DVector prepad(int count) ;
		DVector postpad(int count);
		void insert(const DVector &rhs, int pos);
		void reset(double value=0.0);
        double variance();

		double sum() const;

		void print() const {
			if(vect==NULL) return;
			printf("\n[");
			for(register int i=0;i<taille;i++) {
					printf(" %.3f ", vect[i]);
			}
			printf(" ]\n");
		}
		double DotProd(const DVector &rhs);
		double DotProd(const WVector &rhs);
		void normalizepower();

		DCplx DotProd(const ZVector &rhs);
};

typedef DVector vector; 

typedef DCplx cmplx;

class ZVector{
	public:
		cmplx * vect;
		int taille;

		ZVector(const double *real, const double* imag, const int &size);
		ZVector(const int &size);
		ZVector(const int &size, const DCplx &initial) ;
		ZVector(const cmplx *v, const int size) ;
		ZVector(const ZVector &c);
		ZVector(const DVector &c);
		ZVector(const WVector &c);
		ZVector() ;
		~ZVector() ;
		ZVector conj() const ;
		DVector real()  const ; 
		DVector imag()  const ; 
		DVector mag()  const ; 
		DVector mod2()  const ; 
		DVector phase()  const ;
		ZVector operator+(const ZVector &c)  const ; 
		ZVector operator+(const DCplx &c)  const ; 
		ZVector operator+(const double c)  const ; 
		ZVector operator+(const DVector &c)  const ; 
		ZVector operator+(const WVector &c)  const ; 
		ZVector operator-() ; 
		ZVector operator-(const ZVector &c)  const ; 
		ZVector operator-(const DCplx &c)  const ; 
		ZVector operator-(const double c) const  ; 
		ZVector operator-(const DVector &c)  const ; 
		ZVector operator-(const WVector &c)  const ; 
		ZVector operator*(const ZVector &c)  const ; 
		ZVector operator*(const DCplx &c) const  ; 
		ZVector operator*(const double c)  const ; 
		ZVector operator*(const DVector &c) const  ; 
		ZVector operator*(const WVector &c)  const ; 
		ZVector operator*(const DMatrix &d) const;
		ZVector operator*(const ZMatrix &d) const;
		ZVector operator/(const double c)  const ; 
		ZVector operator/(const DCplx &c)  const ;
		ZVector operator<<(int n) const;
		ZVector operator>>(int n) const;


		ZVector &operator=(const ZVector &c); 
		ZVector &operator=(const DVector &c);
		ZVector &operator=(const WVector &c);
		ZVector &operator=(int d) ; 


		ZVector& operator +=(const ZVector &zv);
		ZVector& operator +=(const DVector &dv);
		ZVector& operator +=(const WVector &wv);
		ZVector& operator +=(const DCplx &c);
		ZVector& operator +=(double c);

		ZVector& operator -=(const ZVector &zv);
		ZVector& operator -=(const DVector &dv);
		ZVector& operator -=(const WVector &wv);
		ZVector& operator -=(const DCplx &c);
		ZVector& operator -=(double c);

		ZVector& operator<<=(int n) ;
		ZVector& operator>>=(int n) ;

		
		bool operator==(const DVector &c) ; 
		bool operator==(const WVector &c) ; 
		bool operator==(const ZVector &c) ; 
		bool operator!=(const DVector &c) ; 
		bool operator!=(const WVector &c) ; 
		bool operator!=(const ZVector &c) ; 

		DCplx DotProd(const DVector &rhs);
		DCplx DotProd(const WVector &rhs);
		DCplx DotProd(const ZVector &rhs);

		void fixedconv(const ZVector &rhs) const;
		void conv(const ZVector &rhs);
		ZVector conv_nip(const ZVector &rhs) const;
		void fixedconv(const DVector &rhs);
		void conv(const DVector &rhs);

		ZVector corr_nip(const ZVector &b);
		void corr(const ZVector &b);

		ZVector xcorr_nip() const;

		void xcorr();

		ZVector conv_nip(const DVector &rhs) const;
		void insert(const ZVector &rhs, int pos);
		void insert(const DVector &rhs, int pos);

		ZVector copy(int start, int stop) const;
		ZVector rev_copy(int start, int stop) const;
		ZVector prepad(int count) ;
		ZVector postpad(int count);

		DCplx sum() const;
		void reset(const DCplx &c);

		void print() const {
			if(vect==NULL) return;
			printf("\n[");
			for(register int i=0;i<taille;i++) {
					vect[i].print();
			}
			printf(" ]\n");
		}

};

inline ZVector operator+(const double c1, const ZVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return ZVector();
	ZVector n(c2.taille);
	for(register int i = 0; i < n.taille;i++) n.vect[i]=c2.vect[i]+c1;
	return n;
}

inline ZVector operator-(const double c1, const ZVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return ZVector();
	ZVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1-c2.vect[i];
	return n;
}

inline ZVector operator*(const double c1, const ZVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return ZVector();
	ZVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1*c2.vect[i];
	return n;
}

inline ZVector operator+(const DCplx &c1, const ZVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return ZVector();
	ZVector n(c2.taille);
	for(register int i = 0; i < n.taille;i++) n.vect[i]=c2.vect[i]+c1;
	return n;
}

inline ZVector operator-(const DCplx  &c1, const ZVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return ZVector();
	ZVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1-c2.vect[i];
	return n;
}

inline ZVector operator*(const DCplx  &c1, const ZVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return ZVector();
	ZVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1*c2.vect[i];
	return n;
}

inline DVector operator+(const double c1, const DVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return DVector();
	DVector n(c2.taille);
	for(register int i = 0; i < n.taille;i++) n.vect[i]=c2.vect[i]+c1;
	return n;
}

inline DVector operator-(const double c1, const DVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return DVector();
	DVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1-c2.vect[i];
	return n;
}

inline DVector operator*(const double c1, const DVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return DVector();
	DVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1*c2.vect[i];
	return n;
}

inline ZVector operator+(const DCplx &c1, const DVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return ZVector();
	ZVector n(c2.taille);
	for(register int i = 0; i < n.taille;i++) n.vect[i]=c2.vect[i]+c1;
	return n;
}

inline ZVector operator-(const DCplx  &c1, const DVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return ZVector();
	ZVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1-c2.vect[i];
	return n;
}

inline ZVector operator*(const DCplx  &c1, const DVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return ZVector();
	ZVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1*c2.vect[i];
	return n;
}


inline DVector operator+(const double c1, const WVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return DVector();
	DVector n(c2.taille);
	for(register int i = 0; i < n.taille;i++) n.vect[i]=c2.vect[i]+c1;
	return n;
}

inline DVector operator-(const double c1, const WVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return DVector();
	DVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1-c2.vect[i];
	return n;
}

inline DVector operator*(const double c1, const WVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return DVector();
	DVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1*c2.vect[i];
	return n;
}

inline WVector operator+(const int &c1, const WVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return WVector();
	WVector n(c2.taille);
	for(register int i = 0; i < n.taille;i++) n.vect[i]=c2.vect[i]+c1;
	return n;
}

inline WVector operator-(const int &c1, const WVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return WVector();
	WVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1-c2.vect[i];
	return n;
}

inline WVector operator*(const int &c1, const WVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return WVector();
	WVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1*c2.vect[i];
	return n;
}

inline ZVector operator+(const DCplx &c1, const WVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return ZVector();
	ZVector n(c2.taille);
	for(register int i = 0; i < n.taille;i++) n.vect[i]=c2.vect[i]+c1;
	return n;
}

inline ZVector operator-(const DCplx  &c1, const WVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return ZVector();
	ZVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1-c2.vect[i];
	return n;
}

inline ZVector operator*(const DCplx  &c1, const WVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return ZVector();
	ZVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1*c2.vect[i];
	return n;
}


typedef ZVector cmplx_vect;

class DMatrix{
	public:
		int col, line;	
		double ** mat;
		double trace();


		DMatrix(int line, int col, double initial=0.0);
		DMatrix(const DMatrix& m) ;

		DMatrix(const WMatrix& m);
		DMatrix() ; 
		~DMatrix() ;
		void transpose() ;
		DMatrix transpose_nip();

		ZMatrix ztranspose_nip();

		static DMatrix eye(int size) {
			DMatrix e(size,size);
			for(int i=0;i<size;i++){
				e.mat[i][i]=1.0;
			}
			return e;
		}
		static DMatrix antieye(int size) {
			DMatrix e(size,size);
			for(int i=0;i<size;i++){
				e.mat[i][size-1-i]=1.0;
			}
			return e;
		}

		void eye();
		void antieye();

		DMatrix operator-();

		DMatrix operator+(const DMatrix &rhs);
		DMatrix operator-(const DMatrix &rhs);
		DMatrix operator*(const DMatrix &rhs);
		ZMatrix operator+(const ZMatrix &rhs);
		ZMatrix operator-(const ZMatrix &rhs);
		ZMatrix operator*(const ZMatrix &rhs);
		DMatrix operator*(const double rhs);
		ZMatrix operator*(const DCplx &rhs);
		DMatrix operator/(const double rhs);
		ZMatrix operator/(const DCplx &rhs);
		DMatrix &operator=(const DMatrix &rhs);

		DVector operator*(const DVector &d) const;
		DVector operator*(const WVector &d) const;
		ZVector operator*(const ZVector &d) const;

		bool operator==(const DMatrix &rhs);
		bool operator!=(const DMatrix &rhs);
		void print() const {
			if(mat==NULL) return;
			printf("\n[");
			for(register int j=0;j<line;j++) {
				printf("\n [");
				for(register int i=0;i<col;i++) {
					printf(" %.3f ", mat[i][j]);
				}
				printf(" ]");
			}
			printf("\n]");
		}
		void swapcol(int c1, int c2);
		void swapline(int l1, int l2);
		DMatrix inv_nip();
		void inv();
		void reset(double value=0.0);
		double sum();
		DMatrix linesum();
		DMatrix colsum();

};


typedef DMatrix matrice;

class ZMatrix {
	public: 
		DCplx trace();
		int col, line;		
		cmplx **mat;
		ZMatrix (int line, int col);
		ZMatrix (int line, int col, const DCplx &initial);
		ZMatrix(const ZMatrix& m);
		ZMatrix(const DMatrix& mat);
		ZMatrix(); 
		~ZMatrix() ;
		void transpose();
		ZMatrix  transpose_nip() const ;
		ZMatrix  operator+(const ZMatrix &rhs) const ;
		ZMatrix  operator+(const DMatrix &rhs) const ;
		ZMatrix operator+(const DCplx &c) const;
		ZMatrix operator+(double c) const;

		ZMatrix operator-();
		ZMatrix  operator-(const ZMatrix &rhs) const ;
		ZMatrix  operator-(const DMatrix &rhs) const ;
		ZMatrix operator-(const DCplx &c) const;
		ZMatrix operator-(double c) const;
		ZMatrix operator*(const ZMatrix &rhs) const ;
		ZMatrix operator*(const DMatrix &rhs) const ;
		ZMatrix operator*(const double rhs) const ;
		ZMatrix operator*(const DCplx &rhs) const ;

		ZVector operator*(const DVector &d) const;
		ZVector operator*(const WVector &d) const;
		ZVector operator*(const ZVector &d) const;

		ZMatrix operator/(const double rhs) const ;
		ZMatrix operator/(const DCplx &rhs) const ;
		ZMatrix &operator=(const ZMatrix &rhs);
		ZMatrix &operator=(const DMatrix &rhs);
		ZMatrix &operator+=(const ZMatrix &rhs);
		ZMatrix &operator+=(const DMatrix &rhs);
		ZMatrix &operator+=(const WMatrix &rhs);
		ZMatrix &operator+=(const DCplx &rhs);
		ZMatrix &operator+=(double r);
		ZMatrix &operator-=(const ZMatrix &rhs);
		ZMatrix &operator-=(const DMatrix &rhs);
		ZMatrix &operator-=(const WMatrix &rhs);
		ZMatrix &operator-=(const DCplx &rhs);
		ZMatrix &operator-=(double r);

		bool operator==(const ZMatrix &rhs) const ;
		bool operator!=(const ZMatrix &rhs) const ;
		DMatrix real() const ;
		DMatrix imag() const ;
		DMatrix mag() const ;
		DMatrix phase() const ;
		ZMatrix conj() const ;
		ZMatrix hermit() const;
		void swapcol(int c1, int c2);
		void swapline(int l1, int l2);
		ZMatrix inv_nip();
		void inv();
		DCplx sum();
		ZMatrix linesum();
		ZMatrix colsum();

		static ZMatrix eye(int size){
			ZMatrix e(size,size);
			for(int i=0;i<size;i++){
				e.mat[i][i]=1.0;
			}
			return e;
		}

		static ZMatrix Q_DFT(int size) {
			ZMatrix e(size, size);
			double size_1 = 1.0/size;
			double coef = sqrt(size_1);
			for(int i=0;i<size;i++) {
				for(int j=0;j<size;j++) {
					e.mat[i][j]=DCplx::exp2RI(coef, -UTILITIS_M_2_PI*i*j*size_1);
				}
			}
			return e;
		}

		static ZMatrix Q_IDFT(int size) {
			ZMatrix e(size, size);
			double size_1 = 1.0/size;
			double coef = sqrt(size_1);
			for(int i=0;i<size;i++) {
				for(int j=0;j<size;j++) {
					e.mat[i][j]=DCplx::exp2RI(coef, UTILITIS_M_2_PI*i*j*size_1);
				}
			}
			return e;
		}

		static ZMatrix Q_DFT(int xline, int xcol, double coef=1.0) {
			if (xline <1 || xcol < 1) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
			ZMatrix e(xline, xcol);
			double size_1 = 1.0/xline;
			//double coef = sqrt(size_1);
			for(int i=0;i<xline;i++) {
				for(int j=0;j<xcol;j++) {
					e.mat[i][j]=DCplx::exp2RI(coef, -UTILITIS_M_2_PI*i*j*size_1);
				}
			}
			return e;
		}

		static ZMatrix Q_IDFT(int xline, int xcol, double coef=1.0) {
			if (xline <1 || xcol < 1) throw CUtilitisException(__FILE__, __LINE__, INVALID_MATRIX_DIMENSION);
			ZMatrix e(xline, xcol);
			double size_1 = 1.0/xline;
			//double coef = sqrt(size_1);
			for(int i=0;i<xline;i++) {
				for(int j=0;j<xcol;j++) {
					e.mat[i][j]=DCplx::exp2RI(coef, UTILITIS_M_2_PI*i*j*size_1);
				}
			}
			return e;
		}

		static ZMatrix antieye(int size) {
			ZMatrix e(size,size);
			for(int i=0;i<size;i++){
				e.mat[i][size-1-i]=1.0;
			}
			return e;
		}

		void eye();
		void antieye();

		void print() const {
			if(mat==NULL) return;
			printf("\n[");
			for(register int j=0;j<line;j++) {
				printf("\n [");
				for(register int i=0;i<col;i++) {
					mat[i][j].print();
				}
				printf(" ]");
			}
			printf("\n]");
		}
		void reset(DCplx val);
		ZMatrix ZMatrix::mod2() const ;
		ZMatrix ZMatrix::mod2r() const;

};

typedef ZMatrix cmplx_mat;

inline ZMatrix operator*(const DCplx &lhs, const ZMatrix &rhs){
	return rhs*lhs;
}

inline ZMatrix operator*(const double lhs, const ZMatrix &rhs){
	return rhs*lhs;
}


DMatrix operator*(const double lhs, const DMatrix &rhs);
ZMatrix operator*(const DCplx &lhs, const DMatrix &rhs);

DMatrix operator/(const double lhs, const DMatrix &rhs);

typedef WVector long_vector;

class WMatrix{
	public:
		int max();
		int min();
		int col, line;	
		int ** mat;
		int trace();

		WMatrix(int line, int col, int initial=0);
		WMatrix(const WMatrix& m) ;
		WMatrix() ;
		~WMatrix() ;
		void transpose() ;
		WMatrix transpose_nip();
		ZMatrix ztranspose_nip();

		static WMatrix eye(int size) {
			WMatrix e(size,size);
			for(int i=0;i<size;i++){
				e.mat[i][i]=1;
			}
			return e;
		}

		void eye();
		WMatrix operator-();
		WMatrix operator+(const WMatrix &rhs);
		WMatrix operator-(const WMatrix &rhs);
		WMatrix operator*(const WMatrix &rhs);
		DMatrix operator+(const DMatrix &rhs);
		DMatrix operator-(const DMatrix &rhs);
		DMatrix operator*(const DMatrix &rhs);
		ZMatrix operator+(const ZMatrix &rhs);
		ZMatrix operator-(const ZMatrix &rhs);
		ZMatrix operator*(const ZMatrix &rhs);
		WMatrix operator*(const int &rhs);
		DMatrix operator*(const double rhs);
		ZMatrix operator*(const DCplx &rhs);

		DVector operator*(const DVector &d) const;
		WVector operator*(const WVector &d) const;
		ZVector operator*(const ZVector &d) const;

		DMatrix operator/(const double rhs);
		ZMatrix operator/(const DCplx &rhs);
		WMatrix &operator=(const WMatrix &rhs);
		bool operator==(const WMatrix &rhs);
		bool operator!=(const WMatrix &rhs);
		bool operator==(const DMatrix &rhs);
		bool operator!=(const DMatrix &rhs);
		void print() const {
			if(mat==NULL) return;
			printf("\n[");
			for(register int j=0;j<line;j++) {
				printf("\n [");
				for(register int i=0;i<col;i++) {
					printf(" %d ", mat[i][j]);
				}
				printf(" ]");
			}
			printf("\n]");
		}
		void swapcol(int c1, int c2);
		void swapline(int l1, int l2);
		void reset(int val=0);
		int sum();
		WMatrix linesum();
		WMatrix colsum();
	
};


typedef WMatrix long_matrice;

ostream& operator<< ( ostream& os, const DCplx& dt)  ;
ostream& operator<< ( ostream& os, const DVector& dt)  ;
ostream& operator<< ( ostream& os, const WVector& dt)  ;
ostream& operator<< ( ostream& os, const ZVector& dt)  ;
ostream& operator<< ( ostream& os, const DMatrix& dt)  ;
ostream& operator<< ( ostream& os, const WMatrix& dt)  ;
ostream& operator<< ( ostream& os, const ZMatrix& dt)  ;
//ostream& operator<< (ostream& os, const ZMatrix& dt );//  {

istream& operator>>(istream &is, DCplx &dt);
istream &operator>> (istream &is, DVector &v);
istream &operator>> (istream &is, WVector &v);
istream &operator>> (istream &is, ZVector &v);
istream& operator>>(istream& is, DMatrix& dt );
istream& operator>>(istream& is, WMatrix& dt ); 
istream& operator>>(istream& is, ZMatrix& dt ); 


vector       col_mat2vect(matrice a, long col);

void         Err_Message(char Msg[]);

double       n_gaussien( double moy , double var );
void         random_bit(long_vector a);
void         rand_M(long_vector a, long M);
void         QPSK_mod(long_vector data,cmplx_vect x);
void         QPSK_demod(cmplx_vect x, long_vector data_hat);
void         convol (vector x, vector h, vector y);
void         convol (cmplx_vect x, vector h, cmplx_vect y);
void         convol (cmplx_vect x, cmplx_vect h, cmplx_vect y);

cmplx_vect   correlation(cmplx_vect a, cmplx_vect b);
vector       File2Vect (char Fname[]);
long_vector  File2long_vect (char Fname[]);
void         File2long_matrice (char Fname[],long_matrice a);
cmplx_vect   File2Vect_cmplx (char Fname[]);
void         File_Print (char Fname[],long_vector x);
void         File_Print (char Fname[],vector x);
long         File_Print(char * nom,cmplx_vect a);
void         Fill_vect(vector a, double val);
void         zero_pad(vector x, vector y, long rate);
void         zero_pad(cmplx_vect x, cmplx_vect y, long rate);
void         decimate (vector x, vector y, long delai,long rate);
void         decimate (cmplx_vect x, cmplx_vect y, long delai,long rate);
long         indice_max(vector a);
long         indice_min(vector a);
void         shift_right(long_vector a, long decal);
void         shift_right(vector a, long decal);
void         shift_right(cmplx_vect a, long decal);
void         shift_left(long_vector a,long decal);
void         shift_left(vector a,long decal);
void         shift_left(cmplx_vect a,long decal);
void         lin_vect(vector a,double initial, double step);
cmplx_vect   cmplx_exp(vector a);
vector       copy(vector a);
cmplx_vect   copy(cmplx_vect a);
long_vector  copy(long_vector a);
matrice      copy(matrice a);
long_matrice copy(long_matrice a);
void         awgn(vector a, double SNR, double factor);
void         awgn(cmplx_vect a, double SNR, double factor);
vector       BPSK_MOD(long_vector a);
long_vector	 BPSK_DEM(vector a);
long         comparer(long_vector a, long_vector b);

double convertsnrtosigma2(double snr);
double convertsigma25tosnr(double sigma2);
double convertsnrtosigma2(double Eb, double snr);
double convertsigma2tosnr(double Eb, double sigma2);

void gaussj(const DMatrix &a, const DMatrix &b);

void gaussj(const ZMatrix &a, const ZMatrix &b);

void ludcmp(const DMatrix &a, const WVector &indx, double *d);
void lubksb(const DMatrix &a, const WVector &indx, const DVector &b);
DMatrix ludcmp2(const DMatrix &x, const WVector &xi, double *d);

#endif
