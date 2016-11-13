#ifndef _UTILITIS_H_
#define _UTILITIS_H_

#ifndef _NSP_INTEL_

#include <stdio.h>
#include <string.h>
#include <iostream.h>
#include <math.h>

#ifndef M_2_PI 
#define M_2_PI 6.283185307
#endif

#ifndef M_PI
#define M_PI   3.141592653
#endif
#ifndef M_PI_2
#define M_PI_2 1,570796327
#endif

#define min(x,y)  (((x)<(y)) ? (x):(y))
#define max(x,y)  (((x)>(y)) ? (x):(y))
#define puiss2(a) ((a)*(a))


class WVector;
class ZVector;

class WMatrix;
class ZMatrix;
class DVector;

class DCplx {
	public:
		double re, im;
		// cmplx        init_cmplx(double a, double b);
		DCplx(const double re, const double im);
		DCplx();
		DCplx(const DCplx& c);

		//cmplx        conjug(cmplx a);
		inline DCplx conj() const ;
		
		//double       modul2(cmplx a);
		inline double mag() const ;
		inline double mod2() const ;

		inline double phase() const ;

		//cmplx        cmplx_add(cmplx a,cmplx b);
		inline DCplx operator+(const DCplx &c) const;
		inline DCplx operator+(const double &c) const;

		// cmplx        cmplx_sub(cmplx a,cmplx b);
		inline DCplx operator-(const DCplx &c) const;
		inline DCplx operator-(const double &c) const;

		//cmplx        cmplx_mul(cmplx a,cmplx b);
		inline DCplx operator*(const DCplx &c) const;
		inline DCplx operator*(const double &c) const;

		inline DCplx operator/(const double &c) const;
		
		inline DCplx operator/(const DCplx &c) const {
			double t=c.mod2();
			if (t==0) return DCplx();
			//DCplx ct=c.conj();
			return DCplx((*this)*c.conj()/t);
		}
		
		inline DCplx operator-(void) const { return DCplx(-re, -im); }

		DCplx &operator=(const DCplx &c) ;
		DCplx &operator=(const double &c) ;

		bool operator==(const DCplx &c) const;
		bool operator==(const double &c) const;

		bool operator!=(const DCplx &c) const;
		bool operator!=(const double &c) const;

		DCplx &operator+=(const DCplx &c) ;

		DCplx &operator+=(const double &c) ;

		DCplx &operator-=(const DCplx &c) ;
		DCplx &operator-=(const double &c) ;

		DCplx &operator*=(const DCplx &c) ;
		DCplx &operator*=(const double  &c);

		void print() const { printf(" (%.3f, %.3f) ", re,im); }
};


inline DCplx operator+(const double &c1, const DCplx &c2) ;

inline DCplx operator-(const double &c1, const DCplx &c2) ;

inline DCplx operator*(const double &c1, const DCplx &c2) ;

inline DCplx operator/(const double &c1, const DCplx &c2) ;

inline bool operator==(const double &c1, const DCplx &c2);

inline bool operator!=(const double &c1, const DCplx &c2);

class WVector{ 
	public:
        int * vect;
		int taille;

		WVector(const int &size);
		WVector(const int &size, const int &initial) ;
		WVector(const int *v, const int size) ;
		WVector(const WVector &c);
		WVector() ;
		~WVector() ;
		WVector operator+(const WVector &c) ; 
		DVector operator+(const DVector &c) ; 
		ZVector operator+(const ZVector &c) ; 
		WVector operator+(const int &c) ; 
		DVector operator+(const double &c) ; 
		ZVector operator+(const DCplx &c) ; 

		WVector operator-(const WVector &c) ; 
		DVector operator-(const DVector &c) ; 
		ZVector operator-(const ZVector &c) ; 
		WVector operator-(const int &c) ; 
		DVector operator-(const double &c) ; 
		ZVector operator-(const DCplx &c) ; 

		WVector operator*(const WVector &c) ; 
		DVector operator*(const DVector &c) ; 
		ZVector operator*(const ZVector &c) ; 
		WVector operator*(const int &c) ; 
		DVector operator*(const double &c) ; 
		ZVector operator*(const DCplx &c) ; 

		DVector operator/(const double &c) ; 
		ZVector operator/(const DCplx &c) ;

		WVector &operator=(const WVector &c) ; 

		bool operator==(const DVector &c); 
		bool operator==(const WVector &c) ; 
		bool operator==(const ZVector &c) ; 
		bool operator!=(const WVector &c) ; 
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
		WVector copy(int start, int stop);

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
		DVector operator+(const double &c)  const; 
		DVector operator-(const DVector &c)  const; 
		DVector operator-(const WVector &c)  const; 
		ZVector operator-(const ZVector &c)  const; 
		ZVector operator-(const DCplx &c)  const; 
		DVector operator-(const double &c)  const; 
		ZVector operator*(const ZVector &c)  const; 
		DVector operator*(const DVector &c)  const; 
		DVector operator*(const WVector &c)  const; 
		ZVector operator*(const DCplx &c)  const; 
		DVector operator*(const double &c)  const; 
		DVector operator/(const double &c)  const; 
		ZVector operator/(const DCplx &c)  const;
		DVector &operator=(const DVector &c) ; 
		DVector &operator=(const WVector &c) ; 

		bool operator==(const DVector &c)  const; 
		bool operator==(const WVector &c)  const; 
		bool operator==(const ZVector &c)  const; 
		bool operator!=(const DVector &c)  const; 
		bool operator!=(const WVector &c)  const; 
		bool operator!=(const ZVector &c)  const; 

		void fixedconv(const DVector &rhs);
		void conv(const DVector &rhs);
		DVector conv_nip(const DVector &rhs) const;
		ZVector conv_nip(const ZVector &rhs) const;
		DVector copy(int start, int stop) ;
		DVector prepad(int count) ;
		DVector postpad(int count);
		void insert(const DVector &rhs, int pos);


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
		ZVector operator+(const double &c)  const ; 
		ZVector operator+(const DVector &c)  const ; 
		ZVector operator+(const WVector &c)  const ; 
		ZVector operator-(const ZVector &c)  const ; 
		ZVector operator-(const DCplx &c)  const ; 
		ZVector operator-(const double &c) const  ; 
		ZVector operator-(const DVector &c)  const ; 
		ZVector operator-(const WVector &c)  const ; 
		ZVector operator*(const ZVector &c)  const ; 
		ZVector operator*(const DCplx &c) const  ; 
		ZVector operator*(const double &c)  const ; 
		ZVector operator*(const DVector &c) const  ; 
		ZVector operator*(const WVector &c)  const ; 
		ZVector operator/(const double &c)  const ; 
		ZVector operator/(const DCplx &c)  const ;
		ZVector &operator=(const ZVector &c); 
		ZVector &operator=(const DVector &c);
		ZVector &operator=(const WVector &c);

		bool operator==(const DVector &c) ; 
		bool operator==(const WVector &c) ; 
		bool operator==(const ZVector &c) ; 
		bool operator!=(const DVector &c) ; 
		bool operator!=(const WVector &c) ; 
		bool operator!=(const ZVector &c) ; 

		DCplx DotProd(const DVector &rhs);
		DCplx DotProd(const WVector &rhs);
		DCplx DotProd(const ZVector &rhs);

		void fixedconv(const ZVector &rhs);
		void conv(const ZVector &rhs);
		ZVector conv_nip(const ZVector &rhs) const;
		void fixedconv(const DVector &rhs);
		void conv(const DVector &rhs);
		ZVector conv_nip(const DVector &rhs) const;
		void insert(const ZVector &rhs, int pos);
		void insert(const DVector &rhs, int pos);

		ZVector copy(int start, int stop) ;
		ZVector prepad(int count) ;
		ZVector postpad(int count);

		void print() const {
			if(vect==NULL) return;
			printf("\n[");
			for(register int i=0;i<taille;i++) {
					vect[i].print();
			}
			printf(" ]\n");
		}

};

inline ZVector operator+(const double &c1, const ZVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return ZVector();
	ZVector n(c2.taille);
	for(register int i = 0; i < n.taille;i++) n.vect[i]=c2.vect[i]+c1;
	return n;
}

inline ZVector operator-(const double &c1, const ZVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return ZVector();
	ZVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1-c2.vect[i];
	return n;
}

inline ZVector operator*(const double &c1, const ZVector &c2) {
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

// DVector
inline DVector operator+(const double &c1, const DVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return DVector();
	DVector n(c2.taille);
	for(register int i = 0; i < n.taille;i++) n.vect[i]=c2.vect[i]+c1;
	return n;
}

inline DVector operator-(const double &c1, const DVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return DVector();
	DVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1-c2.vect[i];
	return n;
}

inline DVector operator*(const double &c1, const DVector &c2) {
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


// WVector
inline DVector operator+(const double &c1, const WVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return DVector();
	DVector n(c2.taille);
	for(register int i = 0; i < n.taille;i++) n.vect[i]=c2.vect[i]+c1;
	return n;
}

inline DVector operator-(const double &c1, const WVector &c2) {
	if (c2.taille==0 || c2.vect==NULL) return DVector();
	DVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1-c2.vect[i];
	return n;
}

inline DVector operator*(const double &c1, const WVector &c2) {
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

/*
inline ZVector operator+(const double &c1, const ZVector &c2) ;
inline ZVector operator-(const double &c1, const ZVector &c2) ;
inline ZVector operator*(const double &c1, const ZVector &c2) ;
inline ZVector operator+(const DCplx &c1, const ZVector &c2);
inline ZVector operator-(const DCplx  &c1, const ZVector &c2);
inline ZVector operator*(const DCplx  &c1, const ZVector &c2);

// DVector
inline DVector operator+(const double &c1, const DVector &c2) ;
inline DVector operator-(const double &c1, const DVector &c2) ;
inline DVector operator*(const double &c1, const DVector &c2) ;
inline ZVector operator+(const DCplx &c1, const DVector &c2) ;
inline ZVector operator-(const DCplx  &c1, const DVector &c2) ;
inline ZVector operator*(const DCplx  &c1, const DVector &c2) ;

// WVector
inline DVector operator+(const double &c1, const WVector &c2) ;
inline DVector operator-(const double &c1, const WVector &c2) ;
inline DVector operator*(const double &c1, const WVector &c2) ;
inline WVector operator+(const int &c1, const WVector &c2) ;
inline WVector operator-(const int &c1, const WVector &c2) ;
inline WVector operator*(const int &c1, const WVector &c2) ;
inline ZVector operator+(const DCplx &c1, const WVector &c2) ;
inline ZVector operator-(const DCplx  &c1, const WVector &c2) ;
inline ZVector operator*(const DCplx  &c1, const WVector &c2) ;*/


typedef ZVector cmplx_vect;

class DMatrix{
	public:
		int col, line;	
		double ** mat;

		DMatrix(int line, int col, double initial=0.0);
		DMatrix(const DMatrix& m) ;
		DMatrix() ; 
		~DMatrix() ;
		void transpose() ;
		DMatrix transpose_nip();
		ZMatrix ztranspose_nip();
		DMatrix operator+(const DMatrix &rhs);
		DMatrix operator-(const DMatrix &rhs);
		DMatrix operator*(const DMatrix &rhs);
		ZMatrix operator+(const ZMatrix &rhs);
		ZMatrix operator-(const ZMatrix &rhs);
		ZMatrix operator*(const ZMatrix &rhs);
		DMatrix operator*(const double &rhs);
		ZMatrix operator*(const DCplx &rhs);
		DMatrix operator/(const double &rhs);
		ZMatrix operator/(const DCplx &rhs);
		DMatrix &operator=(const DMatrix &rhs);
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
};

typedef DMatrix matrice;

class ZMatrix {
	public: 
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
		ZMatrix  operator-(const ZMatrix &rhs) const ;
		ZMatrix  operator-(const DMatrix &rhs) const ;
		ZMatrix operator*(const ZMatrix &rhs) const ;
		ZMatrix operator*(const DMatrix &rhs) const ;
		ZMatrix operator*(const double &rhs) const ;
		ZMatrix operator*(const DCplx &rhs) const ;
		ZMatrix operator/(const double &rhs) const ;
		ZMatrix operator/(const DCplx &rhs) const ;
		ZMatrix &operator=(const ZMatrix &rhs);
		ZMatrix &operator=(const DMatrix &rhs);
		bool operator==(const ZMatrix &rhs) const ;
		bool operator!=(const ZMatrix &rhs) const ;
		DMatrix real() const ;
		DMatrix imag() const ;
		DMatrix mag() const ;
		DMatrix phase() const ;
		ZMatrix conj() const ;
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

};

typedef ZMatrix cmplx_mat;

inline ZMatrix operator*(const DCplx &lhs, const ZMatrix &rhs);

inline ZMatrix operator*(const double &lhs, const ZMatrix &rhs);

DMatrix operator*(const double &lhs, const DMatrix &rhs);

ZMatrix operator*(const DCplx &lhs, const DMatrix &rhs);

typedef WVector long_vector;

class WMatrix{
	public:
		int col, line;	
		int ** mat;

		WMatrix(int line, int col, int initial=0);
		WMatrix(const WMatrix& m) ;
		WMatrix() ;
		~WMatrix() ;
		void transpose() ;
		WMatrix transpose_nip();
		ZMatrix ztranspose_nip();
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
		DMatrix operator*(const double &rhs);
		ZMatrix operator*(const DCplx &rhs);
		DMatrix operator/(const double &rhs);
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

};

typedef WMatrix long_matrice;

ostream& operator<< ( ostream& os, DCplx& dt);
ostream& operator<< ( ostream& os, DVector& dt);
ostream& operator<< ( ostream& os, WVector& dt);
ostream& operator<< ( ostream& os, ZVector& dt);
ostream& operator<< ( ostream& os, DMatrix& dt);
ostream& operator<< ( ostream& os, WMatrix& dt);
ostream& operator<< ( ostream& os, ZMatrix& dt);



/*cmplx        init_cmplx(double a, double b);
cmplx        conjug(cmplx a);
double       modul2(cmplx a);
cmplx        cmplx_add(cmplx a,cmplx b);
cmplx        cmplx_sub(cmplx a,cmplx b);
cmplx        cmplx_mul(cmplx a,cmplx b);
double       argum(cmplx x);
vector       init_vect(long ,double);
vector       init_vect(long);
long_vector  init_long_vect(long,long);
void         liberer(vector a);
void         liberer(long_vector);
cmplx_vect   init_cmplx_vect(long , cmplx);
cmplx_vect   init_cmplx_vect(long);
void         liberer(cmplx_vect a);
long_matrice init_long_mat(long LIN, long COL, long valeur);
matrice      init_matrice(long LIN, long COL, double valeur);
void         liberer(matrice a);
void         liberer(long_matrice a);*/

vector       col_mat2vect(matrice a, long col);
//cmplx_vect   conjug(cmplx_vect a);
void         Err_Message(char Msg[]);
/*void	     cmplx_s_mul_vect (cmplx_vect a, cmplx_vect b, cmplx c);
void	     cmplx_mul_vect(cmplx_vect a, cmplx_vect b, cmplx_vect c,double );
void	     cmplx_mul_vect(cmplx_vect a, cmplx_vect b, cmplx_vect c);
void         sub_matrice(matrice a, matrice b, matrice c);
void	     print(cmplx a);
void         print(cmplx_vect a);
void         print(vector a);
void         print(long_vector a);
void         print(matrice a);
void         print(long_matrice a);
*/
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

#endif //_NSP_INTEL_

#endif
 