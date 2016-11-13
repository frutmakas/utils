#ifndef _UTILITIS_H_
#define _UTILITIS_H_

#ifndef _NSP_INTEL_

#include <stdio.h>
#include <iostream.h>
#include <math.h>

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
class DCplx;
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
		inline DCplx conj();
		
		//double       modul2(cmplx a);
		inline double mag() ;
		inline double mod2() ;

		inline double phase();

		//cmplx        cmplx_add(cmplx a,cmplx b);
		inline DCplx operator+(const DCplx &c) ;
		inline DCplx operator+(const double &c) ;

		// cmplx        cmplx_sub(cmplx a,cmplx b);
		inline DCplx operator-(const DCplx &c) ;
		inline DCplx operator-(const double &c) ;

		//cmplx        cmplx_mul(cmplx a,cmplx b);
		inline DCplx operator*(const DCplx &c) ;
		inline DCplx operator*(const double &c) ;

		inline DCplx operator/(const double &c) ;
		inline DCplx operator/(DCplx &c) ;
		DCplx &operator=(const DCplx &c) ;
		DCplx &operator=(const double &c) ;

		bool operator==(const DCplx &c) ;
		bool operator==(const double &c) ;

		bool operator!=(const DCplx &c) ;
		bool operator!=(const double &c) ;

		DCplx &operator+=(const DCplx &c) ;

		DCplx &operator+=(const double &c) ;

		DCplx &operator-=(const DCplx &c) ;
		DCplx &operator-=(const double &c) ;

		DCplx &operator*=(const DCplx &c) ;
};


class WVector{ 
	public:
        int * vect;
		int taille;

		WVector(const int &size){
			if (size==0) { vect=NULL; taille=0; return; }
			vect=new int[size];
			taille=size;
			for(register int i=0;i<size;i++) vect[i]=0;
		}
		
		WVector(const int &size, const int &initial) {
			if (size==0) { vect=NULL; taille=0; return; }
			vect=new int[size];
			taille=size;
			for(register int i=0;i<size;i++) vect[i]=initial;
		}

		WVector(const int *v, const int size) {
			if (size==0 || vect==NULL) { vect=NULL; taille=0; return; }
			vect=new int[size];
			taille=size;			
			for(register int i=0;i<size;i++) vect[i]=v[i];
		}

		WVector(const WVector &c){
			if(c.taille==0) {
				taille=0; vect=NULL; 
				return;
			}
			taille=c.taille;
			vect=new int[taille];
			for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
		}
		
		WVector() {
			taille=0; vect=NULL;
		}
		
		~WVector() {
			if (vect!=NULL) delete[] vect;
			taille=0;
		}

		//cmplx        cmplx_add(cmplx a,cmplx b);
		WVector operator+(const WVector &c) { 
			if (c.taille != taille) return WVector();
			WVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]+c.vect[i];
			}
			return n;
		}

		ZVector operator+(const DCplx &c) { 
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]+c;
			}
			return n;
		}

		DVector operator+(const double &c) { 
			DVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]+(double)c;
			}
			return n;
		}

		WVector operator+(const int &c) { 
			WVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]+c;
			}
			return n;
		}

		ZVector operator+(const ZVector &c) { 
			if (c.taille != taille) return ZVector();
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=c.vect[i]+vect[i];
			}
			return n;
		}

		// cmplx        cmplx_sub(cmplx a,cmplx b);
		WVector operator-(const WVector &c) { 
			if (c.taille != taille) return WVector();
			WVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]-c.vect[i];
			}
			return n;
		}

		DVector operator-(const DVector &c) { 
			if (c.taille != taille) return DVector();
			DVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=((double)vect[i])-c.vect[i];
			}
			return n;
		}

		ZVector operator-(const DCplx &c) { 
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]-(double)c;
			}
			return n;
		}

		WVector operator-(const int &c) { 
			WVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]-c;
			}
			return n;
		}

		DVector operator-(const double &c) { 
			DVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]-(double)c;
			}
			return n;
		}

		ZVector operator-(const ZVector &c) { 
			if (c.taille != taille) return ZVector();
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]-c.vect[i];
			}
			return n;
		}

		//cmplx        cmplx_mul(cmplx a,cmplx b);
		ZVector operator*(const ZVector &c) { 
			if (c.taille != taille) return ZVector();
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]*c.vect[i];
			}
			return n;
		}

		ZVector operator*(const DCplx &c) { 
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]*c;
			}
			return n;
		}

		WVector operator*(const double &c) { 
			WVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]*c;
			}
			return n;
		}

		WVector operator*(const WVector &c) { 
			if (c.taille != taille) return WVector();
			WVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]*c.vect[i];
			}
			return n;
		}

		
		DVector operator*(const DVector &c) { 
			if (c.taille != taille) return DVector();
			WVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]*c.vect[i];
			}
			return n;
		}

		DVector operator/(const double &c) { 
			if (c==0.0) return DVector();
			double d=1/c;
			return DVector((*this)*d);
		}
		
		ZVector operator/(const DCplx &c) {
			double t=c.mod2();
			if (t==0) return ZVector();
			DCplx ct=c.conj()/t;
			return ZVector((*this)*ct);
		}

		WVector &operator=(const WVector &c) { 
			if(taille!=0 && vect!=NULL) delete[] vect;
			taille=c.taille;
			vect = new int[taille];
			for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
			return *this;
		}
};

class DVector { 
	public:
		double * vect;
		int taille;

		DVector(const int &size,  double initial=0.0) {
			if (size==0) { vect=NULL; taille=0; return; }
			vect=new double[size];
			taille=size;
			for(register int i=0;i<size;i++) vect[i]=initial;
		}


		DVector(const double *v, const int size) {
			if (size==0 || vect==NULL) { vect=NULL; taille=0; return; }
			vect=new double[size];
			taille=size;			
			for(register int i=0;i<size;i++) vect[i]=v[i];
		}

		DVector(const WVector &c){
			if(c.taille==0) {
				taille=0; vect=NULL; 
				return;
			}
			taille=c.taille;
			vect=new double[taille];
			for(register int i=0;i<size;i++) vect[i]=(double)c.vect[i];
		}
		
		DVector(const DVector &c){
			if(c.taille==0) {
				taille=0; vect=NULL; 
				return;
			}
			taille=c.taille;
			vect=new double[taille];
			for(register int i=0;i<size;i++) vect[i]=c.vect[i];
		}
		
		DVector() {
			taille=0; vect=NULL;
		}
		
		~DVector() {
			if (vect!=NULL) delete[] vect;
			taille=0;
		}

		//cmplx        cmplx_add(cmplx a,cmplx b);
		DVector operator+(const DVector &c) { 
			if (c.taille != taille) return DVector();
			DVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]+c.vect[i];
			}
			return n;
		}

		DVector operator+(const WVector &c) { 
			if (c.taille != taille) return DVector();
			DVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]+(double)c.vect[i];
			}
			return n;
		}

		ZVector operator+(const ZVector &c) { 
			if (c.taille != taille) return ZVector();
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=c.vect[i]+vect[i];
			}
			return n;
		}

		ZVector operator+(const DCplx &c) { 
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]+c;
			}
			return n;
		}

		DVector operator+(const double &c) { 
			DVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]+c;
			}
			return n;
		}


		// cmplx        cmplx_sub(cmplx a,cmplx b);
		DVector operator-(const DVector &c) { 
			if (c.taille != taille) return DVector();
			DVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]-c.vect[i];
			}
			return n;
		}

		DVector operator-(const WVector &c) { 
			if (c.taille != taille) return DVector();
			DVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]-(double)c.vect[i];
			}
			return n;
		}

		ZVector operator-(const ZVector &c) { 
			if (c.taille != taille) return ZVector();
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]-c.vect[i];
			}
			return n;
		}

		ZVector operator-(const DCplx &c) { 
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]-c;
			}
			return n;
		}

		DVector operator-(const double &c) { 
			DVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]-c;
			}
			return n;
		}

		//cmplx        cmplx_mul(cmplx a,cmplx b);
		ZVector operator*(const ZVector &c) { 
			if (c.taille != taille) return ZVector();
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]*c.vect[i];
			}
			return n;
		}

		DVector operator*(const DVector &c) { 
			if (c.taille != taille) return ZVector();
			DVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]*c.vect[i];
			}
			return n;
		}

		DVector operator*(const WVector &c) { 
			if (c.taille != taille) return ZVector();
			DVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]*(double)c.vect[i];
			}
			return n;
		}

		ZVector operator*(const DCplx &c) { 
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]*c;
			}
			return n;
		}

		DVector operator*(const double &c) { 
			DVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]*c;
			}
			return n;
		}

		DVector operator/(const double &c) { 
			if (c==0.0) return DVector();
			double d=1/c;
			return DVector((*this)*d);
		}
		
		ZVector operator/(const DCplx &c) {
			double t=c.mod2();
			if (t==0) return ZVector();
			DCplx ct=c.conj()/t;
			return ZVector((*this)*ct);
		}

		DVector &operator=(const DVector &c) { 
			if(taille!=0 && vect!=NULL) delete[] vect;
			taille=c.taille;
			vect = new double[taille];
			for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
			return *this;
		}

		DVector &operator=(const WVector &c) { 
			if(taille!=0 && vect!=NULL) delete[] vect;
			taille=c.taille;
			vect = new double[taille];
			for(register int i=0;i<taille;i++) vect[i]=(double)c.vect[i];
			return *this;
		}

		bool operator==(const DVector &c) { 
			if(taille!=c.taille) return false;
			if(taille==0 && c.taille==0) return true;
			if(c.vect==NULL && vect!=NULL) || (vect==NULL && c.vect!=NULL)) retun false;
			for(register int i=0;i<taille;i++) if (vect[i]!=c.vect[i]) return false;
			return true;
		}

		bool operator==(const WVector &c) { 
			if(taille!=c.taille) return false;
			if(taille==0 && c.taille==0) return true;
			if(c.vect==NULL && vect!=NULL) || (vect==NULL && c.vect!=NULL)) retun false;
			for(register int i=0;i<taille;i++) if (vect[i]!=(double)c.vect[i]) return false;
			return true;
		}
		
		bool operator!=(const DVector &c) { 
			if(taille!=c.taille) return true;
			if(taille==0 && c.taille==0) return false;
			if(c.vect==NULL && vect!=NULL) || (vect==NULL && c.vect!=NULL)) retun true;
			for(register int i=0;i<taille;i++) if (vect[i]!=c.vect[i]) return true;
			return false;
		}

		bool operator!=(const WVector &c) { 
			if(taille!=c.taille) return true;
			if(taille==0 && c.taille==0) return false;
			if(c.vect==NULL && vect!=NULL) || (vect==NULL && c.vect!=NULL)) retun true;
			for(register int i=0;i<taille;i++) if (vect[i]!=(double)c.vect[i]) return true;
			return false;
		}
};

typedef DVector vector; 

typedef DCplx cmplx;

class ZVector{
	public:
		cmplx * vect;
		int taille;

		ZVector(const double *real, const double* imag, const int &size){
			if (size==0 || real==NULL || imag==NULL) { vect=NULL; taille=0; return; }
			taille=size;
			vect=new DCplx[taille];
			for(register int i=0;i<size;i++) {
				vect[i].re=real;
				vect[i].imag=imag;
			}
		}

		ZVector(const int &size){
			if (size==0) { vect=NULL; taille=0; return; }
			vect=new cmplx[size];
			taille=size;
		}
		
		ZVector(const int &size, const DCplx &initial) {
			if (size==0) { vect=NULL; taille=0; return; }
			vect=new cmplx(initial)[size];
			taille=size;
			for(register int i=0;i<size;i++) vect[i]=initial;
		}

		ZVector(const cmplx *v, const int size) {
			if (size==0) { vect=NULL; taille=0; return; }
			vect=new cmplx[size];
			taille=size;			
			for(register int i=0;i<size;i++) vect[i]=v[i];
		}

		ZVector(const ZVector &c){
			if(c.taille==0) {
				taille=0; vect=NULL; 
				return;
			}
			taille=c.taille;
			vect=new cmplx[taille];
			for(register int i=0;i<size;i++) vect[i]=c.vect[i];
		}

		ZVector(const DVector &c){
			if(c.taille==0) {
				taille=0; vect=NULL; 
				return;
			}
			taille=c.taille;
			vect=new cmplx[taille];
			for(register int i=0;i<size;i++) vect[i].re=c.vect[i];
		}

		ZVector(const WVector &c){
			if(c.taille==0) {
				taille=0; vect=NULL; 
				return;
			}
			taille=c.taille;
			vect=new cmplx[taille];
			for(register int i=0;i<size;i++) vect[i].re=(double)c.vect[i];
		}

		
		ZVector() {
			taille=0; vect=NULL;
		}
		
		~ZVector() {
			if (vect!=NULL) delete[] vect;
			taille=0;
		}

		//cmplx        conjug(cmplx a);
		ZVector conj(){
			ZVector t(taille);
			for(register int i=0;i<taille;i++) {
				t.vect[i]=vect[i].conj();
			}
			return t;
		}

		DVector real() { 
			if(taille==0 || vect==NULL) return NULL;
			DVector n=new DVector(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i].re;
			}
			return n;
		}

		DVector imag() { 
			if(taille==0 || vect==NULL) return NULL;
			DVector n=new DVector(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i].im;
			}
			return n;
		}

		//double       modul2(cmplx a);
		DVector mag() { 
			if(taille==0 || vect==NULL) return NULL;
			DVector n=new DVector(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i].mag();
			}
			return n;
		}

		DVector mod2() { 
			if(taille==0 || vect==NULL) return NULL;
			DVector n=new DVector(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i].mod2();
			}
					
			return n;
		}

		DVector phase() {
			if(taille==0 || vect==NULL) return NULL;
			DVector n=new DVector(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i].phase();
			}
			return n;
		}

		//cmplx        cmplx_add(cmplx a,cmplx b);
		ZVector operator+(const ZVector &c) { 
			if (c.taille != taille) return ZVector();
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]+c.vect[i];
			}
			return n;
		}

		ZVector operator+(const DCplx &c) { 
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]+c;
			}
			return n;
		}

		ZVector operator+(const double &c) { 
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i].re+c;
			}
			return n;
		}

		ZVector operator+(const DVector &c) { 
			if (c.taille != taille) return ZVector();
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i].re+c.vect[i];
			}
			return n;
		}

		ZVector operator+(const WVector &c) { 
			if (c.taille != taille) return ZVector();
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i].re+(double)c.vect[i];
			}
			return n;
		}

		// cmplx        cmplx_sub(cmplx a,cmplx b);
		ZVector operator-(const ZVector &c) { 
			if (c.taille != taille) return ZVector();
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]-c.vect[i];
			}
			return n;
		}

		ZVector operator-(const DCplx &c) { 
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]-c;
			}
			return n;
		}

		ZVector operator-(const double &c) { 
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]-c;
			}
			return n;
		}

		ZVector operator-(const DVector &c) { 
			if (c.taille != taille) return ZVector();
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i].re-c.vect[i];
			}
			return n;
		}

		ZVector operator-(const WVector &c) { 
			if (c.taille != taille) return ZVector();
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i].re-(double)c.vect[i];
			}
			return n;
		}


		//cmplx        cmplx_mul(cmplx a,cmplx b);
		ZVector operator*(const ZVector &c) { 
			if (c.taille != taille) return ZVector();
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]*c.vect[i];
			}
			return n;
		}

		ZVector operator*(const DCplx &c) { 
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]*c;
			}
			return n;
		}

		ZVector operator*(const double &c) { 
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]*c;
			}
			return n;
		}

		ZVector operator*(const DVector &c) { 
			if (c.taille != taille) return ZVector();
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]*c.vect[i];
			}
			return n;
		}

		ZVector operator*(const WVector &c) { 
			if (c.taille != taille) return ZVector();
			ZVector n(taille);
			for(register int i=0;i<taille;i++){
				n.vect[i]=vect[i]*(double)c.vect[i];
			}
			return n;
		}

		ZVector operator/(const double &c) { 
			if (c==0.0) return ZVector();
			double d=1/c;
			return ZVector((*this)*d);
		}
		
		ZVector operator/(const DCplx &c) {
			double t=c.mod2();
			if (t==0) return ZVector();
			DCplx ct=c.conj()/t;
			return ZVector((*this)*ct);
		}

		ZVector &operator=(const ZVector &c) { 
			if(taille!=0 && vect!=NULL) delete[] vect;
			taille=c.taille;
			vect = new cmplx[taille];
			for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
			return *this;
		}

		ZVector &operator=(const DVector &c){
			if(taille!=0 && vect!=NULL) delete[] vect;
			taille=c.taille;
			vect = new cmplx[taille];
			for(register int i=0;i<taille;i++) {
				vect[i].re=c.vect[i];
			}
			return *this;
		}

		ZVector &operator=(const WVector &c){
			if(taille!=0 && vect!=NULL) delete[] vect;
			taille=c.taille;
			vect = new cmplx[taille];
			for(register int i=0;i<taille;i++) {
				vect[i].re=(double)c.vect[i];
			}
			return *this;
		}
};

inline ZVector operator+(const double &c1, const ZVector &c2) {
	return c2+c1;
}

inline ZVector operator-(const double &c1, const ZVector &c2) {
	ZVector n(c2.taille);
	for(register int i=0;i<c2.taille;i++) n.vect[i]=c1-c2.vect[i];
	return n;
}

inline ZVector operator*(const double &c1, const ZVector &c2) {
	return c2*c1;
}

typedef ZVector cmplx_vect;

class DMatrix{
	public:
		int col, line;	
		double ** mat;

		DMatrix(int col, int line, double initial=0.0){
			if (col==0 || line==0) {
				this->col=this->line=0; mat==NULL; return;
			}
			this->col=col;this->line=line;
			mat=new (double*)[col];
			for(register int i =0;i<col;i++) {
				mat[i]=new double[line];
				for(register int j=0;j<line;j++) mat[i][j]=initial;
			}
		}

		DMatrix(const DMatrix& m) {
			if (col==0 || line==0) {
				this->col=this->line=0; mat==NULL; return;
			}
			this->col=col;this->line=line;
			mat=new (double*)[col];
			size_t line_size=sizeof(double)*line;
			for(register int i =0;i<col;i++) {
				mat[i]=new double[line];
				memcpy(mat[i], m.mat[i], line_size);
			}
		}

		DMatrix() { col=line=0;mat=NULL; }
		//void Free();
		~DMatrix() {
			for(register int i=0;i<col;i++) 
				delete[] mat[i];
			delete[] m;
		}

		void transpose() {
			if (mat==NULL ||col==0 ||line==0) return;
			double **tmp = new (double*)[line];
			for(register int i=0; i<line;i++) {
				tmp[i]=new double[col];
				for(register int j=0;j<col;j++)
					tmp[i][j]=mat[j][i];
			}
			for(i=0;i<col;i++) delete[] mat[i];
			delete[] mat;

			int x=col;col=line;line=x;
			mat=tmp;
		}

		DMatrix transpose_nip(){
			if (mat==NULL ||col==0 ||line==0) return;
			DMatrix t; t.col=line;t.line=col;

			t.mat = new (double*)[line];
			for(register int i=0; i<line;i++) {
				t.mat[i]=new double[col];
				for(register int j=0;j<col;j++)
					t.mat[i][j]=mat[j][i];
			}
			return t;
		}

		ZMatrix ztranspose_nip(){
			if (mat==NULL ||col==0 ||line==0) return;
			ZMatrix t; t.col=line;t.line=col;

			t.mat = new (DCplx*)[line];
			for(register int i=0; i<line;i++) {
				t.mat[i]=new DCplx[col];
				for(register int j=0;j<col;j++)
					t.mat[i][j].re=mat[j][i];
			}
			return t;
		}


		DMatrix operator+(const DMatrix &rhs){
			if (col==0 || line==0 || mat==NULL) return DMatrix();
			if (rhs.col!=col || rhs.line!=line) return DMatrix();
			DMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=mat[i][j]+rhs.mat[i][j];
			}
			return n;
		}

		DMatrix operator-(const DMatrix &rhs){
			if (col==0 || line==0 || mat==NULL) return DMatrix();
			if (rhs.col!=col || rhs.line!=line) return DMatrix();
			DMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=mat[i][j]-rhs.mat[i][j];
			}
			return n;
		}

		DMatrix operator*(const DMatrix &rhs){
			if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return DMatrix();
			if (rhs.line!=col) return DMatrix();
			DMatrix n(rhs.col, line);

			for(register int j=0;j<n.line;j++) {
				for(register int i=0;i<n.col;i++) {
					register double sum=0.0;
					for(register k=0;k<col;k++) 
						sum+=mat[k][j]*rhs.mat[i][k];
					n.mat[i][j]=sum;
				}
			}
			return n;
		}

		ZMatrix operator+(const ZMatrix &rhs){
			if (col==0 || line==0 || mat==NULL) return ZMatrix();
			if (rhs.col!=col || rhs.line!=line) return ZMatrix();
			ZMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=mat[i][j]+rhs.mat[i][j];
			}
			return n;
		}

		ZMatrix operator-(const ZMatrix &rhs){
			if (col==0 || line==0 || mat==NULL) return ZMatrix();
			if (rhs.col!=col || rhs.line!=line) return ZMatrix();
			ZMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=mat[i][j]-rhs.mat[i][j];
			}
			return n;
		}

		ZMatrix operator*(const ZMatrix &rhs){
			if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return ZMatrix();
			if (rhs.line!=col) return ZMatrix();
			ZMatrix n(rhs.col, line);

			for(register int j=0;j<n.line;j++) {
				for(register int i=0;i<n.col;i++) {
					register double sum=0.0;
					for(register k=0;k<col;k++) 
						sum+=mat[k][j]*rhs.mat[i][k];
					n.mat[i][j]=sum;
				}
			}
			return n;
		}
		
		DMatrix operator*(const double &rhs){
			if (col==0 || line==0 || mat==NULL) return DMatrix();
			DMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=mat[i][j]*rhs;
			}
			return n;
		}

		ZMatrix operator*(const DCplx &rhs){
			if (col==0 || line==0 || mat==NULL) return ZMatrix();
			ZMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=mat[i][j]*rhs;
			}
			return n;
		}
		
		DMatrix operator/(const double &rhs){
			if(rhs==0.0) return DMatrix();
			double irhs=1/rhs;
			return (*this)*irhs;
		}

		ZMatrix operator/(const DCplx &rhs){
			if(rhs==0.0) return ZMatrix();
			DCplx irhs=1/rhs;
			return (*this)*irhs;
		}

		DMatrix &operator=(const DMatrix &rhs){
			if(mat!=NULL) {
				for(register int i=0;i<col;i++) delete[] mat[i];
				delete[] mat;
			}

			line=rhs.line;col=rhs.col;

			if(if line==0|| col==0 || rhs.mat==NULL) {
				this->mat=NULL. return
			}

			mat=new (double*)[col];
			for(register int i = 0;i<col;i++) {
				mat[i]=new double[line];
				for(register int j=0;j<line;j++)
					mat[i][j]=rhs.mat[i][j];
			}
			return *this;

		}

		bool operator==(const DMatrix &rhs){
			if(rhs.col!=col || rhs.line!=line) return false;
			for(int i=0;i<col;i++){
				for(int j=0;j<line;j++)
					if (rhs.mat[i][j]!=mat[i][j]) return false;
			}
			return true;
		}

		bool operator!=(const DMatrix &rhs){
			if(rhs.col!=col || rhs.line!=line) return true;
			for(int i=0;i<col;i++){
				for(int j=0;j<line;j++)
					if (rhs.mat[i][j]!=mat[i][j]) return true;
			}
			return false;
		}
};

typedef DMatrix matrice;

class ZMatrix {
	public: 
		int col, line;		
		cmplx **mat;
		ZMatrix (int col, int line){
			if (col==0 || line==0) {
				this->col=this->line=0; mat==NULL; return;
			}
			this->col=col;this->line=line;
			mat=new (DCplx*)[col];
			DCplx initial;
			for(register int i =0;i<col;i++) {
				mat[i]=new DCplx[line];
				for(register int j=0;j<line;j++) mat[i][j]=initial;
			}
		}

		ZMatrix (int col, int line, const DCplx &initial){
			if (col==0 || line==0) {
				this->col=this->line=0; mat==NULL; return;
			}
			this->col=col;this->line=line;
			mat=new (DCplx*)[col];
			for(register int i =0;i<col;i++) {
				mat[i]=new DCplx[line];
				for(register int j=0;j<line;j++) mat[i][j]=initial;
			}
		}

		ZMatrix(const ZMatrix& m){
			if (col==0 || line==0) {
				this->col=this->line=0; mat==NULL; return;
			}
			this->col=col;this->line=line;
			mat=new (DCplx*)[col];
			
			for(register int i =0;i<col;i++) {
				mat[i]=new DCplx[line];
				for(register int j=0;j<line;j++) mat[i][j]=m.mat[i][j];
			}
		}

		ZMatrix(const DMatrix& mat){
			if (col==0 || line==0) {
				this->col=this->line=0; mat==NULL; return;
			}
			this->col=col;this->line=line;
			mat=new (DCplx*)[col];
			
			for(register int i =0;i<col;i++) {
				mat[i]=new DCplx[line];
				for(register int j=0;j<line;j++) mat[i][j]=m.mat[i][j];
			}
		}

		ZMatrix() { col=line=0;mat=NULL;}

		~ZMatrix() {
			if(col==0 || mat==0) return;
			if(line==0) {delete[] mat; return;}
			for(register int i =0; i<col;i++) delete[] mat[i];
			delete[] mat;
		}
		//void Free();
		void transpose(){
			if (mat==NULL ||col==0 ||line==0) return;
			DCplx **tmp = new (DCplx *)[line];
			for(register int i=0; i<line;i++) {
				tmp[i]=new DCplx[col];
				for(register int j=0;j<col;j++)
					tmp[i][j]=mat[j][i];
			}
			for(i=0;i<col;i++) delete[] mat[i];
			delete[] mat;

			int x=col;col=line;line=x;
			mat=tmp;
		}

		ZMatrix  transpose_nip(){
			if (mat==NULL ||col==0 ||line==0) return;
			DMatrix t; t.col=line;t.line=col;

			t.mat = new (DCplx*)[line];
			for(register int i=0; i<line;i++) {
				t.mat[i]=new DCplx[col];
				for(register int j=0;j<col;j++)
					t.mat[i][j]=mat[j][i];
			}
			return t;
		}

		ZMatrix  operator+(const ZMatrix &rhs){
			if (col==0 || line==0 || mat==NULL) return ZMatrix();
			if (rhs.col!=col || rhs.line!=line) return ZMatrix();
			ZMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=mat[i][j]+rhs.mat[i][j];
			}
			return n;
		}

		ZMatrix  operator+(const DMatrix &rhs){
			if (col==0 || line==0 || mat==NULL) return ZMatrix();
			if (rhs.col!=col || rhs.line!=line) return ZMatrix();
			ZMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=mat[i][j]+rhs.mat[i][j];
			}
			return n;
		}

		ZMatrix  operator-(const ZMatrix &rhs){
			if (col==0 || line==0 || mat==NULL) return ZMatrix();
			if (rhs.col!=col || rhs.line!=line) return ZMatrix();
			ZMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=mat[i][j]-rhs.mat[i][j];
			}
			return n;
		}

		ZMatrix  operator-(const DMatrix &rhs){
			if (col==0 || line==0 || mat==NULL) return ZMatrix();
			if (rhs.col!=col || rhs.line!=line) return ZMatrix();
			ZMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=mat[i][j]-rhs.mat[i][j];
			}
			return n;
		}

		ZMatrix operator*(const ZMatrix &rhs){
			if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return ZMatrix();
			if (rhs.line!=col) return ZMatrix();
			ZMatrix n(rhs.col, line);

			for(register int j=0;j<n.line;j++) {
				for(register int i=0;i<n.col;i++) {
					register double sum=0.0;
					for(register k=0;k<col;k++) 
						sum+=mat[k][j]*rhs.mat[i][k];
					n.mat[i][j]=sum;
				}
			}
			return n;
		}

		ZMatrix operator*(const DMatrix &rhs){
			if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return ZMatrix();
			if (rhs.line!=col) return ZMatrix();
			ZMatrix n(rhs.col, line);

			for(register int j=0;j<n.line;j++) {
				for(register int i=0;i<n.col;i++) {
					register double sum=0.0;
					for(register k=0;k<col;k++) 
						sum+=mat[k][j]*rhs.mat[i][k];
					n.mat[i][j]=sum;
				}
			}
			return n;
		}

		ZMatrix operator*(const double &rhs){
			if (col==0 || line==0 || mat==NULL) return ZMatrix();
			ZMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=mat[i][j]*rhs;
			}
			return n;
		}

		ZMatrix operator*(const DCplx &rhs){
			if (col==0 || line==0 || mat==NULL) return ZMatrix();
			ZMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=mat[i][j]*rhs;
			}
			return n;
		}
		

		ZMatrix operator/(const double &rhs){
			if(rhs==0.0) return ZMatrix();
			double irhs=1/rhs;
			return (*this)*irhs;
		}

		ZMatrix operator/(const DCplx &rhs){
			if(rhs==0.0) return ZMatrix();
			DCplx irhs=1/rhs;
			return (*this)*irhs;
		}

		ZMatrix &operator=(const ZMatrix &rhs){
			if(mat!=NULL) {
				for(register int i=0;i<col;i++) delete[] mat[i];
				delete[] mat;
			}

			line=rhs.line;col=rhs.col;

			if(if line==0|| col==0 || rhs.mat==NULL) {
				this->mat=NULL. return
			}

			mat=new (double*)[col];
			for(register int i = 0;i<col;i++) {
				mat[i]=new double[line];
				for(register int j=0;j<line;j++)
					mat[i][j]=rhs.mat[i][j];
			}
			return *this;
		}

		ZMatrix &operator=(const DMatrix &rhs){
			if(mat!=NULL) {
				for(register int i=0;i<col;i++) delete[] mat[i];
				delete[] mat;
			}

			line=rhs.line;col=rhs.col;

			if(if line==0|| col==0 || rhs.mat==NULL) {
				this->mat=NULL. return
			}

			mat=new (double*)[col];
			for(register int i = 0;i<col;i++) {
				mat[i]=new double[line];
				for(register int j=0;j<line;j++)
					mat[i][j]=rhs.mat[i][j];
			}
			return *this;
		}

		bool operator==(const ZMatrix &rhs){
			if(rhs.col!=col || rhs.line!=line) return false;
			for(int i=0;i<col;i++){
				for(int j=0;j<line;j++)
					if (rhs.mat[i][j]!=mat[i][j]) return false;
			}
			return true;
		}

		bool operator!=(const ZMatrix &rhs){
			if(rhs.col!=col || rhs.line!=line) return true;
			for(int i=0;i<col;i++){
				for(int j=0;j<line;j++)
					if (rhs.mat[i][j]!=mat[i][j]) return true;
			}
			return false;
		}

		DMatrix real(){
			if(col==0 || line==0 || mat==NULL) return DMatrix();
			DMatrix n(col, line);
			n.mat=new (double*)[col];
			for(register int i = 0;i<col;i++) {
				mat[i]=new double[line];
				for(register int j=0;j<line;j++)
					mat[i][j]=rhs.mat[i][j].re;
			}
			return n;
		}

		DMatrix imag(){
			if(col==0 || line==0 || mat==NULL) return DMatrix();
			DMatrix n(col, line);
			n.mat=new (double*)[col];
			for(register int i = 0;i<col;i++) {
				mat[i]=new double[line];
				for(register int j=0;j<line;j++)
					mat[i][j]=rhs.mat[i][j].im;
			}
			return n;
		}

		DMatrix mag(){
			if(col==0 || line==0 || mat==NULL) return DMatrix();
			DMatrix n(col, line);
			n.mat=new (double*)[col];
			for(register int i = 0;i<col;i++) {
				mat[i]=new double[line];
				for(register int j=0;j<line;j++)
					mat[i][j]=rhs.mat[i][j].mag();
			}
			return n;
		}
		DMatrix phase()
			if(col==0 || line==0 || mat==NULL) return DMatrix();
			DMatrix n(col, line);
			n.mat=new (double*)[col];
			for(register int i = 0;i<col;i++) {
				mat[i]=new double[line];
				for(register int j=0;j<line;j++)
					mat[i][j]=rhs.mat[i][j].phase();
			}
			return n;
		}

		ZMatrix conj(){
			if(col==0 || line==0 || mat==NULL) return ZMatrix();
			ZMatrix n(col, line);
			n.mat=new (DCplx*)[col];
			for(register int i = 0;i<col;i++) {
				mat[i]=new DCplx[line];
				for(register int j=0;j<line;j++)
					mat[i][j]=rhs.mat[i][j].conj();
			}
			return n;
		}
};

typedef ZMatrix cmplx_mat;

inline ZMatrix operator*(const DCplx &lhs, const ZMatrix &rhs){
	return rhs*lhs;
}

inline ZMatrix operator*(const double &lhs, const ZMatrix &rhs){
	return rhs*lhs;
}

DMatrix operator*(const double &lhs, const DMatrix &rhs){
	return rhs*lhs;
}

ZMatrix operator*(const DCplx &lhs, const DMatrix &rhs){
	return rhs*lhs;
}


typedef WVector long_vector;

class WMatrix{
	public:
		int col, line;	
		int ** mat;

		WMatrix(int col, int line, double initial=0.0){
			if (col==0 || line==0) {
				this->col=this->line=0; mat==NULL; return;
			}
			this->col=col;this->line=line;
			mat=new (int*)[col];
			for(register int i =0;i<col;i++) {
				mat[i]=new int[line];
				for(register int j=0;j<line;j++) mat[i][j]=initial;
			}
		}

		WMatrix(const WMatrix& m) {
			if (col==0 || line==0) {
				this->col=this->line=0; mat==NULL; return;
			}
			this->col=col;this->line=line;
			mat=new (int*)[col];
			size_t line_size=sizeof(double)*line;
			for(register int i =0;i<col;i++) {
				mat[i]=new int[line];
				memcpy(mat[i], m.mat[i], line_size);
			}
		}

		WMatrix() { col=line=0;mat=NULL; }
		//void Free();
		~WMatrix() {
			for(register int i=0;i<col;i++) 
				delete[] mat[i];
			delete[] m;
		}

		void transpose() {
			if (mat==NULL ||col==0 ||line==0) return;
			double **tmp = new (int*)[line];
			for(register int i=0; i<line;i++) {
				tmp[i]=new int[col];
				for(register int j=0;j<col;j++)
					tmp[i][j]=mat[j][i];
			}
			for(i=0;i<col;i++) delete[] mat[i];
			delete[] mat;

			int x=col;col=line;line=x;
			mat=tmp;
		}

		WMatrix transpose_nip(){
			if (mat==NULL ||col==0 ||line==0) return;
			WMatrix t; t.col=line;t.line=col;

			t.mat = new (int*)[line];
			for(register int i=0; i<line;i++) {
				t.mat[i]=new int[col];
				for(register int j=0;j<col;j++)
					t.mat[i][j]=mat[j][i];
			}
			return t;
		}

		ZMatrix ztranspose_nip(){
			if (mat==NULL ||col==0 ||line==0) return;
			ZMatrix t; t.col=line;t.line=col;

			t.mat = new (DCplx*)[line];
			for(register int i=0; i<line;i++) {
				t.mat[i]=new DCplx[col];
				for(register int j=0;j<col;j++)
					t.mat[i][j].re=(double)mat[j][i];
			}
			return t;
		}


		WMatrix operator+(const WMatrix &rhs){
			if (col==0 || line==0 || mat==NULL) return WMatrix();
			if (rhs.col!=col || rhs.line!=line) return WMatrix();
			WMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=mat[i][j]+rhs.mat[i][j];
			}
			return n;
		}

		WMatrix operator-(const WMatrix &rhs){
			if (col==0 || line==0 || mat==NULL) return WMatrix();
			if (rhs.col!=col || rhs.line!=line) return WMatrix();
			WMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=mat[i][j]-rhs.mat[i][j];
			}
			return n;
		}

		WMatrix operator*(const WMatrix &rhs){
			if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return WMatrix();
			if (rhs.line!=col) return WMatrix();
			WMatrix n(rhs.col, line);

			for(register int j=0;j<n.line;j++) {
				for(register int i=0;i<n.col;i++) {
					register double sum=0.0;
					for(register k=0;k<col;k++) 
						sum+=mat[k][j]*rhs.mat[i][k];
					n.mat[i][j]=sum;
				}
			}
			return n;
		}

		DMatrix operator+(const DMatrix &rhs){
			if (col==0 || line==0 || mat==NULL) return DMatrix();
			if (rhs.col!=col || rhs.line!=line) return DMatrix();
			DMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=(double)mat[i][j]+rhs.mat[i][j];
			}
			return n;
		}

		DMatrix operator-(const DMatrix &rhs){
			if (col==0 || line==0 || mat==NULL) return DMatrix();
			if (rhs.col!=col || rhs.line!=line) return DMatrix();
			DMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=(double)mat[i][j]-rhs.mat[i][j];
			}
			return n;
		}

		DMatrix operator*(const DMatrix &rhs){
			if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return DMatrix();
			if (rhs.line!=col) return DMatrix();
			DMatrix n(rhs.col, line);

			for(register int j=0;j<n.line;j++) {
				for(register int i=0;i<n.col;i++) {
					register double sum=0.0;
					for(register k=0;k<col;k++) 
						sum+=(double)mat[k][j]*rhs.mat[i][k];
					n.mat[i][j]=sum;
				}
			}
			return n;
		}

		ZMatrix operator+(const ZMatrix &rhs){
			if (col==0 || line==0 || mat==NULL) return ZMatrix();
			if (rhs.col!=col || rhs.line!=line) return ZMatrix();
			ZMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j].re=mat[i][j].re+rhs.mat[i][j].re;
			}
			return n;
		}

		ZMatrix operator-(const ZMatrix &rhs){
			if (col==0 || line==0 || mat==NULL) return ZMatrix();
			if (rhs.col!=col || rhs.line!=line) return ZMatrix();
			ZMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j].re=mat[i][j].re-rhs.mat[i][j].re;
			}
			return n;
		}

		ZMatrix operator*(const ZMatrix &rhs){
			if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return ZMatrix();
			if (rhs.line!=col) return ZMatrix();
			ZMatrix n(rhs.col, line);

			for(register int j=0;j<n.line;j++) {
				for(register int i=0;i<n.col;i++) {
					register double sum=0.0;
					for(register k=0;k<col;k++) 
						sum+=mat[k][j]*rhs.mat[i][k];
					n.mat[i][j]=sum;
				}
			}
			return n;
		}
		
		WMatrix operator*(const int &rhs){
			if (col==0 || line==0 || mat==NULL) return WMatrix();
			WMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=mat[i][j]*rhs;
			}
			return n;
		}

		DMatrix operator*(const double &rhs){
			if (col==0 || line==0 || mat==NULL) return DMatrix();
			DMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=double(mat[i][j])*rhs;
			}
			return n;
		}

		ZMatrix operator*(const DCplx &rhs){
			if (col==0 || line==0 || mat==NULL) return ZMatrix();
			ZMatrix n(col, line);
			for(register int i=0;i<col;i++) {
				for(register int j=0;j<line;j++)
					n.mat[i][j]=mat[i][j]*rhs;
			}
			return n;
		}
		
		DMatrix operator/(const double &rhs){
			if(rhs==0.0) return DMatrix();
			double irhs=1/rhs;
			return (*this)*irhs;
		}

		ZMatrix operator/(const DCplx &rhs){
			if(rhs==0.0) return ZMatrix();
			DCplx irhs=1/rhs;
			return (*this)*irhs;
		}

		WMatrix &operator=(const WMatrix &rhs){
			if(mat!=NULL) {
				for(register int i=0;i<col;i++) delete[] mat[i];
				delete[] mat;
			}

			line=rhs.line;col=rhs.col;

			if(if line==0|| col==0 || rhs.mat==NULL) {
				this->mat=NULL. return
			}

			mat=new (double*)[col];
			for(register int i = 0;i<col;i++) {
				mat[i]=new double[line];
				for(register int j=0;j<line;j++)
					mat[i][j]=rhs.mat[i][j];
			}
			return *this;

		}

		bool operator==(const WMatrix &rhs){
			if(rhs.col!=col || rhs.line!=line) return false;
			for(int i=0;i<col;i++){
				for(int j=0;j<line;j++)
					if (rhs.mat[i][j]!=mat[i][j]) return false;
			}
			return true;
		}

		bool operator!=(const WMatrix &rhs){
			if(rhs.col!=col || rhs.line!=line) return true;
			for(int i=0;i<col;i++){
				for(int j=0;j<line;j++)
					if (rhs.mat[i][j]!=mat[i][j]) return true;
			}
			return false;
		}

		bool operator==(const DMatrix &rhs){
			if(rhs.col!=col || rhs.line!=line) return false;
			for(int i=0;i<col;i++){
				for(int j=0;j<line;j++)
					if (rhs.mat[i][j]!=mat[i][j]) return false;
			}
			return true;
		}

		bool operator!=(const DMatrix &rhs){
			if(rhs.col!=col || rhs.line!=line) return true;
			for(int i=0;i<col;i++){
				for(int j=0;j<line;j++)
					if (rhs.mat[i][j]!=mat[i][j]) return true;
			}
			return false;
		}

};

typedef WMatrix long_matrice;



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
 