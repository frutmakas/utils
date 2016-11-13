#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utilitis.h"

// cmplx        init_cmplx(double a, double b);
DCplx::DCplx(const double re, const double im){this->re=re; this->im=im;}
DCplx::DCplx(){this->re=0.0; this->im=0.0;}
DCplx::DCplx(const DCplx& c){re=c.re; im=c.im;}

//cmplx        conjug(cmplx a);
inline DCplx DCplx::conj() const {return DCplx(re,-im);}

//double       modul2(cmplx a);
inline double DCplx::mag() const { return sqrt(re*re+im*im); }
inline double DCplx::mod2() const { return re*re+im*im; }

inline double DCplx::phase() const {
	if (re==0 && im==0) {
		return 0.0;
	}
	if (re==0 && im!=0) { 
		return (im>0 ) ?  M_PI_2 : -M_PI_2; 
	}
	//double at=atan(im/re);
	return atan(im/re);

}

//cmplx        cmplx_add(cmplx a,cmplx b);
inline DCplx DCplx::operator+(const DCplx &c) const{ return DCplx(re+c.re, im+c.im); }
inline DCplx DCplx::operator+(const double &c)const { return DCplx(re+c, im); }

// cmplx        cmplx_sub(cmplx a,cmplx b);
inline DCplx DCplx::operator-(const DCplx &c) const{ return DCplx(re-c.re, im-c.im); }
inline DCplx DCplx::operator-(const double &c) const{ return DCplx(re-c, im); }

//cmplx        cmplx_mul(cmplx a,cmplx b);
inline DCplx DCplx::operator*(const DCplx &c) const{ return DCplx(re*c.re-c.im*im, im*c.re+re*c.im); }
inline DCplx DCplx::operator*(const double &c) const{ return DCplx(re*c, c*im); }

inline DCplx DCplx::operator/(const double &c) const{ if (c!=0.) return DCplx(re/c, im/c); return DCplx();}

/*inline DCplx DCplx::operator/(const DCplx &c) const{
	double t=c.mod2();
	if (t==0) return DCplx();
	//DCplx ct=c.conj();
	return DCplx((*this)*c.conj()/t);
}*/

DCplx &DCplx::operator=(const DCplx &c) { re=c.re; im=c.im ; return *this; }
DCplx &DCplx::operator=(const double &c) { re=c;im=0.0; return *this; }

bool DCplx::operator==(const DCplx &c) const{ return (re==c.re)&&(im==c.im); }
bool DCplx::operator==(const double &c) const{ return (re==c)&&(im==0.0); }

bool DCplx::operator!=(const DCplx &c) const{ return !((re==c.re)&&(im==c.im)); }
bool DCplx::operator!=(const double &c) const{ return !((re==c)&&(im==0.0)); }

DCplx &DCplx::operator+=(const DCplx &c) {
	re+=c.re;
	im+=c.im;
	return (*this); 
}

DCplx &DCplx::operator+=(const double &c) {
	re+=c;
	return (*this); 
}

DCplx &DCplx::operator-=(const DCplx &c) {
	re-=c.re;
	im-=c.im;
	return (*this); 
}

DCplx &DCplx::operator-=(const double &c) {
	re-=c;
	return (*this); 
}

DCplx &DCplx::operator*=(const DCplx &c) {
	DCplx x((*this)*c);
	re=x.re;
	im=x.im;
	return (*this); 
}

DCplx &DCplx::operator*=(const double  &c) {
	re*=c;
	im*=c;
	return (*this); 
}



inline DCplx operator+(const double &c1, const DCplx &c2) {
	return DCplx(c2.re+c1, c2.im);
}

inline DCplx operator-(const double &c1, const DCplx &c2) {
	return DCplx(c1-c2.re, c2.im);
}

inline DCplx operator*(const double &c1, const DCplx &c2) {
	return DCplx(c2.re*c1, c2.im*c1);
}

inline DCplx operator/(const double &c1, const DCplx &c2) {
	return c2.conj()*c1/c2.mod2();
}


inline bool operator==(const double &c1, const DCplx &c2) {
	return c1==c2.re && c2.im==0.0;
}

inline bool operator!=(const double &c1, const DCplx &c2) {
	return !(c1==c2.re && c2.im==0.0);
}

WVector::WVector(const int &size){
	if (size==0) { vect=NULL; taille=0; return; }
	vect=new int[size];
	taille=size;
	for(register int i=0;i<size;i++) vect[i]=0;
}

WVector::WVector(const int &size, const int &initial) {
	if (size==0) { vect=NULL; taille=0; return; }
	vect=new int[size];
	taille=size;
	for(register int i=0;i<size;i++) vect[i]=initial;
}

WVector::WVector(const int *v, const int size) {
	if (size==0 || vect==NULL) { vect=NULL; taille=0; return; }
	vect=new int[size];
	taille=size;			
	for(register int i=0;i<size;i++) vect[i]=v[i];
}

WVector::WVector(const WVector &c){
	if(c.taille==0) {
		taille=0; vect=NULL; 
		return;
	}
	taille=c.taille;
	vect=new int[taille];
	for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
}

WVector::WVector() {
	taille=0; vect=NULL;
}

WVector::~WVector() {
	if (vect!=NULL) delete[] vect;
	vect=NULL;taille=0;
}

//cmplx        cmplx_add(cmplx a,cmplx b);
WVector WVector::operator+(const WVector &c) { 
	if (c.taille != taille) return WVector();
	WVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c.vect[i];
	}
	return n;
}

DVector WVector::operator+(const DVector &c) { 
	if (c.taille != taille) return DVector();
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c.vect[i];
	}
	return n;
}

ZVector WVector::operator+(const DCplx &c) { 
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c;
	}
	return n;
}

DVector WVector::operator+(const double &c) { 
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+(double)c;
	}
	return n;
}

WVector WVector::operator+(const int &c) { 
	WVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c;
	}
	return n;
}

ZVector WVector::operator+(const ZVector &c) { 
	if (c.taille != taille) return ZVector();
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=c.vect[i]+vect[i];
	}
	return n;
}

// cmplx        cmplx_sub(cmplx a,cmplx b);
WVector WVector::operator-(const WVector &c) { 
	if (c.taille != taille) return WVector();
	WVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c.vect[i];
	}
	return n;
}

DVector WVector::operator-(const DVector &c) { 
	if (c.taille != taille) return DVector();
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=((double)vect[i])-c.vect[i];
	}
	return n;
}

ZVector WVector::operator-(const DCplx &c) { 
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=(double)vect[i]-c;
	}
	return n;
}

WVector WVector::operator-(const int &c) { 
	WVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c;
	}
	return n;
}

DVector WVector::operator-(const double &c) { 
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-(double)c;
	}
	return n;
}

ZVector WVector::operator-(const ZVector &c) { 
	if (c.taille != taille) return ZVector();
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c.vect[i];
	}
	return n;
}

//cmplx        cmplx_mul(cmplx a,cmplx b);
ZVector WVector::operator*(const ZVector &c) { 
	if (c.taille != taille) return ZVector();
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c.vect[i];
	}
	return n;
}

ZVector WVector::operator*(const DCplx &c) { 
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c;
	}
	return n;
}


WVector WVector::operator*(const int &c) { 
	WVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c;
	}
	return n;
}

DVector WVector::operator*(const double &c) { 
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c;
	}
	return n;
}


WVector WVector::operator*(const WVector &c) { 
	if (c.taille != taille) return WVector();
	WVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c.vect[i];
	}
	return n;
}


DVector WVector::operator*(const DVector &c) { 
	if (c.taille != taille) return DVector();
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c.vect[i];
	}
	return n;
}

DVector WVector::operator/(const double &c) { 
	if (c==0.0) return DVector();
	double d=1/c;
	return DVector((*this)*d);
}

ZVector WVector::operator/(const DCplx &c) {
	double t=c.mod2();
	if (t==0) return ZVector();
	DCplx ct=c.conj()/t;
	return ZVector((*this)*ct);
}

WVector &WVector::operator=(const WVector &c) { 
	if(taille!=0 && vect!=NULL) delete[] vect;
	taille=c.taille;
	vect = new int[taille];
	for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
	return *this;
}

bool WVector::operator==(const DVector &c){
	if(taille!=c.taille) return false;
	register int i=0;
	while(i<taille && vect[i]==c.vect[i]) i++;
	return (i==taille);
}

bool WVector::operator==(const WVector &c){
	if(taille!=c.taille) return false;
	register int i=0;
	while(i<taille && vect[i]==c.vect[i]) i++;
	return (i==taille);
}

bool WVector::operator==(const ZVector &c) {
	if(taille!=c.taille) return false;
	register int i=0;
	while(i<taille && vect[i]==c.vect[i]) i++;
	return (i==taille);
}

bool WVector::operator!=(const DVector &c) {
	if(taille!=c.taille) return true;
	register int i=0;
	while(vect[i]!=c.vect[i] && i<taille) i++;
	return !(i==taille);
}

bool WVector::operator!=(const WVector &c) {
	if(taille!=c.taille) return true;
	register int i=0;
	while(vect[i]!=c.vect[i] && i<taille) i++;
	return !(i==taille);
}

bool WVector::operator!=(const ZVector &c) {
	if(taille!=c.taille) return true;
	register int i=0;
	while(vect[i]!=c.vect[i] && i<taille) i++;
	return !(i==taille);
}


double WVector::DotProd(const DVector &rhs){
	if (taille!=rhs.taille) return 0.0;
	register double sum=0.0;
	for(register int i=0;i<taille;i++) {
		sum+=vect[i]*rhs.vect[i];
	}
	return sum;
}

int WVector::DotProd(const WVector &rhs){
	if (taille!=rhs.taille) return 0;
	register int sum=0;
	for(register int i=0;i<taille;i++) {
		sum+=vect[i]*rhs.vect[i];
	}
	return sum;
}

DCplx  WVector::DotProd(const ZVector &rhs){
	if (taille!=rhs.taille) return DCplx();
	register DCplx sum(0.0,0.0);
	for(register int i=0;i<taille;i++) {
		sum+=rhs.vect[i]*vect[i];
	}
	return sum;
}

WVector WVector::copy(int start, int stop) {
	stop = start > stop ? start : stop;
	start = start > stop ? stop : start;
	if (stop>taille) stop=taille; 
	if (start<0) start=0;
	if (!(stop-start)) return WVector();
	WVector n(stop-start);
	for(register int i=start, j=0;i<stop;i++,j++)
		n.vect[j]=vect[i];
	return n;
}

DVector::DVector(const int &size,  double initial) {
	if (size==0) { vect=NULL; taille=0; return; }
	vect=new double[size];
	taille=size;
	for(register int i=0;i<size;i++) vect[i]=initial;
}


DVector::DVector(const double *v, const int size) {
	if (size==0 || vect==NULL) { vect=NULL; taille=0; return; }
	vect=new double[size];
	taille=size;			
	for(register int i=0;i<size;i++) vect[i]=v[i];
}

DVector::DVector(const WVector &c){
	if(c.taille==0) {
		taille=0; vect=NULL; 
		return;
	}
	taille=c.taille;
	vect=new double[taille];
	for(register int i=0;i<taille;i++) vect[i]=(double)c.vect[i];
}

DVector::DVector(const DVector &c){
	if(c.taille==0) {
		taille=0; vect=NULL; 
		return;
	}
	taille=c.taille;
	vect=new double[taille];
	for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
}

DVector::DVector() {
	taille=0; vect=NULL;
}

DVector::~DVector() {
	if (vect!=NULL) delete[] vect;
	vect=NULL; taille=0;
}

//cmplx        cmplx_add(cmplx a,cmplx b);
DVector DVector::operator+(const DVector &c) const{ 
	if (c.taille != taille) return DVector();
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c.vect[i];
	}
	return n;
}

DVector DVector::operator+(const WVector &c)  const{ 
	if (c.taille != taille) return DVector();
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+(double)c.vect[i];
	}
	return n;
}

ZVector DVector::operator+(const ZVector &c)  const{ 
	if (c.taille != taille) return ZVector();
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=c.vect[i]+vect[i];
	}
	return n;
}

ZVector DVector::operator+(const DCplx &c)  const{ 
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c;
	}
	return n;
}

DVector DVector::operator+(const double &c)  const{ 
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c;
	}
	return n;
}


// cmplx        cmplx_sub(cmplx a,cmplx b);
DVector DVector::operator-(const DVector &c)  const{ 
	if (c.taille != taille) return DVector();
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c.vect[i];
	}
	return n;
}

DVector DVector::operator-(const WVector &c)  const{ 
	if (c.taille != taille) return DVector();
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-(double)c.vect[i];
	}
	return n;
}

ZVector DVector::operator-(const ZVector &c)  const{ 
	if (c.taille != taille) return ZVector();
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c.vect[i];
	}
	return n;
}

ZVector DVector::operator-(const DCplx &c)  const{ 
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c;
	}
	return n;
}

DVector DVector::operator-(const double &c)  const{ 
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c;
	}
	return n;
}

//cmplx        cmplx_mul(cmplx a,cmplx b);
ZVector DVector::operator*(const ZVector &c)  const{ 
	if (c.taille != taille) return ZVector();
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c.vect[i];
	}
	return n;
}

DVector DVector::operator*(const DVector &c)  const{ 
	if (c.taille != taille) return DVector();
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c.vect[i];
	}
	return n;
}

DVector DVector::operator*(const WVector &c)  const{ 
	if (c.taille != taille) return DVector();
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*(double)c.vect[i];
	}
	return n;
}

ZVector DVector::operator*(const DCplx &c)  const{ 
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c;
	}
	return n;
}

DVector DVector::operator*(const double &c)  const{ 
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c;
	}
	return n;
}

DVector DVector::operator/(const double &c)  const{ 
	if (c==0.0) return DVector();
	double d=1/c;
	return DVector((*this)*d);
}

ZVector DVector::operator/(const DCplx &c)  const{
	double t=c.mod2();
	if (t==0) return ZVector();
	DCplx ct=c.conj()/t;
	return ZVector((*this)*ct);
}

DVector &DVector::operator=(const DVector &c){ 
	if(taille!=0 && vect!=NULL) delete[] vect;
	taille=c.taille;
	vect = new double[taille];
	for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
	return *this;
}

DVector &DVector::operator=(const WVector &c){ 
	if(taille!=0 && vect!=NULL) delete[] vect;
	taille=c.taille;
	vect = new double[taille];
	for(register int i=0;i<taille;i++) vect[i]=(double)c.vect[i];
	return *this;
}

bool DVector::operator==(const DVector &c)  const{ 
	if(taille!=c.taille) return false;
	if(taille==0 && c.taille==0) return true;
	if((c.vect==NULL && vect!=NULL) || (vect==NULL && c.vect!=NULL)) return false;
	for(register int i=0;i<taille;i++) if (vect[i]!=c.vect[i]) return false;
	return true;
}

bool DVector::operator==(const WVector &c)  const{ 
	if(taille!=c.taille) return false;
	if(taille==0 && c.taille==0) return true;
	if((c.vect==NULL && vect!=NULL) || (vect==NULL && c.vect!=NULL)) return false;
	for(register int i=0;i<taille;i++) if (vect[i]!=(double)c.vect[i]) return false;
	return true;
}

bool DVector::operator!=(const DVector &c)  const{ 
	if(taille!=c.taille) return true;
	if(taille==0 && c.taille==0) return false;
	if((c.vect==NULL && vect!=NULL) || (vect==NULL && c.vect!=NULL)) return true;
	for(register int i=0;i<taille;i++) if (vect[i]!=c.vect[i]) return true;
	return false;
}

bool DVector::operator!=(const WVector &c)  const{ 
	if(taille!=c.taille) return true;
	if(taille==0 && c.taille==0) return false;
	if((c.vect==NULL && vect!=NULL) || (vect==NULL && c.vect!=NULL)) return true;
	for(register int i=0;i<taille;i++) if (vect[i]!=(double)c.vect[i]) return true;
	return false;
}

double DVector::DotProd(const DVector &rhs){
	if (taille!=rhs.taille) return 0.0;
	register double sum=0.0;
	for(register int i=0;i<taille;i++) {
		sum+=vect[i]*rhs.vect[i];
	}
	return sum;
}

double DVector::DotProd(const WVector &rhs){
	if (taille!=rhs.taille) return 0.0;
	register double sum=0.0;
	for(register int i=0;i<taille;i++) {
		sum+=vect[i]*(double)rhs.vect[i];
	}
	return sum;
}

DCplx  DVector::DotProd(const ZVector &rhs){
	if (taille!=rhs.taille) return DCplx();
	register DCplx sum(0.0,0.0);
	for(register int i=0;i<taille;i++) {
		sum+=rhs.vect[i]*vect[i];
	}
	return sum;
}

void DVector::fixedconv(const DVector &rhs) {
  register double sum;
  int cache;

  /* Loop over output points */
  for (register int m = taille-1; m >= 0 ; m--) {
	 /* Convolution */
	 sum = 0.0;
	 for (register int j = 0; j < rhs.taille; ++j) {
		 cache=m-j;
		if (cache>=0 && cache<taille) sum += rhs.vect[j] * vect[cache];
	 }
	 vect[m] = sum;
  }
}

void DVector::conv(const DVector &rhs) {
  register double sum;
  register int cache;
  int imax=taille+rhs.taille-1;
  double *tmp=new double[imax];

  /* Loop over output points */
  for (register int i = 0; i<imax ; i++) {
	 /* Convolution */
	 sum = 0.0;
	 for (register int j = 0; j < rhs.taille; ++j) {
		cache=i-j;
		if (cache>=0 && i<taille) 
			sum += rhs.vect[j] * vect[cache];
	 }
	 tmp[i] = sum;
  }
  delete[] vect;
  taille=imax;
  vect=tmp;
}

DVector DVector::conv_nip(const DVector &rhs) const{
  register double sum;
  register int cache;
  int imax=taille+rhs.taille-1;
  DVector n(imax);

  /* Loop over output points */
  for (register int i = 0; i<imax ; i++) {
	 /* Convolution */
	 sum = 0.0;
	 for (register int j = 0; j < rhs.taille; ++j) {
		cache=i-j;
		if (cache>=0 && i<taille) 
			sum += rhs.vect[j] * vect[cache];
	 }
	 n.vect[i] = sum;
  }
  return n;
}

ZVector DVector::conv_nip(const ZVector &rhs) const{
  register DCplx sum;
  register int cache;
  int imax=taille+rhs.taille-1;
  ZVector n(imax);

  /* Loop over output points */
  for (register int i = 0; i<imax ; i++) {
	 /* Convolution */
	 sum = 0.0;
	 for (register int j = 0; j < rhs.taille; ++j) {
		cache=i-j;
		if (cache>=0 && i<taille) 
			sum += rhs.vect[j] * vect[cache];
	 }
	 n.vect[i] = sum;
  }
  return n;
}

DVector DVector::copy(int start, int stop) {
	stop = start > stop ? start : stop;
	start = start > stop ? stop : start;
	if (stop>taille) stop=taille; 
	if (start<0) start=0;
	if (stop-start) return DVector();
	DVector n(stop-start);
	for(register int i=start, j=0;i<stop;i++,j++)
		n.vect[j]=vect[i];
	return n;
}

DVector DVector::prepad(int count) {
	if (count<=0) return DVector();
	DVector n(count+taille);
	for(register int i=count, j=0; j<taille;i++,j++)
		n.vect[i]=vect[j];
	return n;
}

DVector DVector::postpad(int count) {
	if (count<=0) return DVector();
	DVector n(count+taille);
	for(register int j=0; j<taille;j++)
		n.vect[j]=vect[j];
	return n;
}

void DVector::insert(const DVector &rhs, int pos) {
	if(vect==NULL || rhs.vect==NULL) return;
	for(register int i=0, j=pos; i<pos && j<taille;i++, j++)
		vect[j]=rhs.vect[i];
}


ZVector::ZVector(const double *real, const double* imag, const int &size){
	if (size==0 || real==NULL || imag==NULL) { vect=NULL; taille=0; return; }
	taille=size;
	vect=new DCplx[taille];
	for(register int i=0;i<size;i++) {
		vect[i].re=real[i];
		vect[i].im=imag[i];
	}
}

ZVector::ZVector(const int &size){
	if (size==0) { vect=NULL; taille=0; return; }
	vect=new cmplx[size];
	taille=size;
}

ZVector::ZVector(const int &size, const DCplx &initial) {
	if (size==0) { vect=NULL; taille=0; return; }
	vect=new cmplx[size];
	taille=size;
	for(register int i=0;i<size;i++) vect[i]=initial;
}

ZVector::ZVector(const cmplx *v, const int size) {
	if (size==0) { vect=NULL; taille=0; return; }
	vect=new cmplx[size];
	taille=size;			
	for(register int i=0;i<size;i++) vect[i]=v[i];
}

ZVector::ZVector(const ZVector &c){
	if(c.taille==0) {
		taille=0; vect=NULL; 
		return;
	}
	taille=c.taille;
	vect=new cmplx[taille];
	for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
}

ZVector::ZVector(const DVector &c){
	if(c.taille==0) {
		taille=0; vect=NULL; 
		return;
	}
	taille=c.taille;
	vect=new cmplx[taille];
	for(register int i=0;i<taille;i++) vect[i].re=c.vect[i];
}

ZVector::ZVector(const WVector &c){
	if(c.taille==0) {
		taille=0; vect=NULL; 
		return;
	}
	taille=c.taille;
	vect=new cmplx[taille];
	for(register int i=0;i<taille;i++) vect[i].re=(double)c.vect[i];
}


ZVector::ZVector() {
	taille=0; vect=NULL;
}

ZVector::~ZVector() {
	if (vect!=NULL) delete[] vect;
	vect=NULL;taille=0;
}

//cmplx        conjug(cmplx a);
ZVector ZVector::conj() const {
	ZVector t(taille);
	for(register int i=0;i<taille;i++) {
		t.vect[i]=vect[i].conj();
	}
	return t;
}

DVector ZVector::real() const  { 
	if(taille==0 || vect==NULL) return NULL;
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i].re;
	}
	return n;
}

DVector ZVector::imag()  const { 
	if(taille==0 || vect==NULL) return NULL;
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i].im;
	}
	return n;
}

//double       modul2(cmplx a);
DVector ZVector::mag() const  { 
	if(taille==0 || vect==NULL) return NULL;
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i].mag();
	}
	return n;
}

DVector ZVector::mod2() const  { 
	if(taille==0 || vect==NULL) return NULL;
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i].mod2();
	}
			
	return n;
}

DVector ZVector::phase() const  {
	if(taille==0 || vect==NULL) return NULL;
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i].phase();
	}
	return n;
}

void ZVector::fixedconv(const ZVector &rhs) {
  register DCplx sum;
  int cache;

  /* Loop over output points */
  for (register int m = taille-1; m >= 0 ; m--) {
	 /* Convolution */
	 sum = 0.0;
	 for (register int j = 0; j < rhs.taille; ++j) {
		 cache=m-j;
		if (cache>=0 && cache<taille) sum += rhs.vect[j] * vect[cache];
	 }
	 vect[m] = sum;
  }
}

void ZVector::conv(const ZVector &rhs) {
  register DCplx sum;
  register int cache;
  int imax=taille+rhs.taille-1;
  DCplx *tmp=new DCplx[imax];

  /* Loop over output points */
  for (register int i = 0; i<imax ; i++) {
	 /* Convolution */
	 sum = 0.0;
	 for (register int j = 0; j < rhs.taille; ++j) {
		cache=i-j;
		if (cache>=0 && i<taille) 
			sum += rhs.vect[j] * vect[cache];
	 }
	 tmp[i] = sum;
  }
  delete[] vect;
  taille=imax;
  vect=tmp;
}

ZVector ZVector::conv_nip(const ZVector &rhs) const{
  register DCplx sum;
  register int cache;
  int imax=taille+rhs.taille-1;
  ZVector n(imax);

  /* Loop over output points */
  for (register int i = 0; i<imax ; i++) {
	 /* Convolution */
	 sum = 0.0;
	 for (register int j = 0; j < rhs.taille; ++j) {
		cache=i-j;
		if (cache>=0 && i<taille) 
			sum += rhs.vect[j] * vect[cache];
	 }
	 n.vect[i] = sum;
  }
  return n;
}

void ZVector::fixedconv(const DVector &rhs) {
  register DCplx sum;
  int cache;

  /* Loop over output points */
  for (register int m = taille-1; m >= 0 ; m--) {
	 /* Convolution */
	 sum = 0.0;
	 for (register int j = 0; j < rhs.taille; ++j) {
		 cache=m-j;
		if (cache>=0 && cache<taille) sum += rhs.vect[j] * vect[cache];
	 }
	 vect[m] = sum;
  }
}

void ZVector::conv(const DVector &rhs) {
  register DCplx sum;
  register int cache;
  int imax=taille+rhs.taille-1;
  DCplx *tmp=new DCplx[imax];

  /* Loop over output points */
  for (register int i = 0; i<imax ; i++) {
	 /* Convolution */
	 sum = 0.0;
	 for (register int j = 0; j < rhs.taille; ++j) {
		cache=i-j;
		if (cache>=0 && i<taille) 
			sum += rhs.vect[j] * vect[cache];
	 }
	 tmp[i] = sum;
  }
  delete[] vect;
  taille=imax;
  vect=tmp;
}

ZVector ZVector::conv_nip(const DVector &rhs) const{
  register DCplx sum;
  register int cache;
  int imax=taille+rhs.taille-1;
  ZVector n(imax);

  /* Loop over output points */
  for (register int i = 0; i<imax ; i++) {
	 /* Convolution */
	 sum = 0.0;
	 for (register int j = 0; j < rhs.taille; ++j) {
		cache=i-j;
		if (cache>=0 && i<taille) 
			sum += rhs.vect[j] * vect[cache];
	 }
	 n.vect[i] = sum;
  }
  return n;
}

ZVector ZVector::copy(int start, int stop) {
	stop = start > stop ? start : stop;
	start = start > stop ? stop : start;
	if (stop>taille) stop=taille; 
	if (start<0) start=0;
	if (!(stop-start)) return ZVector();
	ZVector n(stop-start);
	for(register int i=start, j=0;i<stop;i++,j++)
		n.vect[j]=vect[i];
	return n;
}

ZVector ZVector::prepad(int count) {
	if (count<=0) return ZVector();
	ZVector n(count+taille);
	for(register int i=count, j=0; j<taille;i++,j++)
		n.vect[i]=vect[j];
	return n;
}

ZVector ZVector::postpad(int count) {
	if (count<=0) return ZVector();
	ZVector n(count+taille);
	for(register int j=0; j<taille;j++)
		n.vect[j]=vect[j];
	return n;
}

void ZVector::insert(const ZVector &rhs, int pos) {
	if(vect==NULL || rhs.vect==NULL) return;
	for(register int i=0, j=pos; i<pos && j<taille;i++, j++)
		vect[j]=rhs.vect[i];
}

void ZVector::insert(const DVector &rhs, int pos) {
	if(vect==NULL || rhs.vect==NULL) return;
	for(register int i=0, j=pos; i<pos && j<taille;i++, j++)
		vect[j]=rhs.vect[i];
}

//cmplx        cmplx_add(cmplx a,cmplx b);
ZVector ZVector::operator+(const ZVector &c)  const { 
	if (c.taille != taille) return ZVector();
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c.vect[i];
	}
	return n;
}

ZVector ZVector::operator+(const DCplx &c)  const { 
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c;
	}
	return n;
}

ZVector ZVector::operator+(const double &c)  const { 
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c;
	}
	return n;
}

ZVector ZVector::operator+(const DVector &c)  const { 
	if (c.taille != taille) return ZVector();
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c.vect[i];
	}
	return n;
}

ZVector ZVector::operator+(const WVector &c)  const { 
	if (c.taille != taille) return ZVector();
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+(double)c.vect[i];
	}
	return n;
}

// cmplx        cmplx_sub(cmplx a,cmplx b);
ZVector ZVector::operator-(const ZVector &c)  const { 
	if (c.taille != taille) return ZVector();
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c.vect[i];
	}
	return n;
}

ZVector ZVector::operator-(const DCplx &c)  const { 
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c;
	}
	return n;
}

ZVector ZVector::operator-(const double &c)  const { 
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c;
	}
	return n;
}

ZVector ZVector::operator-(const DVector &c) const  { 
	if (c.taille != taille) return ZVector();
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c.vect[i];
	}
	return n;
}

ZVector ZVector::operator-(const WVector &c)  const { 
	if (c.taille != taille) return ZVector();
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-(double)c.vect[i];
	}
	return n;
}


//cmplx        cmplx_mul(cmplx a,cmplx b);
ZVector ZVector::operator*(const ZVector &c)  const { 
	if (c.taille != taille) return ZVector();
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c.vect[i];
	}
	return n;
}

ZVector ZVector::operator*(const DCplx &c)  const { 
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c;
	}
	return n;
}

ZVector ZVector::operator*(const double &c)  const { 
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c;
	}
	return n;
}

ZVector ZVector::operator*(const DVector &c)  const { 
	if (c.taille != taille) return ZVector();
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c.vect[i];
	}
	return n;
}

ZVector ZVector::operator*(const WVector &c)  const { 
	if (c.taille != taille) return ZVector();
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*(double)c.vect[i];
	}
	return n;
}

ZVector ZVector::operator/(const double &c)  const { 
	if (c==0.0) return ZVector();
	double d=1/c;
	return ZVector((*this)*d);
}

ZVector ZVector::operator/(const DCplx &c)  const {
	double t=c.mod2();
	if (t==0) return ZVector();
	DCplx ct=c.conj()/t;
	return ZVector((*this)*ct);
}

ZVector &ZVector::operator=(const ZVector &c) { 
	if(taille!=0 && vect!=NULL) delete[] vect;
	taille=c.taille;
	vect = new cmplx[taille];
	for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
	return *this;
}

ZVector &ZVector::operator=(const DVector &c){
	if(taille!=0 && vect!=NULL) delete[] vect;
	taille=c.taille;
	vect = new cmplx[taille];
	for(register int i=0;i<taille;i++) {
		vect[i].re=c.vect[i];
	}
	return *this;
}

ZVector &ZVector::operator=(const WVector &c){
	if(taille!=0 && vect!=NULL) delete[] vect;
	taille=c.taille;
	vect = new cmplx[taille];
	for(register int i=0;i<taille;i++) {
		vect[i].re=(double)c.vect[i];
	}
	return *this;
}

bool ZVector::operator==(const DVector &c){
	if(taille!=c.taille) return false;
	register int i=0;
	while(i<taille && vect[i]==c.vect[i]) i++;
	return (i==taille);
}

bool ZVector::operator==(const WVector &c){
	if(taille!=c.taille) return false;
	register int i=0;
	while(i<taille && vect[i]==c.vect[i]) i++;
	return (i==taille);
}

bool ZVector::operator==(const ZVector &c) {
	if(taille!=c.taille) return false;
	register int i=0;
	while(i<taille && vect[i]==c.vect[i]) i++;
	return (i==taille);
}

bool ZVector::operator!=(const DVector &c) {
	if(taille!=c.taille) return true;
	register int i=0;
	while(vect[i]!=c.vect[i] && i<taille) i++;
	return !(i==taille);
}

bool ZVector::operator!=(const WVector &c) {
	if(taille!=c.taille) return true;
	register int i=0;
	while(vect[i]!=c.vect[i] && i<taille) i++;
	return !(i==taille);
}

bool ZVector::operator!=(const ZVector &c) {
	if(taille!=c.taille) return true;
	register int i=0;
	while(vect[i]!=c.vect[i] && i<taille) i++;
	return !(i==taille);
}

DCplx ZVector::DotProd(const DVector &rhs){
	if (taille!=rhs.taille) return DCplx();
	DCplx sum(0.0,0.0);
	for(register int i=0;i<taille;i++) {
		sum+=vect[i]*rhs.vect[i];
	}
	return sum;
}

DCplx ZVector::DotProd(const WVector &rhs){
	if (taille!=rhs.taille) return DCplx();
	DCplx sum(0.0,0.0);
	for(register int i=0;i<taille;i++) {
		sum+=vect[i]*(double)rhs.vect[i];
	}
	return sum;
}

DCplx ZVector::DotProd(const ZVector &rhs){
	if (taille!=rhs.taille) return DCplx();
	DCplx sum(0.0,0.0);
	for(register int i=0;i<taille;i++) {
		sum+=rhs.vect[i]*vect[i];
	}
	return sum;
}

DMatrix::DMatrix(int line, int col, double initial){
	if (col==0 || line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=col;this->line=line;
	mat=new double*[line];
	for(register int i =0;i<line;i++) {
		mat[i]=new double[col];
		for(register int j=0;j<col;j++) mat[i][j]=initial;
	}
}

DMatrix::DMatrix(const DMatrix& m) {
	if (m.col==0 || m.line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=m.col;this->line=m.line;
	mat=new double*[line];
	size_t col_size=sizeof(double)*col;
	for(register int i =0;i<line;i++) {
		mat[i]=new double[col];
		memcpy(mat[i], m.mat[i], col_size);
	}
}

DMatrix::DMatrix() { col=line=0;mat=NULL; }
//void Free();
DMatrix::~DMatrix() {
	for(register int i=0;i<line;i++) 
		delete[] mat[i];
	delete[] mat;
	line=col=0;
	mat=NULL;
}

void DMatrix::transpose() {
	if (mat==NULL ||col==0 ||line==0) return;
	double **tmp = new double*[col];
	for(register int i=0; i<col;i++) {
		tmp[i]=new double[line];
		for(register int j=0;j<line;j++)
			tmp[i][j]=mat[j][i];
	}
	for(i=0;i<line;i++) delete[] mat[i];
	delete[] mat;

	int x=col;col=line;line=x;
	mat=tmp;
}

DMatrix DMatrix::transpose_nip(){
	if (mat==NULL ||col==0 ||line==0) return DMatrix();
	DMatrix t; 
	t.col=line;t.line=col;

	t.mat = new double*[col];
	for(register int i=0; i<col;i++) {
		t.mat[i]=new double[line];
		for(register int j=0;j<line;j++)
			t.mat[i][j]=mat[j][i];
	}
	return t;
}

ZMatrix DMatrix::ztranspose_nip(){
	if (mat==NULL ||col==0 ||line==0) return ZMatrix();
	ZMatrix t; t.col=line;t.line=col;

	t.mat = new DCplx*[col];
	for(register int i=0; i<col;i++) {
		t.mat[i]=new DCplx[line];
		for(register int j=0;j<line;j++)
			t.mat[i][j].re=mat[j][i];
	}
	return t;
}


DMatrix DMatrix::operator+(const DMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) return DMatrix();
	if (rhs.col!=col || rhs.line!=line) return DMatrix();
	DMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]+rhs.mat[i][j];
	}
	return n;
}

DMatrix DMatrix::operator-(const DMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) return DMatrix();
	if (rhs.col!=col || rhs.line!=line) return DMatrix();
	DMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]-rhs.mat[i][j];
	}
	return n;
}

DMatrix DMatrix::operator*(const DMatrix &rhs){
	if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return DMatrix();
	if (rhs.line!=col) return DMatrix();
	//DMatrix n(rhs.col, line);
	DMatrix n(line, rhs.col);

	for(register int j=0;j<n.line;j++) { //this.line
		for(register int i=0;i<n.col;i++) { // rhs.col
			register double sum=0.0;
			for(register k=0;k<col;k++)  // this.col
				sum+=mat[j][k]*rhs.mat[k][i];
			n.mat[j][i]=sum;
		}
	}

	return n;
}

ZMatrix DMatrix::operator+(const ZMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) return ZMatrix();
	if (rhs.col!=col || rhs.line!=line) return ZMatrix();
	ZMatrix n(line,  col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]+rhs.mat[i][j];
	}
	return n;
}

ZMatrix DMatrix::operator-(const ZMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) return ZMatrix();
	if (rhs.col!=col || rhs.line!=line) return ZMatrix();
	ZMatrix n(line,  col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]-rhs.mat[i][j];
	}
	return n;
}

ZMatrix DMatrix::operator*(const ZMatrix &rhs){
	if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return ZMatrix();
	if (rhs.line!=col) return ZMatrix();
	ZMatrix n(rhs.col, line);

	for(register int j=0;j<n.line;j++) {
		for(register int i=0;i<n.col;i++) {
			register DCplx sum(0.0,0.0);
			for(register k=0;k<col;k++) 
				sum+=mat[j][k]*rhs.mat[k][i];
			n.mat[j][i]=sum;
		}
	}
	return n;
}

DMatrix DMatrix::operator*(const double &rhs){
	if (col==0 || line==0 || mat==NULL) return DMatrix();
	DMatrix n(line,  col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]*rhs;
	}
	return n;
}

ZMatrix DMatrix::operator*(const DCplx &rhs){
	if (col==0 || line==0 || mat==NULL) return ZMatrix();
	ZMatrix n(line,  col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]*rhs;
	}
	return n;
}

DMatrix DMatrix::operator/(const double &rhs){
	if(rhs==0.0) return DMatrix();
	double irhs=1/rhs;
	return (*this)*irhs;
}

ZMatrix DMatrix::operator/(const DCplx &rhs){
	if(rhs==0.0) return ZMatrix();
	DCplx irhs=1/rhs;
	return (*this)*irhs;
}

DMatrix &DMatrix::operator=(const DMatrix &rhs){
	if(mat!=NULL) {
		for(register int i=0;i<line;i++) delete[] mat[i];
		delete[] mat;
	}

	line=rhs.line;col=rhs.col;

	if(line==0|| col==0 || rhs.mat==NULL) {
		this->mat=NULL; return *this;
	}

	mat=new double*[line];
	for(register int i = 0;i<line;i++) {
		mat[i]=new double[col];
		for(register int j=0;j<col;j++)
			mat[i][j]=rhs.mat[i][j];
	}
	return *this;

}

bool DMatrix::operator==(const DMatrix &rhs){
	if(rhs.col!=col || rhs.line!=line) return false;
	for(int i=0;i<line;i++){
		for(int j=0;j<col;j++)
			if (rhs.mat[i][j]!=mat[i][j]) return false;
	}
	return true;
}

bool DMatrix::operator!=(const DMatrix &rhs){
	if(rhs.col!=col || rhs.line!=line) return true;
	for(int i=0;i<line;i++){
		for(int j=0;j<col;j++)
			if (rhs.mat[i][j]!=mat[i][j]) return true;
	}
	return false;
}


ZMatrix::ZMatrix (int line, int col){
	if (col==0 || line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=col;this->line=line;
	mat=new DCplx*[line];
	DCplx initial;
	for(register int i =0;i<line;i++) {
		mat[i]=new DCplx[col];
		for(register int j=0;j<col;j++) mat[i][j]=initial;
	}
}

ZMatrix::ZMatrix (int line, int col, const DCplx &initial)  {
	if (col==0 || line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=col;this->line=line;
	mat=new DCplx*[line];
	for(register int i =0;i<line;i++) {
		mat[i]=new DCplx[col];
		for(register int j=0;j<col;j++) mat[i][j]=initial;
	}
}

ZMatrix::ZMatrix(const ZMatrix& m){
	if (m.col==0 || m.line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=m.col;this->line=m.line;
	mat=new DCplx*[line];
	
	for(register int i =0;i<line;i++) {
		mat[i]=new DCplx[col];
		for(register int j=0;j<col;j++) mat[i][j]=m.mat[i][j];
	}
}

ZMatrix::ZMatrix(const DMatrix& m){
	if (m.col==0 || m.line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=m.col;this->line=m.line;
	mat=new DCplx*[line];
	
	for(register int i =0;i<line;i++) {
		mat[i]=new DCplx[col];
		for(register int j=0;j<col;j++) mat[i][j]=m.mat[i][j];
	}
}

ZMatrix::ZMatrix() { col=line=0;mat=NULL;}

ZMatrix::~ZMatrix() {
	if(col==0 || mat==0) return;
	if(line==0) {delete[] mat; return;}
	for(register int i =0; i<line;i++) delete[] mat[i];
	delete[] mat;
	line=col=0;
	mat=NULL;
}
//void Free();
void ZMatrix::transpose(){
	if (mat==NULL ||col==0 ||line==0) return;
	DCplx **tmp = new DCplx *[col];
	for(register int i=0; i<col;i++) {
		tmp[i]=new DCplx[line];
		for(register int j=0;j<line;j++)
			tmp[i][j]=mat[j][i];
	}
	for(i=0;i<line;i++) delete[] mat[i];
	delete[] mat;

	int x=col;col=line;line=x;
	mat=tmp;
}

ZMatrix  ZMatrix::transpose_nip() const {
	if (mat==NULL ||col==0 ||line==0) return ZMatrix();
	ZMatrix t; t.col=line;t.line=col;

	t.mat = new DCplx*[col];
	for(register int i=0; i<col;i++) {
		t.mat[i]=new DCplx[line];
		for(register int j=0;j<line;j++)
			t.mat[i][j]=mat[j][i];
	}
	return t;
}

ZMatrix  ZMatrix::operator+(const ZMatrix &rhs) const {
	if (col==0 || line==0 || mat==NULL) return ZMatrix();
	if (rhs.col!=col || rhs.line!=line) return ZMatrix();
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]+rhs.mat[i][j];
	}
	return n;
}

ZMatrix  ZMatrix::operator+(const DMatrix &rhs) const {
	if (col==0 || line==0 || mat==NULL) return ZMatrix();
	if (rhs.col!=col || rhs.line!=line) return ZMatrix();
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]+rhs.mat[i][j];
	}
	return n;
}

ZMatrix  ZMatrix::operator-(const ZMatrix &rhs) const {
	if (col==0 || line==0 || mat==NULL) return ZMatrix();
	if (rhs.col!=col || rhs.line!=line) return ZMatrix();
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]-rhs.mat[i][j];
	}
	return n;
}

ZMatrix  ZMatrix::operator-(const DMatrix &rhs) const {
	if (col==0 || line==0 || mat==NULL) return ZMatrix();
	if (rhs.col!=col || rhs.line!=line) return ZMatrix();
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]-rhs.mat[i][j];
	}
	return n;
}

ZMatrix ZMatrix::operator*(const ZMatrix &rhs) const {
	if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return ZMatrix();
	if (rhs.line!=col) return ZMatrix();

	ZMatrix n(line , rhs.col);

	for(register int j=0;j<n.line;j++) {
		for(register int i=0;i<n.col;i++) {
			register DCplx sum(0.0,0.0);
			for(register k=0;k<col;k++) 
				sum+=mat[j][k]*rhs.mat[k][i];
			n.mat[j][i]=sum;
		}
	}
	return n;
}

ZMatrix ZMatrix::operator*(const DMatrix &rhs) const {
	if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return ZMatrix();
	if (rhs.line!=col) return ZMatrix();
	ZMatrix n(rhs.col, line);

	for(register int j=0;j<n.line;j++) {
		for(register int i=0;i<n.col;i++) {
			register DCplx sum(0.0,0.0);
			for(register k=0;k<col;k++) 
				sum+=mat[j][k]*rhs.mat[k][i];
			n.mat[j][i]=sum;
		}
	}
	return n;
}

ZMatrix ZMatrix::operator*(const double &rhs) const {
	if (col==0 || line==0 || mat==NULL) return ZMatrix();
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]*rhs;
	}
	return n;
}

ZMatrix ZMatrix::operator*(const DCplx &rhs) const {
	if (col==0 || line==0 || mat==NULL) return ZMatrix();
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]*rhs;
	}
	return n;
}


ZMatrix ZMatrix::operator/(const double &rhs) const {
	if(rhs==0.0) return ZMatrix();
	double irhs=1/rhs;
	return (*this)*irhs;
}

ZMatrix ZMatrix::operator/(const DCplx &rhs) const {
	if(rhs==0.0) return ZMatrix();
	DCplx irhs=1/rhs;
	return (*this)*irhs;
}

ZMatrix &ZMatrix::operator=(const ZMatrix &rhs) {
	if(mat!=NULL) {
		for(register int i=0;i<col;i++) delete[] mat[i];
		delete[] mat;
	}

	line=rhs.line;col=rhs.col;

	if(line==0|| col==0 || rhs.mat==NULL) {
		this->mat=NULL; line=col=0; return *this;
	}

	mat=new DCplx*[line];
	for(register int i = 0;i<line;i++) {
		mat[i]=new DCplx[col];
		for(register int j=0;j<col;j++)
			mat[i][j]=rhs.mat[i][j];
	}
	return *this;
}

ZMatrix &ZMatrix::operator=(const DMatrix &rhs){
	if(mat!=NULL) {
		for(register int i=0;i<col;i++) delete[] mat[i];
		delete[] mat;
	}

	line=rhs.line;col=rhs.col;

	if(line==0|| col==0 || rhs.mat==NULL) {
		this->mat=NULL; line=col=0; return *this;
	}

	mat=new DCplx*[line];
	for(register int i = 0;i<line;i++) {
		mat[i]=new DCplx[col];
		for(register int j=0;j<col;j++)
			mat[i][j]=rhs.mat[i][j];
	}
	return *this;
}

bool ZMatrix::operator==(const ZMatrix &rhs) const {
	if(rhs.col!=col || rhs.line!=line) return false;
	for(int i=0;i<line;i++){
		for(int j=0;j<col;j++)
			if (rhs.mat[i][j]!=mat[i][j]) return false;
	}
	return true;
}

bool ZMatrix::operator!=(const ZMatrix &rhs) const {
	if(rhs.col!=col || rhs.line!=line) return true;
	for(int i=0;i<line;i++){
		for(int j=0;j<col;j++)
			if (rhs.mat[i][j]!=mat[i][j]) return true;
	}
	return false;
}

DMatrix ZMatrix::real() const {
	if(col==0 || line==0 || mat==NULL) return DMatrix();
	DMatrix n(line, col);
	n.mat=new double*[line];
	for(register int i = 0;i<line;i++) {
		n.mat[i]=new double[col];
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j].re;
	}
	return n;
}

DMatrix ZMatrix::imag() const {
	if(col==0 || line==0 || mat==NULL) return DMatrix();
	DMatrix n(line, col);
	n.mat=new double*[line];
	for(register int i = 0;i<line;i++) {
		n.mat[i]=new double[col];
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j].im;
	}
	return n;
}

DMatrix ZMatrix::mag() const {
	if(col==0 || line==0 || mat==NULL) return DMatrix();
	DMatrix n(line, col);
	n.mat=new double*[line];
	for(register int i = 0;i<line;i++) {
		n.mat[i]=new double[col];
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j].mag();
	}
	return n;
}

DMatrix ZMatrix::phase() const {
	if(col==0 || line==0 || mat==NULL) return DMatrix();
	DMatrix n(line, col);
	n.mat=new double*[line];
	for(register int i = 0;i<line;i++) {
		n.mat[i]=new double[col];
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j].phase();
	}
	return n;
}

ZMatrix ZMatrix::conj() const {
	if(col==0 || line==0 || mat==NULL) return ZMatrix();
	ZMatrix n(line, col);
	n.mat=new DCplx*[line];
	for(register int i = 0;i<line;i++) {
		n.mat[i]=new DCplx[col];
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j].conj();
	}
	return n;
}

inline ZMatrix operator*(const DCplx &lhs, const ZMatrix &rhs){
	return rhs*lhs;
}

inline ZMatrix operator*(const double &lhs, const ZMatrix &rhs){
	return rhs*lhs;
}

DMatrix operator*(const double &lhs, const DMatrix &rhs){
	if (rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return DMatrix();
	DMatrix n(rhs.line, rhs.col);
	for(register int i=0;i<rhs.line;i++) {
		for(register int j=0;j<rhs.col;j++)
			n.mat[i][j]=rhs.mat[i][j]*lhs;
	}
	return n;
}

ZMatrix operator*(const DCplx &lhs, const DMatrix &rhs){
	if (rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return ZMatrix();
	ZMatrix n(rhs.line, rhs.col);
	for(register int i=0;i<rhs.line;i++) {
		for(register int j=0;j<rhs.col;j++)
			n.mat[i][j]=rhs.mat[i][j]*lhs;
	}
	return n;
}

WMatrix::WMatrix(int line, int col, int initial){
	if (col==0 || line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=col;this->line=line;
	mat=new int*[line];
	for(register int i =0;i<line;i++) {
		mat[i]=new int[col];
		for(register int j=0;j<col;j++) mat[i][j]=initial;
	}
}

WMatrix::WMatrix(const WMatrix& m) {
	if (m.col==0 || m.line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=m.col;this->line=m.line;
	mat=new int*[line];
	size_t col_size=sizeof(double)*col;
	for(register int i =0;i<line;i++) {
		mat[i]=new int[col];
		memcpy(mat[i], m.mat[i], col_size);
	}
}

WMatrix::WMatrix() { col=line=0;mat=NULL; }

WMatrix::~WMatrix() {
	for(register int i=0;i<line;i++) 
		delete[] mat[i];
	delete[] mat;
	line=col=0;
	mat=NULL;
}

void WMatrix::transpose() {
	if (mat==NULL ||col==0 ||line==0) return;
	int **tmp = new int*[col];
	for(register int i=0; i<col;i++) {
		tmp[i]=new int[line];
		for(register int j=0;j<line;j++)
			tmp[i][j]=mat[j][i];
	}
	for(i=0;i<line;i++) delete[] mat[i];
	delete[] mat;

	int x=col;col=line;line=x;
	mat=tmp;
}

WMatrix WMatrix::transpose_nip(){
	if (mat==NULL ||col==0 ||line==0) return WMatrix();
	WMatrix t; t.col=line;t.line=col;

	t.mat = new int*[col];
	for(register int i=0; i<col;i++) {
		t.mat[i]=new int[line];
		for(register int j=0;j<line;j++)
			t.mat[i][j]=mat[j][i];
	}
	return t;
}

ZMatrix WMatrix::ztranspose_nip(){
	if (mat==NULL ||col==0 ||line==0) return ZMatrix();
	ZMatrix t; t.col=line;t.line=col;

	t.mat = new DCplx*[col];
	for(register int i=0; i<col;i++) {
		t.mat[i]=new DCplx[line];
		for(register int j=0;j<line;j++)
			t.mat[i][j].re=(double)mat[j][i];
	}
	return t;
}


WMatrix WMatrix::operator+(const WMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) return WMatrix();
	if (rhs.col!=col || rhs.line!=line) return WMatrix();
	WMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]+rhs.mat[i][j];
	}
	return n;
}

WMatrix WMatrix::operator-(const WMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) return WMatrix();
	if (rhs.col!=col || rhs.line!=line) return WMatrix();
	WMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]-rhs.mat[i][j];
	}
	return n;
}

WMatrix WMatrix::operator*(const WMatrix &rhs){
	if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return WMatrix();
	if (rhs.line!=col) return WMatrix();
	WMatrix n(line, rhs.col);

	for(register int j=0;j<n.line;j++) {
		for(register int i=0;i<n.col;i++) {
			register int sum=0;
			for(register k=0;k<col;k++) 
				sum+=mat[j][k]*rhs.mat[k][i];
			n.mat[j][i]=sum;
		}
	}
	return n;
}

DMatrix WMatrix::operator+(const DMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) return DMatrix();
	if (rhs.col!=col || rhs.line!=line) return DMatrix();
	DMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=(double)mat[i][j]+rhs.mat[i][j];
	}
	return n;
}

DMatrix WMatrix::operator-(const DMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) return DMatrix();
	if (rhs.col!=col || rhs.line!=line) return DMatrix();
	DMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=(double)mat[i][j]-rhs.mat[i][j];
	}
	return n;
}

DMatrix WMatrix::operator*(const DMatrix &rhs){
	if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return DMatrix();
	if (rhs.line!=col) return DMatrix();
	DMatrix n(rhs.col, line);

	for(register int j=0;j<n.line;j++) {
		for(register int i=0;i<n.col;i++) {
			register double sum=0.0;
			for(register k=0;k<col;k++) 
				sum+=(double)mat[j][k]*rhs.mat[k][i];
			n.mat[j][i]=sum;
		}
	}
	return n;
}

ZMatrix WMatrix::operator+(const ZMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) return ZMatrix();
	if (rhs.col!=col || rhs.line!=line) return ZMatrix();
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j].re=mat[i][j]+rhs.mat[i][j].re;
	}
	return n;
}

ZMatrix WMatrix::operator-(const ZMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) return ZMatrix();
	if (rhs.col!=col || rhs.line!=line) return ZMatrix();
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j].re=mat[i][j]-rhs.mat[i][j].re;
	}
	return n;
}

ZMatrix WMatrix::operator*(const ZMatrix &rhs){
	if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) return ZMatrix();
	if (rhs.line!=col) return ZMatrix();
	ZMatrix n(rhs.col, line);

	for(register int j=0;j<n.line;j++) {
		for(register int i=0;i<n.col;i++) {
			register DCplx sum(0.0,0.0);
			for(register k=0;k<col;k++) 
				sum+=mat[j][k]*rhs.mat[k][i];
			n.mat[j][i]=sum;
		}
	}
	return n;
}

WMatrix WMatrix::operator*(const int &rhs){
	if (col==0 || line==0 || mat==NULL) return WMatrix();
	WMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]*rhs;
	}
	return n;
}

DMatrix WMatrix::operator*(const double &rhs){
	if (col==0 || line==0 || mat==NULL) return DMatrix();
	DMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=double(mat[i][j])*rhs;
	}
	return n;
}

ZMatrix WMatrix::operator*(const DCplx &rhs){
	if (col==0 || line==0 || mat==NULL) return ZMatrix();
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]*rhs;
	}
	return n;
}

DMatrix WMatrix::operator/(const double &rhs){
	if(rhs==0.0) return DMatrix();
	double irhs=1/rhs;
	return (*this)*irhs;
}

ZMatrix WMatrix::operator/(const DCplx &rhs){
	if(rhs==0.0) return ZMatrix();
	DCplx irhs=1/rhs;
	return (*this)*irhs;
}

WMatrix &WMatrix::operator=(const WMatrix &rhs){
	if(mat!=NULL) {
		for(register int i=0;i<line;i++) delete[] mat[i];
		delete[] mat;
	}

	line=rhs.line;col=rhs.col;

	if(line==0|| col==0 || rhs.mat==NULL) {
		this->mat=NULL; line=col=0; return *this;
	}

	mat=new int*[line];
	for(register int i = 0;i<line;i++) {
		mat[i]=new int[col];
		for(register int j=0;j<col;j++)
			mat[i][j]=rhs.mat[i][j];
	}
	return *this;

}

bool WMatrix::operator==(const WMatrix &rhs){
	if(rhs.col!=col || rhs.line!=line) return false;
	for(int i=0;i<line;i++){
		for(int j=0;j<col;j++)
			if (rhs.mat[i][j]!=mat[i][j]) return false;
	}
	return true;
}

bool WMatrix::operator!=(const WMatrix &rhs){
	if(rhs.col!=col || rhs.line!=line) return true;
	for(int i=0;i<line;i++){
		for(int j=0;j<col;j++)
			if (rhs.mat[i][j]!=mat[i][j]) return true;
	}
	return false;
}

bool WMatrix::operator==(const DMatrix &rhs){
	if(rhs.col!=col || rhs.line!=line) return false;
	for(int i=0;i<line;i++){
		for(int j=0;j<col;j++)
			if (rhs.mat[i][j]!=mat[i][j]) return false;
	}
	return true;
}

bool WMatrix::operator!=(const DMatrix &rhs){
	if(rhs.col!=col || rhs.line!=line) return true;
	for(int i=0;i<line;i++){
		for(int j=0;j<col;j++)
			if (rhs.mat[i][j]!=mat[i][j]) return true;
	}
	return false;
}

ostream& operator<< ( ostream& os, DCplx& dt) {
   os << '(' << dt.re << ',' << dt.im << ')';
   return os;
}

ostream& operator<< ( ostream& os, DVector& dt ) {
	if (dt.vect==NULL || dt.taille==0) {
		os << "DVector is empty !";
		return os;
	}
	os << endl << "[";
	for(register int i=0;i<dt.taille;i++) {
		os << dt.vect[i]<< ' ';
	}
	os << "]" << endl;
	return os;
}

ostream& operator<< ( ostream& os, WVector& dt ) {
	if (dt.vect==NULL || dt.taille==0) {
		os << "WVector is empty !";
		return os;
	}
	os << endl << "[";

	for(register int i=0;i<dt.taille;i++) {
		os << dt.vect[i]<< ' ';
	}
	os << "]" << endl;

	return os;
}

ostream& operator<< ( ostream& os, ZVector& dt ) {
	if (dt.vect==NULL || dt.taille==0) {
		os << "ZVector is empty !";
		return os;
	}
	os << endl << "[";
	for(register int i=0;i<dt.taille;i++) {
		os << dt.vect[i]<< ' ';
	}
	os << "]" << endl;
	return os;
}

ostream& operator<< ( ostream& os, DMatrix& dt ) {
	if (dt.mat==NULL || dt.col==0 || dt.line==0) {
		os << "DMatrix is empty !";
		return os;
	}
	os << endl << "[";
	for(register int i=0;i<dt.line;i++) {
		os << "[";
		for(register int j=0;j<dt.col;j++) {
			os << dt.mat[i][j]<< ' ';
		}
		os << ']';
	}
	os << "]" << endl;
	return os;
}

ostream& operator<< ( ostream& os, WMatrix& dt ) {
	if (dt.mat==NULL || dt.col==0 || dt.line==0) {
		os << "WMatrix is empty !";
		return os;
	}
	os << endl << "[";
	for(register int i=0;i<dt.line;i++) {
		os << "[";
		for(register int j=0;j<dt.col;j++) {
			os << dt.mat[i][j]<< ' ';
		}
		os << ']';
	}
	os << "]" << endl;
	return os;
}

ostream& operator<< ( ostream& os, ZMatrix& dt ) {
	if (dt.mat==NULL || dt.col==0 || dt.line==0) {
		os << "ZMatrix is empty !";
		return os;
	}
	os << endl << "[";
	for(register int i=0;i<dt.line;i++) {
		os << "[";
		for(register int j=0;j<dt.col;j++) {
			os << dt.mat[i][j]<< ' ';
		}
		os << ']';
	}
	os << "]" << endl;
	return os;
}



/********************************************************/
/*cmplx init_cmplx(double a, double b)
{
	cmplx c;
	c.re=a;
	c.im=b;
	return c;
}
 
cmplx conjug(cmplx a)
{
	a.im=-a.im;
	return a;
}
 
double modul2(cmplx a)
{
	return a.re*a.re+a.im*a.im;
}
 
cmplx cmplx_add(cmplx a,cmplx b)
{
	cmplx c;
	c.re=a.re+b.re;
	c.im=a.im+b.im;
	return c;
}
 
cmplx cmplx_sub(cmplx a,cmplx b)
{
	cmplx c;
	c.re=a.re-b.re;
	c.im=a.im-b.im;
	return c;
}
 
cmplx cmplx_mul(cmplx a,cmplx b)
{
	cmplx c;
	double temp;
	temp=a.re*b.re-a.im*b.im;
	c.im=a.re*b.im+a.im*b.re;
	c.re=temp;
	return c;
}
 
double argum(cmplx x)
{
	double phase;
	phase =((x.re==0.0) && (x.im ==0.0)) ?
		0.0 : atan2(x.im,x.re);
	return phase;
}
/*----------------------------------------------*
vector init_vect(long taille ,double valeur)
{
	vector a;
	long i;
	a.vect=(double *)calloc(taille , sizeof(double));
	if (a.vect==NULL) Err_Message("Memory allocation Error in init_vect");
	a.taille = taille;
	for (i=0; i<taille ; i++) a.vect[i]=valeur;
	return a;
}
/*----------------------------------------------*
vector init_vect(long taille)
{
	vector a;
	a.vect=(double *)calloc(taille , sizeof(double));
	if (a.vect==NULL) Err_Message("Memory allocation Error in init_vect");
	a.taille = taille;
	return a;
}


/*----------------------------------------------*
long_vector init_long_vect(long taille ,long valeur=0)
{
	long_vector a;
	long i;
	a.vect=(long *)calloc(taille , sizeof(long));
	if (a.vect==NULL) Err_Message("Memory allocation Error in init_long_vect");
	a.taille = taille;
	for (i=0; i<taille ; i++) a.vect[i]=valeur;
	return a;
}

/*----------------------------------------------*
void liberer(vector a)
{
	free(a.vect);
}

/*----------------------------------------------*
void liberer(long_vector a)
{
	free(a.vect);
}
/*-------------------------------------------------*
cmplx_vect init_cmplx_vect(long taille, cmplx val)
{
	cmplx_vect a;
	long i;
	a.vect=(cmplx *)calloc(taille, sizeof(cmplx));
	if (a.vect==NULL) Err_Message("Memory allocation Error in init_cmplx_vect");
	for (i=0; i<taille ; i++) 
		a.vect[i]=val;
	a.taille=taille;
	return a;
}
/*-------------------------------------------------*
cmplx_vect init_cmplx_vect(long taille)
{
	cmplx_vect a;
	a.vect=(cmplx *)calloc(taille, sizeof(cmplx));
	if (a.vect==NULL) Err_Message("Memory allocation Error in init_cmplx_vect");
	a.taille=taille;
	return a;
}

/*------------------------------------------------*
void liberer(cmplx_vect a)
{
	free(a.vect);
}

 
long_matrice init_long_mat(long LIN, long COL, long valeur)
{
	long_matrice a;
	long i,j;
	a.mat=(long * *)calloc(LIN,sizeof(long *));
	if (a.mat==NULL) Err_Message("Memory allocation Error in init_matrice");
	a.mat[0]=(long *)calloc(LIN*COL, sizeof(long));
	if (a.mat[0]==NULL) Err_Message("Memory allocation Error in init_matrice");
	for (i=1; i<LIN; i++)
		a.mat[i]=COL+a.mat[i-1];
	for (i=0; i<LIN; i++)
		for (j=0 ; j<COL ; j++)
			a.mat[i][j]=valeur;
	a.lin=LIN;
	a.col=COL;
	return a;
}

/*---------------------------------------------*
void liberer(long_matrice a)
{
	free (a.mat[0]);
	free(a.mat);
}

 
matrice init_matrice(long LIN, long COL, double valeur)
{
	matrice a;
	long i,j;
	a.mat=(double * *)calloc(LIN,sizeof(double *));
	if (a.mat==NULL) Err_Message("Memory allocation Error in init_matrice");
	a.mat[0]=(double *)calloc(LIN*COL, sizeof(double));
	if (a.mat[0]==NULL) Err_Message("Memory allocation Error in init_matrice");
	for (i=1; i<LIN; i++)
		a.mat[i]=COL+a.mat[i-1];
	for (i=0; i<LIN; i++)
		for (j=0 ; j<COL ; j++)
			a.mat[i][j]=valeur;
	a.lin=LIN;
	a.col=COL;
	return a;
}
/*---------------------------------------------*
void liberer(matrice a)
{
	free (a.mat[0]);
	free(a.mat);
}

/*-----------------------------------------------------------*
vector col_mat2vect(matrice a, long col)
{
	vector b;
	long i;
	b=init_vect(a.lin,0);
	for (i=0; i<a.lin; i++)
		b.vect[i]=a.mat[i][col];
	return b;
}

/*-------------------------------------------------------------*
cmplx_vect conjug(cmplx_vect a)
{
	long	i;
	cmplx_vect	c;
	c=init_cmplx_vect(a.taille,init_cmplx(0,0));
	for(i=0; i<a.taille; i++)
		c.vect[i]=conjug(a.vect[i]);
	return c;
}

/*-------------------------------------------------------------*/
void Err_Message(char Msg[])
{
	cout << Msg << endl;
	exit (2);
}
/*------------------------------------------------*
void cmplx_s_mul_vect (cmplx_vect a, cmplx_vect b, cmplx c)
{
	long i;
	for (i=0 ; i<a.taille ; i++)
		b.vect[i]=cmplx_mul(a.vect[i],c);
}
/*----------------------------------------------*
void cmplx_mul_vect(cmplx_vect a, cmplx_vect b, cmplx_vect c,double conj=1.0)
{
	long i,size;
	size=min(min(a.taille,b.taille),c.taille);
	if (conj==1)
		for (i=0 ; i<size ; i++)
			c.vect[i]=cmplx_mul(a.vect[i],b.vect[i]);
	else
		for (i=0 ; i<size ; i++)
			c.vect[i]=cmplx_mul(a.vect[i],conjug(b.vect[i]));
}
/*----------------------------------------------*
void sub_matrice(matrice a, matrice b, matrice c)
{
	long i,j;
	if(a.lin!=b.lin || a.col!=b.col || a.lin!=c.lin || a.col != c.col)
		Err_Message("taille incompatible dans sub_matrice");
	for (i=0; i<a.lin; i++)
		for (j=0; j<a.col; j++)
			c.mat[i][j]=a.mat[i][j]-c.mat[i][j];
}
/*----------------------------------------------*

void print (cmplx a)
{
	if (a.im < 0)
		cout <<a.re<<"-j"<<-a.im<<endl;
	else
		cout <<a.re<<"+j"<<a.im<<endl;
}
/*----------------------------------------------*
void print(cmplx_vect a)
{
	long i;
	cout << endl;
	for (i=0 ; i<a.taille ; i++)
		print(a.vect[i]);
	cout << endl;
}

/*---------------------------------------------*
void print(vector a)
{
	long i;
	cout << endl;
	for (i=0; i<a.taille ; i++) cout << "   " << a.vect[i] << endl;
	cout << endl;
}

/*---------------------------------------------*
void print(long_vector a)
{
	long i;
	cout << endl;
	for (i=0; i<a.taille ; i++) cout << "  " << a.vect[i] << endl;
	cout << endl;
}
/*----------------------------------------------*
void print(matrice a)
{
	long i,j;
	cout << endl;
	for (i=0; i<a.lin ; i++){
		for (j=0; j<a.col ; j++)
			printf("%10f ",a.mat[i][j]);
		cout << endl;
	}
}
/*----------------------------------------------*
void print(long_matrice a)
{
	long i,j;
	cout << endl;
	for (i=0; i<a.lin ; i++){
		for (j=0; j<a.col ; j++)
			printf("%5d ",a.mat[i][j]);
		cout << endl;
	}
}

/*----------------------------------------------*/
double n_gaussien( double moy , double var ) //___________LOI GAUSSIENNE
//	retourne un evenement suivant une loi gaussienne de moyenne moy
//	et de variance var
{
	static long iset;
	static double gset;
	double fac,rsq,v1,v2;
	const double RAND_MAX_2=2.0/RAND_MAX;

	if (iset==0){
		do{
			v1=rand()*RAND_MAX_2-1.0;
			v2=rand()*RAND_MAX_2-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq ==0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac*sqrt(var)+moy;
	} else {
		iset=0;
		return gset*sqrt(var)+moy;
	}
}
/*----------------------------------------------*/
/* Generation d'un vecteur de bits aleatoires  */
void random_bit(long_vector a)
{
	long i;
	for (i=0; i<a.taille; i++)
		a.vect[i]=rand() & 1;
}

/*-------------------------------------------
*   Génération d'un vecteur de symbole aléatoires entre 0 et M-1
*   Attention : M doit etre 2**k
*/
void rand_M(long_vector a, long M)
{
	long i;
	for (i=0; i<a.taille; i++)
		a.vect[i]=rand() & M-1;
}
/*-------------------------------------------------------*/
/*  Modulation QPSK shifted, codage Gray  */
void QPSK_mod(long_vector data,cmplx_vect x)
{
	long i;
	double amp;
	amp=sqrt(2.0)/2.0;
	if (x.taille*2<data.taille)
		Err_Message("incompatible sizes in QPSK_MOD function");
	for (i=0; i<data.taille/2; i++){
		x.vect[i].re=(data.vect[2*i]==0) ? amp : -amp;
		x.vect[i].im=(data.vect[2*i+1]==0) ? amp : -amp;
	}
}
/*-------------------------------------------------------*/
/*  Demodulation QPSK shifted, codage Gray  */
void QPSK_demod(cmplx_vect x, long_vector data_hat)
{
	long i;
	for (i=0; i<data_hat.taille/2; i++){
		data_hat.vect[2*i] = (x.vect[i].re >= 0) ? 0 : 1;
		data_hat.vect[2*i+1] = (x.vect[i].im >= 0) ? 0 : 1;
	}
}
/*-------------------------------------------------------*/
void convol (vector x, vector h, vector y)
{
  long m;
  long j;
  double sum;

  /* Loop over output points */
  for (m = y.taille-1; m >= 0 ; m--) {

	 /* Convolution */
	 sum = 0.0;
	 for (j = 0; j < h.taille; ++j) {
		if (m-j>=0 && m-j<x.taille) sum += h.vect[j] * x.vect[m-j];
	 }
	 y.vect[m] = sum;
  }
}
/*-------------------------------------------------------*/
/* convolution : entree et sortie complexes, filtre réel  */
/* La sortie peut etre le meme vector que l'entrée */
void convol (cmplx_vect x, vector h, cmplx_vect y)
{
  long m;
  long j;
  double sum_R, sum_I;

  /* Loop over output points */
  for (m = y.taille-1; m >= 0 ; m--) {

	 /* Convolution */
	 sum_R = sum_I = 0.0;
	 for (j = 0; j < h.taille; ++j) {
		if (m-j>=0 && m-j<x.taille){
			sum_R += h.vect[j] * x.vect[m-j].re;
			sum_I += h.vect[j] * x.vect[m-j].im;
		}
	 }
	 y.vect[m].re = sum_R;
	 y.vect[m].im = sum_I;
  }
}
/*-------------------------------------------------------*/
/* convolution : entree, sortie et filtre complexes  */
/* La sortie peut etre le meme vector que l'entrée */
void convol (cmplx_vect x, cmplx_vect h, cmplx_vect y)
{
  long m;
  long j;
  cmplx sum(0.0,0.0);

  /* Loop over output points */
  for (m = y.taille-1; m >= 0 ; m--) {

	 /* Convolution */
	 sum = 0.0; //init_cmplx(0,0);
	 for (j = 0; j < h.taille; ++j)
		 if (m-j>=0 && m-j<x.taille)
//			sum = cmplx_add(sum,cmplx_mul(h.vect[j] , x.vect[m-j]) );
			sum = sum+h.vect[j]*x.vect[m-j];
	 y.vect[m] = sum;
  }
}

/*------------------------------------------------------*/
cmplx_vect	correlation(cmplx_vect a, cmplx_vect b)
{
//	cmplx_vect R;
	long	n,m;
	cmplx	sum;
	//R=init_cmplx_vect(a.taille+b.taille-1, init_cmplx(0,0));
	cmplx_vect R(a.taille+b.taille-1);;
		
	for (m=b.taille-1 ; m>-a.taille ; m--){
		sum=0.0;//init_cmplx(0,0);
		for (n=0;n<a.taille; n++)
			if (m+n>=0 && m+n<b.taille)
				//sum=cmplx_add(sum,cmplx_mul(a.vect[n],conjug(b.vect[m+n])));
				sum+=a.vect[n]*b.vect[m+n].conj();
		R.vect[b.taille-1-m]=sum;
	}
	return R;
}
/*------------------------------------------------------*/
vector File2Vect (char Fname[])
{
  FILE *fp;
  long Nval=0;
  double tmp;
  
  double *xp;
/* Open the data file */
  fp = fopen (Fname, "r");
  if (fp == NULL) {
	 cout << "File ReadData: Cannot open file "<< Fname<<endl;
	 exit (1);
  }
/* Read the values */
  while (!feof(fp)){
	fscanf(fp,"%lf\n",&tmp);
	Nval++;
  }
  fseek(fp,0L,SEEK_SET);
  vector x(Nval,0);
  xp=x.vect;
  while (!feof(fp))
		fscanf(fp,"%lf\n",xp++);
/* Close the file */
  fclose (fp);
  x.taille=Nval;
  return x;
}
/*------------------------------------------------------*/
long_vector File2long_vect (char Fname[])
{
  FILE *fp;
  int Nval=0;
  int  tmp;
  int  *xp;
/* Open the data file */
  fp = fopen (Fname, "r");
  if (fp == NULL) {
	 cout << "File ReadData: Cannot open file "<< Fname<<endl;
	 exit (1);
  }
/* Read the values */
  while (!feof(fp)){
	fscanf(fp,"%d\n",&tmp);
	Nval++;
  }
  fseek(fp,0L,SEEK_SET);
  long_vector x(Nval,0);
  xp=x.vect;
  while (!feof(fp))
		fscanf(fp,"%d\n",xp++);
/* Close the file */
  fclose (fp);
  x.taille=Nval;
  return x;
}
/*------------------------------------------------------*/
void File2long_matrice (char Fname[],long_matrice a)
{
  FILE *fp;
  long Nval=0,i,j;
  long tmp;
/* Open the data file */
  fp = fopen (Fname, "r");
  if (fp == NULL) {
	 cout << "File ReadData: Cannot open file "<< Fname<<endl;
	 exit (1);
  }
/* Read the values */
  while (!feof(fp)){
	fscanf(fp,"%d\n",&tmp);
	Nval++;
  }
  if (Nval != a.col*a.line) Err_Message("File size incompatible!");
  fseek(fp,0L,SEEK_SET);
  for (i=0; i<a.line; i++)
	  for (j=0; j<a.col ; j++)
		  fscanf(fp,"%d\n",&a.mat[i][j]);
/* Close the file */
  fclose (fp);
}
/*------------------------------------------------------*/
cmplx_vect File2Vect_cmplx (char Fname[])
{
  FILE *fp;
  long Nval=0;
  double temp1, temp2;
  cmplx *p;

/* Open the data file */
  fp = fopen (Fname, "r");
  if (fp == NULL) {
	 cout << "File ReadData: Cannot open file "<< Fname<<endl;
	 exit (1);
  }
/* Read the values */
  while (!feof(fp)){
	fscanf(fp,"%lf %lf\n",&temp1, &temp2);
	Nval++;
  }
  fseek(fp,0L,SEEK_SET);
  cmplx_vect x(Nval);
  p=x.vect;
  while (!feof(fp)){
	fscanf(fp,"%lf %lf\n",&p->re, &p->im);
	p++;
  }
/* Close the file */
  fclose (fp);
  x.taille=Nval;
  return x;
}
/*-----------------------------------------------------------*/
void File_Print (char Fname[],vector x)
{
	FILE *fp;
	long i;
	fp=fopen (Fname, "w");
	for (i=0; i<x.taille ; i++)
		fprintf(fp,"%lf\n",*x.vect++);
	fclose(fp);
}
/*-----------------------------------------------------------*/
void File_Print (char Fname[],long_vector x)
{
	FILE *fp;
	long i;
	fp=fopen (Fname, "w");
	for (i=0; i<x.taille ; i++)
		fprintf(fp,"%ld\n",*x.vect++);
	fclose(fp);
}
/* -------------------------------------------------------*/
/* ecrire des données complexes dans  un fichier  */
long File_Print(char nom[],cmplx_vect a)
{
	FILE *fp;
	long i;
	if ((fp = fopen(nom,"w")) ==NULL) return 1;
	else {
		for (i=0; i<a.taille; i++)
			fprintf(fp,"%lf  %lf\n",a.vect[i].re, a.vect[i].im);
		fclose(fp);
		return 0;
	}
}
/*--------------------------------------------------------*/
void Fill_vect(vector a, double val)
{
	long i;
	for (i=0; i<a.taille ; i++)
		*a.vect++ = val;
}
/*--------------------------------------------------------*/
/* Fill a comlex vector */
void Fill_cmplx_vect(cmplx_vect a, cmplx val)
{
	long i;
	for (i=0; i<a.taille ; i++)
		*a.vect++ = val;
}
/*--------------------------------------------------------*/
void zero_pad(vector x, vector y, long rate)
{
	long i;
	if (y.taille<x.taille*rate)
		Err_Message("Error in zero_pad function: incompatible output vector size");
	Fill_vect(y,0);
	for (i=0; i<y.taille ; i++)
		y.vect[i*rate]=x.vect[i];
}
/*----------------------------------------------------*/
/* zero padding for complex signals  */
void zero_pad(cmplx_vect x, cmplx_vect y, long rate)
{
	long i;
	if (y.taille<x.taille*rate)
		Err_Message("Error in zero_pad function: incompatible output vector size");
	Fill_cmplx_vect(y, DCplx(0.0,0.0));
	for (i=0; i<x.taille ; i++)
		y.vect[i*rate]=x.vect[i];
}
/*--------------------------------------------------------*/
void decimate (vector x, vector y, long delai,long rate)
{
	long i,ind;
	for (i=0 ; i<y.taille ; i ++){
		ind=i*rate+delai;
		if (ind<x.taille)
			*y.vect++ = x.vect[ind];
		else *y.vect++ = 0;
	}
}
/*--------------------------------------------------------*/
void decimate (cmplx_vect x, cmplx_vect y, long delai,long rate)
{
	long i,ind;
	for (i=0 ; i<y.taille ; i ++){
		ind=i*rate+delai;
		if (ind<x.taille)
			*y.vect++ = x.vect[i*rate+delai];
		else
			*y.vect++ = 0.0;
	}
}
/*---------------------------------------------------------*/
long indice_min(vector a)
{
	long i, indice;
	double mini;
	mini=a.vect[0];
	indice=0;
	for (i=1 ; i<a.taille ; i++)
		if (a.vect[i]<mini) {
			mini=a.vect[i];
			indice=i;
		}
	return indice;
}
/*---------------------------------------------------------*/
long indice_max(vector a)
{
	long i, indice;
	double maxi;
	maxi=a.vect[0];
	indice=0;
	for (i=1 ; i<a.taille ; i++)
		if (a.vect[i]>maxi) {
			maxi=a.vect[i];
			indice=i;
		}
	return indice;
}
/*-----------------------------------------------------------*/
void shift_right(long_vector a, long decal)
{
	int *des,*src;
	long i;
	des=&a.vect[a.taille-1];
	src=&a.vect[a.taille-decal-1];
	for (i=0;i<a.taille-decal;i++)
		*des-- = *src--;
	for (i=0; i<decal;i++)
		*des-- = 0;
}
/*-----------------------------------------------------------*/
void shift_right(vector a, long decal)
{
	double *des,*src;
	long i;
	des=&a.vect[a.taille-1];
	src=&a.vect[a.taille-decal-1];
	for (i=0;i<a.taille-decal;i++)
		*des-- = *src--;
	for (i=0; i<decal;i++)
		*des-- = 0;
}
/*-----------------------------------------------------------*/
void shift_right(cmplx_vect a, long decal)
{
	cmplx *des,*src;
	long i;
	des=&a.vect[a.taille-1];
	src=&a.vect[a.taille-decal-1];

	for (i=0;i<a.taille-decal;i++)
		*des-- = *src--;
	for (i=0; i<decal;i++)
		*des-- = 0.0;//init_cmplx(0,0);
}
/*-----------------------------------------------------------*/
void shift_left(long_vector a,long decal)
{
	int *des,*src,i;
	des=&a.vect[0];
	src=&a.vect[decal];
	for (i=0;i<a.taille-decal;i++)
		*des++ = *src++;
	for (i=0; i<decal;i++)
		*des++ = 0;
}
/*-----------------------------------------------------------*/
void shift_left(vector a,long decal)
{
	double *des,*src;
	long i;
	des=&a.vect[0];
	src=&a.vect[decal];
	for (i=0;i<a.taille-decal;i++)
		*des++ = *src++;
	for (i=0; i<decal;i++)
		*des++ = 0;
}
/*-----------------------------------------------------------*/
void shift_left(cmplx_vect a,long decal)
{
	cmplx *des,*src;
	long i;
	des=a.vect;
	src=&a.vect[decal];

	for (i=0;i<a.taille-decal;i++)
		*des++ = *src++;
	for (i=0; i<decal;i++)
		*des++ = 0.0;//init_cmplx(0,0);
}
/*------------------------------------*/
void lin_vect(vector a,double initial, double step)
{
	long i;
	a.vect[0]=initial;
	for (i=1; i<a.taille; i++)
		a.vect[i] = a.vect[i-1] + step;
}
/*-----------------------------------------*/
/* cette routine calcule exp(j*input)  */
cmplx_vect cmplx_exp(vector a)
{
	long i;
	cmplx_vect x(a.taille);
	for (i=0 ; i<a.taille ; i++){
		x.vect[i].re=cos(*a.vect);
		x.vect[i].im=sin(*a.vect++);
	}
	x.taille=a.taille;
	return x;
}
/*-----------------------------------------*/
vector copy(vector a)
{
	
	long i;
	vector b(a.taille);
	b.taille=a.taille;
	for (i=0; i<a.taille; i++)
		b.vect[i]=a.vect[i];
	return b;
}
/*-----------------------------------------*/
cmplx_vect copy(cmplx_vect a)
{
	
	long i;
	cmplx_vect b(a.taille);
	b.taille=a.taille;
	for (i=0; i<a.taille; i++)
		b.vect[i]=a.vect[i];
	return b;
}
/*-----------------------------------------*/
long_vector copy(long_vector a)
{
	long i;
	long_vector b(a.taille);
	b.taille=a.taille;
	for (i=0; i<a.taille; i++)
		b.vect[i]=a.vect[i];
	return b;
}
/*-----------------------------------------*/
matrice copy(matrice a)
{
	
	double *p_a,*p_b;
	long i;
	p_a=&a.mat[0][0];
	matrice b(a.line,a.col);
	p_b=&b.mat[0][0];
	b.col=a.col;
	b.line=a.line;
	for (i=0; i<a.line*a.col; i++)
		*p_b++ = *p_a++;
	return b;
}
/*-----------------------------------------*/
long_matrice copy(long_matrice a)
{

	int *p_a,*p_b;
	long i;
	p_a=&a.mat[0][0];
	long_matrice b(a.line,a.col);
	p_b=&b.mat[0][0];
	b.col=a.col;
	b.line=a.line;
	for (i=0; i<a.line*a.col; i++)
		*p_b++ = *p_a++;
	return b;
}
/*-----------------------------------------*/
/* cette fonction présume que la puissance du signal en entrée est 
egale à un */
/*
/*This model adds discrete-time zero-mean white gaussian noise of variance
   FACTOR * 10**(-SNR/10) * (Signal_power)
to the input signal, where Signal_ power is 1.
Note:
The input signal of AWGN must be real,
	FACTOR = Ns/(2*K)
	SNR    = 10*log(Eb/N0)
where Ns   = number of samples per symbol
	K    = T/Tb = (transmitted symbol duration) / (information bit duration)
	Eb   = input signal energy per bit of the original
			information-carrying bit stream of rate 1/Tb
	N0/2 = power spectral density of the continuous-
			time Gaussian noise which is assumed to have
			a doublesided bandwidth equal to the sampling rate
*/
void awgn(vector a, double SNR, double factor)
{
	double var;
	long i;
	var=factor*pow(10.0,-SNR/10.0);
	for (i=0 ; i<a.taille ; i++)
		a.vect[i] += n_gaussien(0,var);
}

/*-----------------------------------------*/
/* cette fonction présume que la puissance du signal en entrée est 
egale à un */
/*
/*This model adds discrete-time zero-mean white gaussian noise of variance
   FACTOR * 10**(-SNR/10) * (Signal_power)
to the input signal, where Signal_ power is 1.
Note:
The input signal of AWGN must be real,
	FACTOR = Ns/2K
	SNR    = 10*log(Eb/N0)
where Ns   = number of samples per symbol
	K    = T/Tb = (transmitted symbol duration) / (information bit duration)
	Eb   = input signal energy per bit of the original
			information-carrying bit stream of rate 1/Tb
	N0/2 = power spectral density of the continuous-
			time Gaussian noise which is assumed to have
			a doublesided bandwidth equal to the sampling rate
*/

void awgn(cmplx_vect a, double SNR, double factor)
{
	double var;
	long i;
	var=factor*pow(10.0,-SNR/10.0);
	for (i=0 ; i<a.taille ; i++){
		a.vect[i].re += n_gaussien(0,var);
		a.vect[i].im += n_gaussien(0,var);
	}
}

/*----------------------------------------*/
vector BPSK_MOD(long_vector a)
{
	
	long i;
	vector b(a.taille);
	for (i=0 ; i<a.taille ; i++)
		b.vect[i]=((a.vect[i]==1) ? 1.0 : -1.0);
	return b;
}
/*----------------------------------------*/
long_vector	BPSK_DEM(vector a)
{
	
   long i;
	long_vector b(a.taille);
	for (i=0 ; i<a.taille ; i++)
		b.vect[i]=((a.vect[i]>0) ? 1 : 0);
	return b;
}
/*----------------------------------------*/
long	comparer(long_vector a, long_vector b)
{
	long cmpt=0,i;
	int *p_a, *p_b;
	p_a=a.vect;
	p_b=b.vect;
	for (i=0 ; i<min(a.taille, b.taille) ; i++)
		if (*p_a++ != *p_b++) cmpt++;
	return cmpt;
}


