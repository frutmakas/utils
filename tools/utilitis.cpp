/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/07/28 11:03:32 $
 * $Revision: 1.1.2.72 $
 * $Id: utilitis.cpp,v 1.1.2.72 2006/07/28 11:03:32 syed Exp $
 ********************************************************************/

/*#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
*/
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

#include "tools/utilitis.h"

#include "tools/tools.h"

#include "globaldef.h"

#ifdef WIN_DOS
#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 
#pragma intrinsic (atan, exp, log10, sqrt, atan2, log, sin, tan, cos, fabs, abs ,pow)
#endif
//#error("Unfinished .. ");


// cmplx        init_cmplx(double a, double b);
DCplx::DCplx(const double re, const double im){this->re=re; this->im=im;}
DCplx::DCplx(){this->re=0.0; this->im=0.0;}
DCplx::DCplx(const DCplx& c){re=c.re; im=c.im;}

DCplx &DCplx::operator=(const DCplx &c) { re=c.re; im=c.im ; return *this; }
DCplx &DCplx::operator=(const double c) { re=c;im=0.0; return *this; }

bool DCplx::operator==(const DCplx &c) const{ return (re==c.re)&&(im==c.im); }
bool DCplx::operator==(const double c) const{ return (re==c)&&(im==0.0); }

bool DCplx::operator!=(const DCplx &c) const{ return !((re==c.re)&&(im==c.im)); }
bool DCplx::operator!=(const double c) const{ return !((re==c)&&(im==0.0)); }

DCplx &DCplx::operator+=(const DCplx &c) {
	re+=c.re;
	im+=c.im;
	return (*this); 
}

DCplx &DCplx::operator+=(const double c) {
	re+=c;
	return (*this); 
}

DCplx &DCplx::operator-=(const DCplx &c) {
	re-=c.re;
	im-=c.im;
	return (*this); 
}

DCplx &DCplx::operator-=(const double c) {
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


DCplx log(DCplx c) {

	return DCplx(log(c.mag()), c.phase());

}

WVector::WVector(const int &size){
	if (size<0) throw CUtilitisException(__FILE__, __LINE__, INVALID_VECTOR_DIMENSION);
    if(size==0) { vect=NULL; taille=0; return; }
	vect=new int[size];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	taille=size;
	for(register int i=0;i<size;i++) vect[i]=0;
}

WVector::WVector(const int &size, const int &initial) {
	if (size<0) throw CUtilitisException(__FILE__, __LINE__, INVALID_VECTOR_DIMENSION);
    if(size==0) { vect=NULL; taille=0; return; }
	vect=new int[size];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	taille=size;
	for(register int i=0;i<size;i++) vect[i]=initial;
}

WVector::WVector(const int &size, const int &initial, const char &mode) {
	if (size<0) throw CUtilitisException(__FILE__, __LINE__, INVALID_VECTOR_DIMENSION);
    if(size==0) { vect=NULL; taille=0; return; }
	vect=new int[size];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	taille=size;
    int init = initial;

    switch(mode) {
        case 'i' : // incremental value
                	for(register int i=0;i<size;i++) vect[i]=i;
                    break;
        case 'd' : // decremental value
                	for(register int i=0,p=size-1;i<size;i++,p--) vect[i]=p;
                    break;
        case 'I' : // incremental value from initial
                	for(register int i=initial, p=0 ;p<size;p++,i++) vect[p]=i;
                    break;
        case 'D' : // decremental value from initial
                	for(register int i=initial+size-1, p=0;p<size;p++,i--) vect[p]=i;
                    break;
        case 'p' : // generate prime number
                    register int i,x,y, l;    bool isprime;
                    vect[0]=1; if(size==1) break;
                    vect[1]=2; if(size==2) break;
                    vect[2]=3; if(size==3) break;

                    for(i=5,x=3;x<size;i+=2) {
                        l = (int) sqrt(i)+1;
                        isprime=true;
                        for(y=1;y<x && vect[y]<=l;y++) {
                            if(i%vect[y]==0){
                                isprime=false;
                                break;
                            }
                        }
                        if(isprime) {
                            vect[x++]=i;
                        }
                    }
                    break;
        case 'o' : // odd value
                	for(register int i=0,p=1;i<size;i++,p+=2) vect[i]=p;
                    break;
        case 'O' : // odd value starting from initial
                    if(initial&1==0) init++;
                	for(register int i=init, p=0 ;p<size;p++,i++) vect[p]=i;
                    break;
        case 'e' : // even value
                	for(register int i=0,p=0;i<size;i++,p+=2) vect[i]=p;
                    break;
        case 'E' : // odd value starting from initial
                    if(initial&1) init++;
                	for(register int i=init, p=0 ;p<size;p++,i++) vect[p]=i;
                    break;
        case 'b' : // fill with random bit randbit1 algorithm
                    for(register int i=0, tseed=initial ? initial : time(NULL);i<size;i++) {
                    	vect[i] = (tseed >> 17) & 1 ^ (tseed >> 4) &1 ^ (tseed >> 4) & 1 ^ (tseed &1);
	                    tseed = (tseed<<1)| vect[i];
                    }
                    break;
        case 'f' : // fibonacci numbers
                    vect[0]=1; if(size==1) break;
                    vect[1]=1; if(size==2) break;
                	for(register int i=2; i<size;i++) vect[i]=vect[i-1]+vect[i-2];
                    break;
        case 'c' : // constant value
        default  :
                	for(register int i=0;i<size;i++) vect[i]=initial;
                    break;
    }
}

WVector::WVector(const int *v, const int size) {
	if (size<0) throw CUtilitisException(__FILE__, __LINE__, INVALID_VECTOR_DIMENSION);
    if(size==0) { vect=NULL; taille=0; return; }
    if (v==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	vect=new int[size];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	taille=size;
	for(register int i=0;i<size;i++) vect[i]=v[i];
}

WVector::WVector(const WVector &c){
    if(c.taille<0)  throw CUtilitisException(__FILE__, __LINE__, INVALID_VECTOR_DIMENSION);
	if(c.taille==0) {
		taille=0; vect=NULL;
		return;
	}
	taille=c.taille;
	vect=new int[taille];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
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
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL)
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille)
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	WVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c.vect[i];
	}
	return n;
}

DVector WVector::operator+(const DVector &c) {
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL)
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille)
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
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

DVector WVector::operator+(const double c) { 
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
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL)
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=c.vect[i]+vect[i];
	}
	return n;
}

// cmplx        cmplx_sub(cmplx a,cmplx b);
WVector WVector::operator-(const WVector &c) { 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	WVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c.vect[i];
	}
	return n;
}

WVector WVector::operator-() { 
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	WVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=-vect[i];//-c.vect[i];
	}
	return n;
}

DVector WVector::operator-(const DVector &c) { 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
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

DVector WVector::operator-(const double c) { 
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-(double)c;
	}
	return n;
}

ZVector WVector::operator-(const ZVector &c) { 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c.vect[i];
	}
	return n;
}

//cmplx        cmplx_mul(cmplx a,cmplx b);
ZVector WVector::operator*(const ZVector &c) { 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
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

DVector WVector::operator*(const double c) { 
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c;
	}
	return n;
}


WVector WVector::operator*(const WVector &c) { 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	WVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c.vect[i];
	}
	return n;
}


DVector WVector::operator*(const DVector &c) { 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c.vect[i];
	}
	return n;
}

DVector WVector::operator/(const double c) { 
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c==0.0) 
		throw CUtilitisException(__FILE__, __LINE__, DIVISION_BY_ZERO);
	double d=1/c;
	return DVector((*this)*d);
}

ZVector WVector::operator/(const DCplx &c) {
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);

	double t=c.mod2();

	if (t==0.0) 
		throw CUtilitisException(__FILE__, __LINE__, DIVISION_BY_ZERO);

	DCplx ct=c.conj()/t;
	return ZVector((*this)*ct);
}


WVector &WVector::operator=(const WVector &c) { 
    if(taille==c.taille && vect!=NULL) {
        for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
        return *this;
    }
	if(taille!=0 && vect!=NULL) delete[] vect;
	taille=c.taille;
	if(c.taille>0) vect = new int[taille];
	if(c.taille>0 && vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
	return *this;
}

/* * *  accelerated code above  * *

WVector &WVector::operator=(const WVector &c) { 
	if(taille!=0 && vect!=NULL) delete[] vect;
	taille=c.taille;
	vect = new int[taille];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
	return *this;
}*/

WVector &WVector::operator=(int d) { 
	throw CUtilitisException(__FILE__, __LINE__, INVALID_OPERATION);
}

DVector &DVector::operator=(int d) { 
	throw CUtilitisException(__FILE__, __LINE__, INVALID_OPERATION);
}

ZVector &ZVector::operator=(int d) { 
	throw CUtilitisException(__FILE__, __LINE__, INVALID_OPERATION);
}

bool WVector::operator==(const DVector &c){
	if(taille!=c.taille) return false;
	register int i=0;
	while(i<taille && vect[i]==c.vect[i]) i++;
	return (i==taille);
}

bool WVector::operator==(const WVector &c) const{
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

bool WVector::operator!=(const WVector &c) const {
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
	if(rhs.taille<=0 || rhs.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (rhs.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	register DCplx sum(0.0,0.0);
	for(register int i=0;i<taille;i++) {
		sum+=rhs.vect[i]*vect[i];
	}
	return sum;
}

WVector WVector::copy(int start, int stop) const{
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(start>stop) swap(start, stop);
	if (stop>taille) stop=taille; 
	if (start<0) start=0;
	if (!(stop-start)) return WVector();
	WVector n(stop-start);
	for(register int i=start, j=0;i<stop;i++,j++)
		n.vect[j]=vect[i];
	return n;
}

WVector WVector::rev_copy(int start, int stop) const{
	if(taille <=0 || vect==NULL)
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(stop>start) swap(stop, start);
	if (start>taille) start=taille; 
	if (stop<0) stop=0;
	if (stop==start) return WVector();
	WVector n(start-stop);
	for(register int i=start, j=0;i>stop;i--,j++)
		n.vect[j]=vect[i];
	return n;
}

void WVector::reset(int value){
	if(vect==NULL ||taille==0) return;
	for(register int i=0;i<taille;i++) vect[i]=value;
}

void WVector::insert(const WVector &rhs, int pos) const{
	if(rhs.taille<=0 || rhs.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(pos<0) 
		throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	for(register int i=0, j=pos; i < rhs.taille && j<taille;i++, j++)
		vect[j]=rhs.vect[i];
}

WVector WVector::RZtoNRZ_nip(void){
	WVector nv(taille);
	for(register int i=0;i<taille;i++){
		nv.vect[i]=vect[i]==1?1:-1;
	}
	return nv;
}

void WVector::RZtoNRZ(void){
	for(register int i=0;i<taille;i++){
		vect[i]=vect[i]==1?1:-1;
	}
}

DVector::DVector(const int &size,  double initial) {
	if (size<0) throw CUtilitisException(__FILE__, __LINE__, INVALID_VECTOR_DIMENSION);
    if(size==0) { vect=NULL; taille=0; return; }
	vect=new double[size];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	taille=size;
	for(register int i=0;i<size;i++) vect[i]=initial;
}


DVector::DVector(const double *v, const int size) {
	if (size<0) throw CUtilitisException(__FILE__, __LINE__, INVALID_VECTOR_DIMENSION);
    if(size==0) { vect=NULL; taille=0; return; }
	vect=new double[size];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
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
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i=0;i<taille;i++) vect[i]=(double)c.vect[i];
}

DVector::DVector(const DVector &c){
	if(c.taille==0) {
		taille=0; vect=NULL;
		return;
	}
	taille=c.taille;
	vect=new double[taille];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
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
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL)
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille)
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c.vect[i];
	}
	return n;
}

DVector DVector::operator+(const WVector &c)  const{ 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+(double)c.vect[i];
	}
	return n;
}

ZVector DVector::operator+(const ZVector &c)  const{ 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
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

DVector DVector::operator+(const double c)  const{ 
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c;
	}
	return n;
}


DVector DVector::operator-() { 
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=-vect[i];//-c.vect[i];
	}
	return n;
}

// cmplx        cmplx_sub(cmplx a,cmplx b);
DVector DVector::operator-(const DVector &c)  const{ 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c.vect[i];
	}
	return n;
}

DVector DVector::operator-(const WVector &c)  const{ 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-(double)c.vect[i];
	}
	return n;
}

ZVector DVector::operator-(const ZVector &c)  const{
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
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

DVector DVector::operator-(const double c)  const{ 
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c;
	}
	return n;
}

//cmplx        cmplx_mul(cmplx a,cmplx b);
ZVector DVector::operator*(const ZVector &c)  const{ 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c.vect[i];
	}
	return n;
}

DVector DVector::operator*(const DVector &c)  const{ 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c.vect[i];
	}
	return n;
}

DVector DVector::operator*(const WVector &c)  const{ 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*(double)c.vect[i];
	}
	return n;
}

ZVector DVector::operator*(const DCplx &c)  const{ 
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c;
	}
	return n;
}

DVector DVector::operator*(const double c)  const{ 
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c;
	}
	return n;
}

DVector DVector::operator/(const double c)  const{ 
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c==0.0) 
		throw CUtilitisException(__FILE__, __LINE__, DIVISION_BY_ZERO);
	double d=1/c;
	return DVector((*this)*d);
}

ZVector DVector::operator/(const DCplx &c)  const{
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	double t=c.mod2();
	if (t==0.0) 
		throw CUtilitisException(__FILE__, __LINE__, DIVISION_BY_ZERO);

	DCplx ct=c.conj()/t;
	return ZVector((*this)*ct);
}

DVector &DVector::operator=(const DVector &c){ 
    if(taille==c.taille && vect!=NULL) {
        for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
        return *this;
    }
	if(taille>0 && vect!=NULL) delete[] vect;
	taille=c.taille;
    vect=NULL;
	if(c.taille>0) if((vect = new double[taille])==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
	return *this;
}

DVector &DVector::operator=(const WVector &c){ 
    if(taille==c.taille && vect!=NULL) {
        for(register int i=0;i<taille;i++) vect[i]=(double)c.vect[i];
        return *this;
    }
	if(taille!=0 && vect!=NULL) delete[] vect;
	taille=c.taille;
	vect = new double[taille];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i=0;i<taille;i++) vect[i]=(double)c.vect[i];
	return *this;
}

DVector& DVector::operator +=(const DVector &dv){
	if (dv.taille != taille || dv.vect==NULL || vect==NULL) return (*this);
	for (register int i=0; i < taille; i++) {
		vect[i]+=dv.vect[i];
	}
	return (*this);
}

DVector& DVector::operator +=(const WVector &wv){
	if (wv.taille != taille || wv.vect==NULL || vect==NULL) return (*this);
	for (register int i=0; i < taille; i++) {
		vect[i]+=wv.vect[i];
	}
	return (*this);
}

DVector& DVector::operator +=(double c){
	if (taille == 0 || vect==NULL) return (*this);
	for (register int i=0; i < taille; i++) {
		vect[i]+=c;
	}
	return (*this);
}

DVector& DVector::operator *=(double c){
	if (taille == 0 || vect==NULL) return (*this);
	for (register int i=0; i < taille; i++) {
		vect[i]*=c;
	}
	return (*this);
}

DVector& DVector::operator -=(const DVector &dv){
	if (dv.taille != taille || dv.vect==NULL || vect==NULL) return (*this);
	for (register int i=0; i < taille; i++) {
		vect[i]-=dv.vect[i];
	}
	return (*this);
}

DVector& DVector::operator -=(const WVector &wv){
	if (wv.taille != taille || wv.vect==NULL || vect==NULL) return (*this);
	for (register int i=0; i < taille; i++) {
		vect[i]-=wv.vect[i];
	}
	return (*this);
}

DVector& DVector::operator -=(double c){
	if (taille == 0 || vect==NULL) return (*this);
	for (register int i=0; i < taille; i++) {
		vect[i]-=c;
	}
	return (*this);
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
	if(rhs.taille<=0 || rhs.vect==NULL || taille <=0 || vect==NULL)
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (rhs.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	register DCplx sum(0.0,0.0);
	for(register int i=0;i<taille;i++) {
		sum+=rhs.vect[i]*vect[i];
	}
	return sum;
}

void DVector::fixedconv(const DVector &rhs) const{
	if(rhs.taille<=0 || rhs.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	
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
	if(rhs.taille<=0 || rhs.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);

	register double sum;
	register int cache;
	int imax=taille+rhs.taille-1;
	double *tmp=new double[imax];
	if(tmp==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);

	/* Loop over output points */
	for (register int i = 0; i<imax ; i++) {
		/* Convolution */
		sum = 0.0;
		for (register int j = 0; j < rhs.taille; ++j) {
			cache=i-j;
			if (cache>=0 && cache<taille) 
				sum += rhs.vect[j] * vect[cache];
		}
		tmp[i] = sum;
	}
	delete[] vect;
	taille=imax;
	vect=tmp;
}

DVector DVector::conv_nip(const DVector &rhs) const{
	if(rhs.taille<=0 || rhs.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
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
			if (cache>=0 && cache<taille) 
				sum += rhs.vect[j] * vect[cache];
		}
		n.vect[i] = sum;
	}
	return n;
}

ZVector DVector::conv_nip(const ZVector &rhs) const{
	if(rhs.taille<=0 || rhs.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);

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
			if (cache>=0 && cache<taille)
				sum += rhs.vect[j] * vect[cache];
		}
		n.vect[i] = sum;
	}
	return n;
}

DVector DVector::copy(int start, int stop) const {
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(start>stop) swap(start,stop);
	if (stop>taille) stop=taille;
	if (start<0) start=0;
    if(start>taille) throw CUtilitisException(__FILE__,__LINE__,INDEX_OUT_OF_RANGE);
	if ((stop-start)==0) return DVector();
	DVector n(stop-start);
	for(register int i=start, j=0;i<stop;i++,j++)
		n.vect[j]=vect[i];
	return n;
}

DVector DVector::rev_copy(int start, int stop) const{
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(stop>start) swap(stop, start);
	if (start>taille) start=taille; 
	if (stop<0) stop=0;
	if (stop==start) return DVector();
	DVector n(start-stop);
	for(register int i=start, j=0;i>stop;i--,j++)
		n.vect[j]=vect[i];
	return n;
}


DVector DVector::prepad(int count) {
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);

	if (count<=0) //return DVector();
		throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);

	DVector n(count+taille);
	for(register int i=count, j=0; j<taille;i++,j++)
		n.vect[i]=vect[j];
	return n;
}

DVector DVector::postpad(int count) {
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (count<=0) //return DVector();
		throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	DVector n(count+taille);
	for(register int j=0; j<taille;j++)
		n.vect[j]=vect[j];
	return n;
}

void DVector::insert(const DVector &rhs, int pos) {
	if(rhs.taille<=0 || rhs.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(pos<0) 
		throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);

	for(register int i=0, j=pos; i < rhs.taille && j<taille;i++, j++)
		vect[j]=rhs.vect[i];
}

void DVector::reset(double value){
	if(vect==NULL ||taille==0) return;
	for(register int i=0;i<taille;i++) vect[i]=value;
}

ZVector::ZVector(const double *real, const double* imag, const int &size){
	if (size<0) throw CUtilitisException(__FILE__, __LINE__, INVALID_VECTOR_DIMENSION);
    if(size==0) { vect=NULL; taille=0; return; }
	if (real==NULL || imag==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	taille=size;
	vect=new DCplx[taille];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i=0;i<size;i++) {
		vect[i].re=real[i];
		vect[i].im=imag[i];
	}
}

ZVector::ZVector(const int &size){
	if (size<0) throw CUtilitisException(__FILE__, __LINE__, INVALID_VECTOR_DIMENSION);
    if(size==0) { vect=NULL; taille=0; return; }
	vect=new cmplx[size];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	taille=size;
}

ZVector::ZVector(const int &size, const DCplx &initial) {
	if (size<0) throw CUtilitisException(__FILE__, __LINE__, INVALID_VECTOR_DIMENSION);
    if(size==0) { vect=NULL; taille=0; return; }
	vect=new cmplx[size];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	taille=size;
	for(register int i=0;i<size;i++) vect[i]=initial;
}

ZVector::ZVector(const cmplx *v, const int size) {
	if (size<0) throw CUtilitisException(__FILE__, __LINE__, INVALID_VECTOR_DIMENSION);
    if(size==0) { vect=NULL; taille=0; return; }
    if(v==NULL) throw CUtilitisException(__FILE__,__LINE__, EMPTY_VECTOR);
	vect=new cmplx[size];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	taille=size;
	for(register int i=0;i<size;i++) vect[i]=v[i];
}

ZVector::ZVector(const ZVector &c){
	if(c.taille<0) throw CUtilitisException(__FILE__, __LINE__, INVALID_VECTOR_DIMENSION);
    if(c.taille==0 || c.vect==NULL) {
		taille=0; vect=NULL;
		return;
	}
    if(c.vect==NULL) throw CUtilitisException(__FILE__,__LINE__, EMPTY_VECTOR);
	taille=c.taille;
	vect=new cmplx[taille];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
}

ZVector::ZVector(const DVector &c){
	if(c.taille<0) throw CUtilitisException(__FILE__, __LINE__, INVALID_VECTOR_DIMENSION);
    if(c.taille==0 || c.vect==NULL) {
		taille=0; vect=NULL;
		return;
	}
    if(c.vect==NULL) throw CUtilitisException(__FILE__,__LINE__, EMPTY_VECTOR);
	taille=c.taille;
	vect=new cmplx[taille];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i=0;i<taille;i++) vect[i].re=c.vect[i];
}

ZVector::ZVector(const WVector &c){
	if(c.taille<0) throw CUtilitisException(__FILE__, __LINE__, INVALID_VECTOR_DIMENSION);
    if(c.taille==0 || c.vect==NULL) {
		taille=0; vect=NULL;
		return;
	}
    if(c.vect==NULL) throw CUtilitisException(__FILE__,__LINE__, EMPTY_VECTOR);

	taille=c.taille;
	vect=new cmplx[taille];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
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
	if(taille==0 || vect==NULL) //return DVector();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i].re;
	}
	return n;
}

DVector ZVector::imag()  const { 
	if(taille==0 || vect==NULL) //return DVector();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i].im;
	}
	return n;
}

//double       modul2(cmplx a);
DVector ZVector::mag() const  { 
	if(taille==0 || vect==NULL) //return DVector();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i].mag();
	}
	return n;
}

DVector ZVector::mod2() const  { 
	if(taille==0 || vect==NULL) // return DVector();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i].mod2();
	}
			
	return n;
}

DVector ZVector::phase() const  {
	if(taille==0 || vect==NULL) //return DVector();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	DVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i].phase();
	}
	return n;
}

void ZVector::fixedconv(const ZVector &rhs) const{
	if(rhs.taille<=0 || rhs.vect==NULL || taille<=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);

	register DCplx sum;
	int cache;

	/* Loop over output points */
	for (register int m = taille-1; m >= 0 ; m--) {
		/* Convolution */
		sum = 0.0;
		for (register int j = 0; j < rhs.taille; ++j) {
			cache=m-j;
			if (cache>=0 && cache<taille) sum += rhs.vect[j]* vect[cache];
		}
		vect[m] = sum;
	}
}

void ZVector::reset(const DCplx& c ){
	if(taille<=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	for(register int i=0;i<taille;i++) {
		vect[i]=c;
	}
}

void ZVector::conv(const ZVector &rhs) {
	if(rhs.taille<=0 || rhs.vect==NULL || taille<=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	register DCplx sum;
	register int cache;
	int imax=taille+rhs.taille-1;
	DCplx *tmp=new DCplx[imax];
	if(tmp==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);

	/* Loop over output points */
	for (register int i = 0; i<imax ; i++) {
		/* Convolution */
		sum = 0.0;
		for (register int j = 0; j < rhs.taille; ++j) {
			cache=i-j;
			if (cache>=0 && cache<taille) 
				sum += rhs.vect[j] * vect[cache];
		}
		tmp[i] = sum;
	}
	delete[] vect;
	taille=imax;
	vect=tmp;
}

ZVector ZVector::corr_nip(const ZVector &b){
	if(taille<=0 || vect==NULL || b.vect==NULL || b.taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	DCplx sum;
	ZVector n(2*nr_max(taille,b.taille)-1);;
	if(taille!=b.taille) cerr << "ZVector::corr : Warning. Vector size are not the same. output size will be 2*max(taille, b.taille)-1" << endl;
		
	for (register int i=b.taille-1 ; i>-taille ; i--){
		sum=0.0;//init_cmplx(0,0);
		for (register int j=0;j<taille; j++){
			register int cache=i+j;
			if (cache>=0 && cache<b.taille)
				sum+=vect[j]*b.vect[cache].conj();
		}
		n.vect[b.taille-1-i]=sum;
	}
	return n;
}

void ZVector::corr(const ZVector &b){
	if(taille<=0 || vect==NULL || b.vect==NULL || b.taille<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	DCplx sum;
	DCplx *n= new DCplx[2*nr_max(taille,b.taille)-1];
	if(taille!=b.taille) cerr << "ZVector::corr : Warning. Vector size are not the same. output size will be 2*max(taille, b.taille)-1" << endl;
		
	for (register int i=b.taille-1 ; i>-taille ; i--){
		sum=0.0;//init_cmplx(0,0);
		for (register int j=0;j<taille; j++) {
			register int cache=i+j;
			if (cache>=0 && cache<b.taille)
				sum+=vect[j]*b.vect[cache].conj();
		}
		n[b.taille-1-i]=sum;
	}
	delete vect;
	vect=n;
	taille=taille+b.taille-1;
}

void ZVector::xcorr() {
	if(taille<=0 || vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	register DCplx sum;
	DCplx *n= new DCplx[taille+taille-1];
		
	for (register int i=taille-1 ; i>-taille ; i--){
		sum=0.0;//init_cmplx(0,0);
		for (register int j=0;j<taille; j++) {
			register int cache=i+j;
			if (cache>=0 && cache<taille)
				sum+=vect[j]*vect[cache].conj();
		}
		n[taille-1-i]=sum;
	}
	delete vect;
	vect=n;
	taille=taille+taille-1;
}

ZVector ZVector::conv_nip(const ZVector &rhs) const{
	if(rhs.taille<=0 || rhs.vect==NULL || taille<=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);

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
			if (cache>=0 && cache<taille) 
				sum += rhs.vect[j] * vect[cache];
		}
		n.vect[i] = sum;
	}
	return n;
}

ZVector ZVector::xcorr_nip() const{
	if(taille<=0 || vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);

	DCplx sum;
	ZVector n(taille+taille-1);;
		
	for (register int i=taille-1 ; i>-taille ; i--){
		sum=0.0;//init_cmplx(0,0);
		for (register int j=0;j<taille; j++){
			register int cache=i+j;
			if (cache>=0 && cache<taille)
				sum+=vect[j]*vect[cache].conj();
		}
		n.vect[taille-1-i]=sum;
	}
	return n;
}

void ZVector::fixedconv(const DVector &rhs) {
	if(rhs.taille<=0 || rhs.vect==NULL || taille<=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
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
	if(rhs.taille<=0 || rhs.vect==NULL || taille<=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);

	register DCplx sum;
	register int cache;
	int imax=taille+rhs.taille-1;
	DCplx *tmp=new DCplx[imax];
	if(tmp==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);

	/* Loop over output points */
	for (register int i = 0; i<imax ; i++) {
		/* Convolution */
		sum = 0.0;
		for (register int j = 0; j < rhs.taille; ++j) {
			cache=i-j;
			if (cache>=0 && cache<taille) 
				sum += rhs.vect[j] * vect[cache];
		}
		tmp[i] = sum;
	}
	delete[] vect;
	taille=imax;
	vect=tmp;
}

ZVector ZVector::conv_nip(const DVector &rhs) const{
	if(rhs.taille<=0 || rhs.vect==NULL || taille<=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
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
			if (cache>=0 && cache<taille) 
				sum += rhs.vect[j] * vect[cache];
		}
		n.vect[i] = sum;
	}
	return n;
}

ZVector ZVector::copy(int start, int stop) const{
	if(taille<=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(start>stop) swap(start, stop);
	if (stop>taille) stop=taille; 
	if (start<0) start=0;
	if (!(stop-start)) return ZVector();
	ZVector n(stop-start);
	for(register int i=start, j=0;i<stop;i++,j++)
		n.vect[j]=vect[i];
	return n;
}

ZVector ZVector::rev_copy(int start, int stop) const{
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(stop>start) swap(stop, start);
	if (start>taille) start=taille; 
	if (stop<0) stop=0;
	if (stop==start) return ZVector();
	ZVector n(start-stop);
	for(register int i=start, j=0;i>stop;i--,j++)
		n.vect[j]=vect[i];
	return n;
}

ZVector ZVector::prepad(int count) {
	if (count<=0) //return ZVector();
		throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	if(taille<=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);

	ZVector n(count+taille);
	for(register int i=count, j=0; j<taille;i++,j++)
		n.vect[i]=vect[j];
	return n;
}

DCplx ZVector::sum() const{
	if(taille<=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	DCplx sum;
	//ZVector n(count+taille);
	for(register int i=0; i<taille;i++)
		sum+=vect[i];
	return sum;
}


double DVector::sum() const{
	if(taille<=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	register double sum=0.0;
	//ZVector n(count+taille);
	for(register int i=0; i<taille;i++)
		sum+=vect[i];
	return sum;
}

void DVector::normalizepower() {
	if(taille<=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	register double sum=0.0;
	//ZVector n(count+taille);
        register int i;
	for(i=0; i<taille;i++)
		sum+=vect[i]*vect[i];
	double norm = 1.0/sqrt(sum);
	for(i=0; i<taille;i++)
		vect[i]*=norm;
	
}


ZVector ZVector::postpad(int count) {
	if (count<=0) //return ZVector();
		throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	if(taille<=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);

	ZVector n(count+taille);
	for(register int j=0; j<taille;j++)
		n.vect[j]=vect[j];
	return n;
}

void ZVector::insert(const ZVector &rhs, int pos) {
	if(rhs.taille<=0 || rhs.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(pos<0) 
		throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	for(register int i=0, j=pos; i<rhs.taille && j<taille;i++, j++)
		vect[j]=rhs.vect[i];
}

void ZVector::insert(const DVector &rhs, int pos) {
	if(rhs.taille<=0 || rhs.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(pos<0) 
		throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);

	for(register int i=0, j=pos; i<pos && j<taille;i++, j++)
		vect[j]=rhs.vect[i];
}

//cmplx        cmplx_add(cmplx a,cmplx b);
ZVector ZVector::operator+(const ZVector &c)  const { 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c.vect[i];
	}
	return n;
}

ZVector ZVector::operator+(const DCplx &c)  const { 
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c;
	}
	return n;
}

ZVector ZVector::operator+(const double c)  const { 
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);

	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c;
	}
	return n;
}

ZVector ZVector::operator+(const DVector &c)  const { 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+c.vect[i];
	}
	return n;
}

ZVector ZVector::operator+(const WVector &c)  const { 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]+(double)c.vect[i];
	}
	return n;
}

ZVector ZVector::operator-() { 
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=-vect[i];//-c.vect[i];
	}
	return n;
}

// cmplx        cmplx_sub(cmplx a,cmplx b);
ZVector ZVector::operator-(const ZVector &c)  const { 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c.vect[i];
	}
	return n;
}

ZVector ZVector::operator-(const DCplx &c)  const { 
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c;
	}
	return n;
}

ZVector ZVector::operator-(const double c)  const { 
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c;
	}
	return n;
}

ZVector ZVector::operator-(const DVector &c) const  { 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-c.vect[i];
	}
	return n;
}

ZVector ZVector::operator-(const WVector &c)  const { 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]-(double)c.vect[i];
	}
	return n;
}


//cmplx        cmplx_mul(cmplx a,cmplx b);
ZVector ZVector::operator*(const ZVector &c)  const { 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c.vect[i];
	}
	return n;
}

ZVector ZVector::operator*(const DCplx &c)  const { 
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c;
	}
	return n;
}

ZVector ZVector::operator*(const double c)  const { 
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);

	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c;
	}
	return n;
}

ZVector ZVector::operator*(const DVector &c)  const { 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*c.vect[i];
	}
	return n;
}

ZVector ZVector::operator*(const WVector &c)  const { 
	if(c.taille<=0 || c.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	ZVector n(taille);
	for(register int i=0;i<taille;i++){
		n.vect[i]=vect[i]*(double)c.vect[i];
	}
	return n;
}

ZVector ZVector::operator/(const double c)  const { 
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (c==0.0) 
		throw CUtilitisException(__FILE__, __LINE__, DIVISION_BY_ZERO);
	double d=1/c;
	return ZVector((*this)*d);
}

ZVector ZVector::operator/(const DCplx &c)  const {
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	double t=c.mod2();
	if (t==0.0) 
		throw CUtilitisException(__FILE__, __LINE__, DIVISION_BY_ZERO);

	DCplx ct=c.conj()/t;
	return ZVector((*this)*ct);
}

ZVector &ZVector::operator=(const ZVector &c) { 

    if(taille==c.taille && vect!=NULL) {
        for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
        return *this;
    }
    
	if(taille!=0 && vect!=NULL) delete[] vect;
	taille=c.taille;
	vect = new cmplx[taille];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i=0;i<taille;i++) vect[i]=c.vect[i];
	return *this;
}

ZVector &ZVector::operator=(const DVector &c){
    if(taille==c.taille && vect!=NULL) {
        for(register int i=0;i<taille;i++) {
            vect[i].re=c.vect[i];
        }
        return *this;
    }
	if(taille!=0 && vect!=NULL) delete[] vect;
	taille=c.taille;
	vect = new cmplx[taille];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i=0;i<taille;i++) {
		vect[i].re=c.vect[i];
	}
	return *this;
}

ZVector &ZVector::operator=(const WVector &c){
    if(taille==c.taille && vect!=NULL) {
        for(register int i=0;i<taille;i++) {
            vect[i].re=c.vect[i];
        }
        return *this;
    }
	if(taille!=0 && vect!=NULL) delete[] vect;
	taille=c.taille;
	vect = new cmplx[taille];
	if(vect==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
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

ZVector ZVector::operator<<(int n) const { // n = 1..taille
	if(taille<=0 || vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(n>taille || n<0) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	if (n==0) return *this;
	ZVector nz(taille, DCplx(0,0));
	for (int i=0, j=n; j<taille; j++,i++) {
		nz.vect[i]=vect[j];
	}
	return nz;
}

ZVector ZVector::operator>>(int n) const { // n = 1..taille
	if(taille<=0 || vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(n>taille || n<0) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	if (n==0) return *this;
	ZVector nz(taille, DCplx(0,0));
	for (int i=0, j=n; j<taille; j++,i++) {
		nz.vect[j]=vect[i];
	}
	return nz;
}

WVector WVector::operator<<(int n) const { // n = 1..taille
	if(taille<=0 || vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(n>taille || n<0) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	if (n==0) return *this;
	WVector nz(taille, 0);
	for (int i=0, j=n; j<taille; j++,i++) {
		nz.vect[i]=vect[j];
	}
	return nz;
}

WVector WVector::operator>>(int n) const { // n = 0..taille-1
	if(taille<=0 || vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(n>taille || n<0) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	if (n==0) return *this;
	WVector nz(taille, 0);
	for (int i=0, j=n; j<taille; j++,i++) {
		nz.vect[j]=vect[i];
	}
	return nz;
}

DVector DVector::operator<<(int n) const { // n = 1..taille
	if(taille<=0 || vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(n>taille || n<0) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	if (n==0) return *this;
	DVector nz(taille, 0);
	for (int i=0, j=n; j<taille; j++,i++) {
		nz.vect[i]=vect[j];
	}
	return nz;
}

DVector DVector::operator>>(int n) const { // n = 1..taille
	if(taille<=0 || vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(n>taille || n<0) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	if (n==0) return *this;
	DVector nz(taille, 0);
	for (int i=0, j=n; j<taille; j++,i++) {
		nz.vect[j]=vect[i];
	}
	return nz;
}


ZVector& ZVector::operator<<=(int n)  { // n = 1..taille
	if(taille<=0 || vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(n>taille || n<0) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	if (n==0) return *this;
	ZVector nz(taille, DCplx(0,0));
	for (int i=0, j=n; j<taille; j++,i++) {
		nz.vect[i]=vect[j];
	}
	*this = nz;
	return *this;
}

ZVector& ZVector::operator>>=(int n)  { // n = 1..taille
	if(taille<=0 || vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(n>taille || n<0) throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	if (n==0) return *this;
	ZVector nz(taille, DCplx(0,0));
	for (int i=0, j=n-1; j<taille; j++,i++) {
		nz.vect[j]=vect[i];
	}
	*this = nz;
	return *this;
}

DCplx ZVector::DotProd(const DVector &rhs){
	if(rhs.taille<=0 || rhs.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (rhs.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	DCplx sum(0.0,0.0);
	for(register int i=0;i<taille;i++) {
		sum+=vect[i]*rhs.vect[i];
	}
	return sum;
}

DCplx ZVector::DotProd(const WVector &rhs){
	if(rhs.taille<=0 || rhs.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (rhs.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	DCplx sum(0.0,0.0);
	for(register int i=0;i<taille;i++) {
		sum+=vect[i]*(double)rhs.vect[i];
	}
	return sum;
}

DCplx ZVector::DotProd(const ZVector &rhs){
	if(rhs.taille<=0 || rhs.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (rhs.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	DCplx sum(0.0,0.0);
	for(register int i=0;i<taille;i++) {
		sum+=rhs.vect[i]*vect[i];
	}
	return sum;
}


ZVector& ZVector::operator +=(const ZVector &zv){
	if(zv.taille<=0 || zv.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (zv.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	for (register int i=0; i < taille; i++) {
		vect[i]+=zv.vect[i];
	}
	return (*this);
}

ZVector& ZVector::operator +=(const DVector &dv){
	if(dv.taille<=0 || dv.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (dv.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	for (register int i=0; i < taille; i++) {
		vect[i]+=dv.vect[i];
	}
	return (*this);
}

ZVector& ZVector::operator +=(const WVector &wv){
	if(wv.taille<=0 || wv.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (wv.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	for (register int i=0; i < taille; i++) {
		vect[i]+=wv.vect[i];
	}
	return (*this);
}

ZVector& ZVector::operator +=(const DCplx &c){
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);

	for (register int i=0; i < taille; i++) {
		vect[i]+=c;
	}
	return (*this);
}

ZVector& ZVector::operator +=(double c){
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);

	for (register int i=0; i < taille; i++) {
		vect[i]+=c;
	}
	return (*this);
}

ZVector& ZVector::operator -=(const ZVector &zv){
	if(zv.taille<=0 || zv.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (zv.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	
	for (register int i=0; i < taille; i++) {
		vect[i]-=zv.vect[i];
	}
	return (*this);
}

ZVector& ZVector::operator -=(const DVector &dv){
	if(dv.taille<=0 || dv.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (dv.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	for (register int i=0; i < taille; i++) {
		vect[i]-=dv.vect[i];
	}
	return (*this);
}

ZVector& ZVector::operator -=(const WVector &wv){
	if(wv.taille<=0 || wv.vect==NULL || taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (wv.taille != taille) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	for (register int i=0; i < taille; i++) {
		vect[i]-=wv.vect[i];
	}
	return (*this);
}

ZVector& ZVector::operator -=(const DCplx &c){
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);

	for (register int i=0; i < taille; i++) {
		vect[i]-=c;
	}
	return (*this);
}

ZVector& ZVector::operator -=(double c){
	if(taille <=0 || vect==NULL) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);

	for (register int i=0; i < taille; i++) {
		vect[i]-=c;
	}
	return (*this);
}

DMatrix::DMatrix(int line, int col, double initial){
	if (col==0 || line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=col;this->line=line;
	int linecol= line*col;
	mat=new double*[line];
	if(mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	mat[0]=new double[linecol];
	if(mat[0]==NULL) {
		delete[] mat;
		throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	}
	for(register int j=0,i=0;j<linecol;j++) {
		mat[0][j]=initial;
		if (j%col==0) { 
			mat[i]=&mat[0][i*col]; 
			i++;
		}
	}
}

DMatrix::DMatrix(const DMatrix& m) {
	if (m.col==0 || m.line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=m.col;this->line=m.line;

	int linecol= line*col;
	mat=new double*[line];
	if(mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	mat[0]=new double[linecol];
	if(mat[0]==NULL) {
		delete[] mat;
		throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	}
	for(register int i=0;i<line;i++) 
		mat[i]=&mat[0][i*col];

	memcpy(mat[0], m.mat[0], sizeof(double)*linecol);
}

DMatrix::DMatrix(const WMatrix& m) {
	if (m.col==0 || m.line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=m.col;this->line=m.line;
	int linecol= line*col;
	mat=new double*[line];
	if(mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	mat[0]=new double[linecol];
	if(mat[0]==NULL) {
		delete[] mat;
		throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	}
	int jp=0;
	for(register int j=0,i=0;j<linecol;j++, jp=j%col) {
		mat[0][j]=m.mat[i][jp];
		if (j!=0&&jp==0) { 
			i++;
			mat[i]=&mat[0][i*col]; 
		}
	}
}

DMatrix::DMatrix() { col=line=0;mat=NULL; }

DMatrix::~DMatrix() {

	if(col==0 || mat==0) return;
	if(line==0) {delete[] mat; return;}

	delete[] mat[0];
	delete[] mat;
	line=col=0;
	mat=NULL;
}

void DMatrix::transpose() {
	if (mat==NULL ||col==0 ||line==0) return;
	DMatrix t(col, line); 
	for(register int i=0;i<col;i++) {
		for(register int j=0;j<line;j++) {
			t.mat[i][j]=mat[j][i];
		}
	}
	double **tmp=mat;
	mat=t.mat;
	t.mat=tmp;
}

DMatrix DMatrix::transpose_nip(){
	if (mat==NULL ||col==0 ||line==0) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix t(col, line); 
	for(register int i=0;i<col;i++) {
		for(register int j=0;j<line;j++) {
			t.mat[i][j]=mat[j][i];
		}
	}
	return t;
}

ZMatrix DMatrix::ztranspose_nip(){
	if (mat==NULL ||col==0 ||line==0) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);

	ZMatrix t(col, line); // ; t.col=line;t.line=col;

	for(register int i=0; i<col;i++) {
		for(register int j=0;j<line;j++)
			t.mat[i][j].re=mat[j][i];
	}
	return t;
}

void DMatrix::eye(){
	if(mat==NULL || line<=0 || col<=0)  // return;
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(line!=col)
		throw CUtilitisException(__FILE__, __LINE__, MATRIX_NOT_SQUARE);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			mat[i][j] = (i!=j) ? 0 : 1;
		}
	}
}

void DMatrix::antieye(){
	if(mat==NULL || line<=0 || col<=0)  // return;
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(line!=col)
		throw CUtilitisException(__FILE__, __LINE__, MATRIX_NOT_SQUARE);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j--) {
			mat[i][j] = (i!=(col-1-j)) ? 0 : 1;
		}
	}
}


DMatrix DMatrix::operator+(const DMatrix &rhs){
	if (col==0 || line==0 || mat==NULL)
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.col!=col || rhs.line!=line) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	DMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]+rhs.mat[i][j];
	}
	return n;
}

DMatrix DMatrix::operator-(){
if (col==0 || line==0 || mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=-mat[i][j];
	}
	return n;
}

DMatrix DMatrix::operator-(const DMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) //return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.col!=col || rhs.line!=line) //return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	DMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]-rhs.mat[i][j];
	}
	return n;
}

DMatrix DMatrix::operator*(const DMatrix &rhs){
	if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) //return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.line!=col) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);

	DMatrix n(line, rhs.col);

	for(register int j=0;j<n.line;j++) { //this.line
		for(register int i=0;i<n.col;i++) { // rhs.col
			register double sum=0.0;
			for(register int k=0;k<col;k++)  // this.col
				sum+=mat[j][k]*rhs.mat[k][i];
			n.mat[j][i]=sum;
		}
	}
	return n;
}

ZMatrix DMatrix::operator+(const ZMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.col!=col || rhs.line!=line) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(line,  col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]+rhs.mat[i][j];
	}
	return n;
}

ZMatrix DMatrix::operator-(const ZMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.col!=col || rhs.line!=line) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(line,  col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]-rhs.mat[i][j];
	}
	return n;
}

ZMatrix DMatrix::operator*(const ZMatrix &rhs){
	if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.line!=col) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(rhs.col, line);

	for(register int j=0;j<n.line;j++) {
		for(register int i=0;i<n.col;i++) {
			register DCplx sum(0.0,0.0);
			for(register int k=0;k<col;k++)
				sum+=mat[j][k]*rhs.mat[k][i];
			n.mat[j][i]=sum;
		}
	}
	return n;
}

DMatrix DMatrix::operator*(const double rhs){
	if (col==0 || line==0 || mat==NULL) //return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(line,  col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]*rhs;
	}
	return n;
}

ZMatrix DMatrix::operator*(const DCplx &rhs){
	if (col==0 || line==0 || mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	ZMatrix n(line,  col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]*rhs;
	}
	return n;
}

DMatrix DMatrix::operator/(const double rhs){
	if(rhs==0.0) //return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, DIVISION_BY_ZERO);
	double irhs=1/rhs;
	return (*this)*irhs;
}

ZMatrix DMatrix::operator/(const DCplx &rhs){
	if(rhs==0.0) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, DIVISION_BY_ZERO);
	DCplx irhs=1.0/rhs;
	return (*this)*irhs;
}

DMatrix &DMatrix::operator=(const DMatrix &rhs){
	if(mat!=NULL&&line>0 && col>0) {
		//for(register int i=0;i<line;i++) delete[] mat[i];
		delete[] mat[0];
		delete[] mat;
	}

	line=rhs.line;col=rhs.col;

	if(line==0|| col==0 || rhs.mat==NULL) {
		this->mat=NULL; return *this;
	}

	int linecol= line*col;
	mat=new double*[line];
	if(mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	mat[0]=new double[linecol];
	if(mat[0]==NULL) {
		delete[] mat;
		throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	}

	//int jp=0;
	for(register int i=0;i<line;i++)
		mat[i] = &mat[0][i*col];
	memcpy(mat[0], rhs.mat[0], sizeof(double)*linecol);

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

void DMatrix::reset(double val) {
	if(col<=0 || line<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			mat[i][j]=val;
		}
	}
}

void DMatrix::inv() {
	if(mat==NULL || line <=0 || col<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(line!=col) throw CMatrixInversionException(INVALID_MATRIX_DIMENSION);;
	DMatrix dummy(line, col);
	gaussj(*this, dummy);
}

DMatrix DMatrix::inv_nip() {
	if(mat==NULL || line <=0 || col<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(line!=col) throw CMatrixInversionException(INVALID_MATRIX_DIMENSION);;
	DMatrix inv(*this);
	DMatrix dummy(line, col);
	gaussj(inv, dummy);
	return inv;
}

double DMatrix::sum() {
	if (col==0 || line==0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	double  sum=0.0;
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			sum+=mat[i][j];
		}
	}
	return sum;
}

DMatrix DMatrix::linesum(){
	if (col==0 || line==0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix sum(line,1,0.0);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			sum.mat[i][1]+=mat[i][j];
		}
	}
	return sum;
}

DMatrix DMatrix::colsum(){
	if (col==0 || line==0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix sum(1,col,0.0);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			sum.mat[0][j]+=mat[i][j];
		}
	}
	return sum;
}


ZMatrix::ZMatrix (int line, int col){ // converted
	if (col==0 || line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=col;this->line=line;
	mat=new DCplx*[line];
	if(mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	DCplx initial;
	mat[0] = new DCplx[line*col];
	if (mat[0]==NULL) { delete[] mat; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }

	for(register int j=0;j<line*col;j++) mat[0][j]=initial;

	for(register int i=1;i<line;i++) {
		mat[i]=&mat[0][i*col];
	}
}

ZMatrix::ZMatrix (int line, int col, const DCplx &initial)  { // converted
	if (col==0 || line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=col;this->line=line;
	mat=new DCplx*[line];
	if(mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	mat[0] = new DCplx[line*col];
	if (mat[0]==NULL) { delete[] mat; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }

	for(register int j=0;j<line*col;j++) mat[0][j]=initial;

	for(register int i=1;i<line;i++) {
		mat[i]=&(mat[0][i*col]);
	}
}

ZMatrix::ZMatrix(const ZMatrix& m){ // converted
	if (m.col==0 || m.line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=m.col;this->line=m.line;
	mat=new DCplx*[line];
	if(mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	mat[0] = new DCplx[line*col];
	if (mat[0]==NULL) { delete[] mat; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }

	memcpy(mat[0], m.mat[0], sizeof(DCplx)*line*col);

	for(register int i=1;i<line;i++) {
		mat[i]=&(mat[0][i*col]);
	}
}

ZMatrix::ZMatrix(const DMatrix& m){ // converted
	if (m.col==0 || m.line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=m.col;this->line=m.line;
	mat=new DCplx*[line];
	if(mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	mat[0] = new DCplx[line*col];
	if (mat[0]==NULL) { delete[] mat; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }

	for(register int i=0;i<line;i++) {
		if(i) mat[i]=&(mat[0][i*col]);
		for(register int j=0;j<col;j++) 
			mat[i][j]=m.mat[i][j];
	}
}

ZMatrix::ZMatrix() { col=line=0;mat=NULL;}  // converted

ZMatrix::~ZMatrix() { // converted
	if(col==0 || mat==0) return;
	if(line==0) {delete[] mat; return;}
	//for(register int i =0; i<line;i++) delete[] mat[i];
	delete[] mat[0];
	delete[] mat;
	line=col=0;
	mat=NULL;
}
//void Free();
void ZMatrix::transpose(){ // converted
	if (mat==NULL ||col==0 ||line==0) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);

	DCplx **tmp = new DCplx *[col];
	if(tmp==NULL) { throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }

	tmp[0] = new DCplx[line*col];
	if (tmp[0]==NULL) { delete[] tmp; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }

	register int i=0;

	for(i=1;i<col;i++) tmp[i]=&(tmp[0][i*col]);

	for(i=0; i<col;i++) {
		for(register int j=0;j<line;j++)
			tmp[i][j]=mat[j][i];
	}

	delete[] mat[0];
	delete[] mat;

	int x=col;col=line;line=x;
	mat=tmp;
}


ZMatrix  ZMatrix::transpose_nip() const { // converted
	if (mat==NULL ||col==0 ||line==0) 
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	ZMatrix t(col, line); 
	for(register int i=0; i<col;i++) {
		for(register int j=0;j<line;j++)
			t.mat[i][j]=mat[j][i];
	}
	return t;
}

ZMatrix  ZMatrix::operator+(const ZMatrix &rhs) const {
	if (col==0 || line==0 || mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.col!=col || rhs.line!=line) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]+rhs.mat[i][j];
	}
	return n;
}

ZMatrix  ZMatrix::operator+(const DMatrix &rhs) const {
	if (col==0 || line==0 || mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.col!=col || rhs.line!=line) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]+rhs.mat[i][j];
	}
	return n;
}

ZMatrix ZMatrix::operator-(){
if (col==0 || line==0 || mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=-mat[i][j];
	}
	return n;
}

ZMatrix  ZMatrix::operator-(const ZMatrix &rhs) const {
	if (col==0 || line==0 || mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.col!=col || rhs.line!=line) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]-rhs.mat[i][j];
	}
	return n;
}

ZMatrix  ZMatrix::operator-(const DMatrix &rhs) const {
	if (col==0 || line==0 || mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.col!=col || rhs.line!=line) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]-rhs.mat[i][j];
	}
	return n;
}

ZMatrix ZMatrix::operator*(const ZMatrix &rhs) const {
	if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.line!=col) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);

	ZMatrix n(line , rhs.col);

	for(register int j=0;j<n.line;j++) {
		for(register int i=0;i<n.col;i++) {
			register DCplx sum(0.0,0.0);
			for(register int k=0;k<col;k++)
				sum+=mat[j][k]*rhs.mat[k][i];
			n.mat[j][i]=sum;
		}
	}
	return n;
}

ZMatrix ZMatrix::operator*(const DMatrix &rhs) const {
	if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.line!=col) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);

	ZMatrix n(line, rhs.col);

	for(register int j=0;j<n.line;j++) {
		for(register int i=0;i<n.col;i++) {
			register DCplx sum(0.0,0.0);
			for(register int k=0;k<col;k++)
				sum+=mat[j][k]*rhs.mat[k][i];
			n.mat[j][i]=sum;
		}
	}
	return n;
}

ZMatrix ZMatrix::operator*(const double rhs) const {
	if (col==0 || line==0 || mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]*rhs;
	}
	return n;
}

ZMatrix ZMatrix::operator*(const DCplx &rhs) const {
	if (col==0 || line==0 || mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]*rhs;
	}
	return n;
}

ZMatrix ZMatrix::operator/(const double rhs) const {
	if(rhs==0.0) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, DIVISION_BY_ZERO);
	double irhs=1.0/rhs;
	return (*this)*irhs;
}

ZMatrix ZMatrix::operator/(const DCplx &rhs) const {
	if(rhs==0.0) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, DIVISION_BY_ZERO);
	DCplx irhs=1.0/rhs;
	return (*this)*irhs;
}

ZMatrix &ZMatrix::operator=(const ZMatrix &rhs) { // converted
	if(mat!=NULL&&line>0&&col>0) {
		delete[] mat[0];
		delete[] mat;
	}

	line=rhs.line;col=rhs.col;

	if(line==0|| col==0 || rhs.mat==NULL) {
		this->mat=NULL; line=col=0; return *this;
	}

	mat=new DCplx*[line];
	if(mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	mat[0] = new DCplx[line*col];
	if (mat[0]==NULL) { delete[] mat; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }

	for(register int i = 0;i<line;i++) {
		mat[i]=&(mat[0][i*col]);
		for(register int j=0;j<col;j++)
			mat[i][j]=rhs.mat[i][j];
	}
	return *this;
}

ZMatrix &ZMatrix::operator=(const DMatrix &rhs){ // converted
	if(mat!=NULL&&line>0&&col>0) {
		delete[] mat[0];
		delete[] mat;
	}

	line=rhs.line;col=rhs.col;

	if(line==0|| col==0 || rhs.mat==NULL) {
		this->mat=NULL; line=col=0; return *this;
	}

	mat=new DCplx*[line];
	if(mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	mat[0] = new DCplx[line*col];
	if (mat[0]==NULL) { delete[] mat; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }

	for(register int i = 0;i<line;i++) {
		mat[i]=&(mat[0][i*col]);
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

ZMatrix &ZMatrix::operator+=(const ZMatrix &rhs){
	if(rhs.col<=0 || rhs.line<=0 || rhs.mat==NULL || col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(rhs.col!=col || rhs.line!=line) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			mat[i][j]+=rhs.mat[i][j];
		}
	}
	return *this;
}

ZMatrix &ZMatrix::operator+=(const DMatrix &rhs){
	if(rhs.col<=0 || rhs.line<=0 || rhs.mat==NULL || col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(rhs.col!=col || rhs.line!=line) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			mat[i][j]+=rhs.mat[i][j];
		}
	}
	return *this;
}

ZMatrix &ZMatrix::operator+=(const WMatrix &rhs){
	if(rhs.col<=0 || rhs.line<=0 || rhs.mat==NULL || col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(rhs.col!=col || rhs.line!=line) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			mat[i][j]+=rhs.mat[i][j];
		}
	}
	return *this;
}

ZMatrix &ZMatrix::operator+=(const DCplx &rhs){
	if(col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			mat[i][j]+=rhs;
		}
	}
	return *this;
}

ZMatrix &ZMatrix::operator+=(double r){
	if(col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			mat[i][j].re+=r;
		}
	}
	return *this;
}

//***
ZMatrix &ZMatrix::operator-=(const ZMatrix &rhs){
	if(rhs.col<=0 || rhs.line<=0 || rhs.mat==NULL || col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(rhs.col!=col || rhs.line!=line) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			mat[i][j]-=rhs.mat[i][j];
		}
	}
	return *this;
}

ZMatrix &ZMatrix::operator-=(const DMatrix &rhs){
	if(rhs.col<=0 || rhs.line<=0 || rhs.mat==NULL || col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(rhs.col!=col || rhs.line!=line) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			mat[i][j]-=rhs.mat[i][j];
		}
	}
	return *this;
}

ZMatrix &ZMatrix::operator-=(const WMatrix &rhs){
	if(rhs.col<=0 || rhs.line<=0 || rhs.mat==NULL || col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(rhs.col!=col || rhs.line!=line) 
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			mat[i][j]-=rhs.mat[i][j];
		}
	}
	return *this;
}

ZMatrix &ZMatrix::operator-=(const DCplx &rhs){
	if(col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			mat[i][j]-=rhs;
		}
	}
	return *this;
}

ZMatrix &ZMatrix::operator-=(double r){
	if(col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			mat[i][j].re-=r;
		}
	}
	return *this;
}

ZMatrix ZMatrix::operator-(const DCplx &c) const {
	if(col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	ZMatrix n(line,col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			n.mat[i][j]+=mat[i][j]-c;
		}
	}
	return n;
}

ZMatrix ZMatrix::operator-(double c) const {
	if(col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	ZMatrix n(line,col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			n.mat[i][j]+=mat[i][j]-c;
		}
	}
	return n;
}

ZMatrix ZMatrix::operator+(const DCplx &c) const {
	if(col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	ZMatrix n(line,col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			n.mat[i][j]+=mat[i][j]+c;
		}
	}
	return n;
}

ZMatrix ZMatrix::operator+(double c) const {
	if(col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	ZMatrix n(line,col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			n.mat[i][j]+=mat[i][j]+c;
		}
	}
	return n;
}

void ZMatrix::swapline(int l1, int l2) {
	if(col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(l1 < 0 || l2 < 0 || l1 >= line || l2 >=line) 
		throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	//swap(mat[l1], mat[l2]);
	for(register int i=0;i<col;i++) swap(mat[l1][i], mat[l2][i]);
}

void DMatrix::swapline(int l1, int l2) {
	if(col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(l1 < 0 || l2 < 0 || l1 >= line || l2 >=line) 
		throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	for(register int i=0;i<col;i++) swap(mat[l1][i], mat[l2][i]);
}

void WMatrix::swapline(int l1, int l2) {
	if(col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(l1 < 0 || l2 < 0 || l1 >= line || l2 >=line) 
		throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	for(register int i=0;i<col;i++) swap(mat[l1][i], mat[l2][i]);
}

void ZMatrix::swapcol(int c1, int c2) {
	if(col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(c1 < 0 || c2 < 0 || c1 >= col || c2 >=col ) 
		throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	for(register int i=0;i<line;i++) swap(mat[i][c1], mat[i][c2]);
}

void DMatrix::swapcol(int c1, int c2) {
	if(col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(c1 < 0 || c2 < 0 || c1 >= col || c2 >=col ) 
		throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	for(register int i=0;i<line;i++) swap(mat[i][c1], mat[i][c2]);
}

void WMatrix::swapcol(int c1, int c2) {
	if(col<=0 || line<=0 || mat==NULL)  
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(c1 < 0 || c2 < 0 || c1 >= col || c2 >=col ) 
		throw CUtilitisException(__FILE__, __LINE__, INDEX_OUT_OF_RANGE);
	for(register int i=0;i<line;i++) swap(mat[i][c1], mat[i][c2]);
}


DMatrix ZMatrix::real() const { // converted
	if(col==0 || line==0 || mat==NULL) // return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(line, col);
	//n.mat=new double*[line];
//	if(n.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i = 0;i<line;i++) {
//		n.mat[i]=new double[col];
//		if(n.mat[i]==NULL) { for(i=i-1;i>=0;i++) delete[] n.mat[i]; delete[] n.mat; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j].re;
	}
	return n;
}

DMatrix ZMatrix::imag() const { // converted
	if(col==0 || line==0 || mat==NULL) //return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(line, col);
//	n.mat=new double*[line];
//	if(n.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i = 0;i<line;i++) {
//		n.mat[i]=new double[col];
//		if(n.mat[i]==NULL) { for(i=i-1;i>=0;i++) delete[] n.mat[i]; delete[] n.mat; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j].im;
	}
	return n;
}

DMatrix ZMatrix::mag() const {
	if(col==0 || line==0 || mat==NULL) // return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(line, col);
//	n.mat=new double*[line];
//	if(n.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i = 0;i<line;i++) {
//		n.mat[i]=new double[col];
//		if(n.mat[i]==NULL) { for(i=i-1;i>=0;i++) delete[] n.mat[i]; delete[] n.mat; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j].mag();
	}
	return n;
}

DMatrix ZMatrix::phase() const {
	if(col==0 || line==0 || mat==NULL) //return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(line, col);
//	n.mat=new double*[line];
//	if(n.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i = 0;i<line;i++) {
//		n.mat[i]=new double[col];
//		if(n.mat[i]==NULL) { for(i=i-1;i>=0;i++) delete[] n.mat[i]; delete[] n.mat; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j].phase();
	}
	return n;
}


ZMatrix ZMatrix::hermit() const {
	if (mat==NULL ||col==0 ||line==0) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
//	ZMatrix t; t.col=line;t.line=col;
	ZMatrix t(col, line);

	//t.mat = new DCplx*[col];
	for(register int i=0; i<col;i++) {
	//	t.mat[i]=new DCplx[line];
		for(register int j=0;j<line;j++)
			t.mat[i][j]=mat[j][i].conj();
	}
	return t;
}

ZMatrix ZMatrix::mod2() const {
	if (mat==NULL ||col==0 ||line==0) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
//	ZMatrix t; t.col=line;t.line=col;
	if (line!=col) cerr << "Warning : Matrix not square for ZMatrix::mod2().. will apply z.hermit*z"<<endl;
	return this->hermit()*(*this);
}

ZMatrix ZMatrix::mod2r() const {
	if (mat==NULL ||col==0 ||line==0) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
//	ZMatrix t; t.col=line;t.line=col;
	if (line!=col) cerr << "Warning : Matrix not square for ZMatrix::mod2().. will apply z*z.hermit"<<endl;
	return (*this)*this->hermit();
}

ZMatrix ZMatrix::conj() const { 
	if(col==0 || line==0 || mat==NULL)// return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	ZMatrix n(line, col);
//	n.mat=new DCplx*[line];
//	if(n.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i = 0;i<line;i++) {
//		n.mat[i]=new DCplx[col];
//		if(n.mat[i]==NULL) { for(i=i-1;i>=0;i++) delete[] n.mat[i]; delete[] n.mat; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j].conj();
	}
	return n;
}

void ZMatrix::eye(){
	if(mat==NULL || line<=0 || col<=0)  // return;
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(line!=col)
		throw CUtilitisException(__FILE__, __LINE__, MATRIX_NOT_SQUARE);
	for(register int i=0;i<line;i++) {
		for(int j=0;j<col;j++) {
			mat[i][j] = (i!=j) ? 0 : 1;
		}
	}
}

void ZMatrix::antieye(){
	if(mat==NULL || line<=0 || col<=0)  // return;
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(line!=col)
		throw CUtilitisException(__FILE__, __LINE__, MATRIX_NOT_SQUARE);

	for(register int i=0;i<line;i++) {
		for(int j=0;j<col;j++) {
			mat[i][j] = (i!=(col-1-j)) ? 0 : 1;
		}
	}
}

void ZMatrix::reset(DCplx val) {
	if(col<=0 || line<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			mat[i][j]=val;
		}
	}
}


DMatrix operator*(const double lhs, const DMatrix &rhs){
	if (rhs.col==0 || rhs.line==0 || rhs.mat==NULL) //return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(rhs.line, rhs.col);
	for(register int i=0;i<rhs.line;i++) {
		for(register int j=0;j<rhs.col;j++)
			n.mat[i][j]=rhs.mat[i][j]*lhs;
	}
	return n;
}

ZMatrix operator*(const DCplx &lhs, const DMatrix &rhs){
	if (rhs.col==0 || rhs.line==0 || rhs.mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
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
	if(mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i =0;i<line;i++) {
		mat[i]=new int[col];
		if(mat[i]==NULL) { for(i=i-1;i>=0;i++) delete[] mat[i]; delete[] mat; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }

		for(register int j=0;j<col;j++) mat[i][j]=initial;
	}
}

WMatrix::WMatrix(const WMatrix& m) {
	if (m.col==0 || m.line==0) {
		this->col=this->line=0; mat=NULL; return;
	}
	this->col=m.col;this->line=m.line;
	mat=new int*[line];
	if(mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	size_t col_size=sizeof(int)*col;
	for(register int i =0;i<line;i++) {
		mat[i]=new int[col];
		if(mat[i]==NULL) { for(i=i-1;i>=0;i++) delete[] mat[i]; delete[] mat; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }
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
	if (mat==NULL ||col==0 ||line==0) //return;
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	
	int **tmp = new int*[col];
	if(mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	if(tmp==NULL) { throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }

	
	register int i=0;
	for(i=0; i<col;i++) {
		tmp[i]=new int[line];
		if(tmp[i]==NULL) { for(i=i-1;i>=0;i++) delete[] tmp[i];  delete[] tmp; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); } 
		for(register int j=0;j<line;j++)
			tmp[i][j]=mat[j][i];
	}
	for(i=0;i<line;i++) delete[] mat[i];
	delete[] mat;

	int x=col;col=line;line=x;
	mat=tmp;
}

WMatrix WMatrix::transpose_nip(){
	if (mat==NULL ||col==0 ||line==0) //return WMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);

	WMatrix t; t.col=line;t.line=col;

	t.mat = new int*[col];
	if(t.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i=0; i<col;i++) {
		t.mat[i]=new int[line];
		if(t.mat[i]==NULL) { for(i=i-1;i>=0;i++) delete[] t.mat[i]; delete[] t.mat; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }
		for(register int j=0;j<line;j++)
			t.mat[i][j]=mat[j][i];
	}
	return t;
}

ZMatrix WMatrix::ztranspose_nip(){
	if (mat==NULL ||col==0 ||line==0) // return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	ZMatrix t(col, line); t.col=line;t.line=col;

//	t.mat = new DCplx*[col];
//	if(t.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i=0; i<col;i++) {
//		t.mat[i]=new DCplx[line];
//		if(t.mat[i]==NULL) { for(i=i-1;i>=0;i++) delete[] t.mat[i]; delete[] t.mat; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }
		for(register int j=0;j<line;j++)
			t.mat[i][j].re=(double)mat[j][i];
	}
	return t;
}


void WMatrix::eye(){
	if(mat==NULL || line<=0 || col<=0 || line!=col) //return;
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);

	for(register int i=0;i<line;i++) {
		for(int j=0;j<col;j++) {
			mat[i][j] = (i!=j) ? 0 : 1;
		}
	}
}

WMatrix WMatrix::operator+(const WMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) //return WMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.col!=col || rhs.line!=line) //return WMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	WMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]+rhs.mat[i][j];
	}
	return n;
}

WMatrix WMatrix::operator-(){
if (col==0 || line==0 || mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	WMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=-mat[i][j];
	}
	return n;
}

WMatrix WMatrix::operator-(const WMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) //return WMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.col!=col || rhs.line!=line) //return WMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	WMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]-rhs.mat[i][j];
	}
	return n;
}

WMatrix WMatrix::operator*(const WMatrix &rhs){
	if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) //return WMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.line!=col) //return WMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	WMatrix n(line, rhs.col);

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

DMatrix WMatrix::operator+(const DMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) //return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.col!=col || rhs.line!=line) //return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);

	DMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=(double)mat[i][j]+rhs.mat[i][j];
	}
	return n;
}

DMatrix WMatrix::operator-(const DMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) //return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.col!=col || rhs.line!=line) //return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);

	DMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=(double)mat[i][j]-rhs.mat[i][j];
	}
	return n;
}

DMatrix WMatrix::operator*(const DMatrix &rhs){
	if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL)// return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.line!=col) //return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);

	DMatrix n(rhs.col, line);
	for(register int j=0;j<n.line;j++) {
		for(register int i=0;i<n.col;i++) {
			register double sum=0.0;
			for(register int k=0;k<col;k++)
				sum+=(double)mat[j][k]*rhs.mat[k][i];
			n.mat[j][i]=sum;
		}
	}
	return n;
}

ZMatrix WMatrix::operator+(const ZMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.col!=col || rhs.line!=line) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j].re=mat[i][j]+rhs.mat[i][j].re;
	}
	return n;
}

ZMatrix WMatrix::operator-(const ZMatrix &rhs){
	if (col==0 || line==0 || mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.col!=col || rhs.line!=line) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j].re=mat[i][j]-rhs.mat[i][j].re;
	}
	return n;
}

ZMatrix WMatrix::operator*(const ZMatrix &rhs){
	if (col==0 || line==0 || mat==NULL || rhs.col==0 || rhs.line==0 || rhs.mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (rhs.line!=col) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZMatrix n(rhs.col, line);

	for(register int j=0;j<n.line;j++) {
		for(register int i=0;i<n.col;i++) {
			register DCplx sum(0.0,0.0);
			for(register int k=0;k<col;k++)
				sum+=mat[j][k]*rhs.mat[k][i];
			n.mat[j][i]=sum;
		}
	}
	return n;
}

WMatrix WMatrix::operator*(const int &rhs){
	if (col==0 || line==0 || mat==NULL) //return WMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	WMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]*rhs;
	}
	return n;
}

DMatrix WMatrix::operator*(const double rhs){
	if (col==0 || line==0 || mat==NULL) //return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=double(mat[i][j])*rhs;
	}
	return n;
}

ZMatrix WMatrix::operator*(const DCplx &rhs){
	if (col==0 || line==0 || mat==NULL) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	ZMatrix n(line, col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++)
			n.mat[i][j]=mat[i][j]*rhs;
	}
	return n;
}

DMatrix WMatrix::operator/(const double rhs){
	if(rhs==0.0) //return DMatrix();
		throw CUtilitisException(__FILE__, __LINE__, DIVISION_BY_ZERO);
	double irhs=1/rhs;
	return (*this)*irhs;
}

ZMatrix WMatrix::operator/(const DCplx &rhs){
	if(rhs==0.0) //return ZMatrix();
		throw CUtilitisException(__FILE__, __LINE__, DIVISION_BY_ZERO);
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
	if(mat==NULL) throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY);
	for(register int i = 0;i<line;i++) {
		mat[i]=new int[col];
		if(mat[i]==NULL) { for(i=i-1;i>=0;i++) delete[] mat[i]; delete[] mat; throw CUtilitisException(__FILE__, __LINE__, OUT_OF_MEMORY); }
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

int WMatrix::sum() {
	if (col==0 || line==0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	int sum=0;
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			sum+=mat[i][j];
		}
	}
	return sum;
}

WMatrix WMatrix::linesum(){
	if (col==0 || line==0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	WMatrix sum(line,1,0);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			sum.mat[i][1]+=mat[i][j];
		}
	}
	return sum;
}

WMatrix WMatrix::colsum(){
	if (col==0 || line==0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	WMatrix sum(1,col,0);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			sum.mat[0][j]+=mat[i][j];
		}
	}
	return sum;
}

ostream& operator<< ( ostream& os, const DCplx& dt)  {
   if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
   if (dt.im>=0.0) 
	   os << ' ' << dt.re << '+' << dt.im << 'i'; 
   else 
	   os << ' ' << dt.re << dt.im<<'i';
   return os;
}

istream& operator>>(istream &is, DCplx &dt) {
   	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	char junk;
	is >> dt.re >> dt.im>> junk ;
	return is;
}

ostream& operator<< ( ostream& os, const DVector& dt )   {
/*	if (dt.vect==NULL || dt.taille==0) {
		os << "DVector is empty !";
		return os;
	}*/
       if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	
	os << dt.taille  << " : [";
	for(register int i=0;i<dt.taille;i++) {
		os << dt.vect[i]<< ' ';
	}
	os << "]" << endl;
	return os;
}

istream &operator>> (istream &is, DVector &v){
   	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	int taille;
	char junk;
	v.~DVector();
	is >> taille >> junk >> junk;
	v=DVector(taille);
	for(register int i=0;i<taille;i++){
		is >> v.vect[i];
	}
	is >> junk;
	return is;
}

ostream& operator<< ( ostream& os, const WVector& dt )  {
/*	if (dt.vect==NULL || dt.taille==0) {
		os << "WVector is empty !";
		return os;
	}
*/
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	os << dt.taille << " : [";

	for(register int i=0;i<dt.taille;i++) {
		os << dt.vect[i]<< ' ';
	}
	os << "]" << endl;

	return os;
}

istream &operator>> (istream &is, WVector &v){
	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	int taille;
	char junk;
	v.~WVector();
	is >> taille >> junk >> junk;
	v=WVector(taille);
	for(register int i=0;i<taille;i++){
		is >> v.vect[i];
	}
	is >> junk;
	return is;
}

ostream& operator<< ( ostream& os, const ZVector& dt )  {
/*	if (dt.vect==NULL || dt.taille==0) {
		os << "ZVector is empty !";
		return os;
	}
*/
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	os << dt.taille << " : [";
	for(register int i=0;i<dt.taille;i++) {
		os << dt.vect[i]<< ' ';
	}
	os << "]" << endl;
	return os;
}

istream &operator>> (istream &is, ZVector &v){
	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	int taille;
	char junk;
	v.~ZVector();
	is >> taille >> junk >> junk;
	v=ZVector(taille);
	for(register int i=0;i<taille;i++){
		is >> v.vect[i];
	}
	is >> junk;
	return is;
}


ostream& operator<< ( ostream& os, const DMatrix& dt )  {
	/*if (dt.mat==NULL || dt.col==0 || dt.line==0) {
		os << "DMatrix is empty !";
		return os;
	}*/
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	
	os << dt.line << "," << dt.col  << " : [" << endl;
	for(register int i=0;i<dt.line;i++) {
		os << "[ ";
		for(register int j=0;j<dt.col;j++) {
			os << dt.mat[i][j]<< ' ';
		}
		os << " ]" << endl;
	}
	os << "]" << endl;
	return os;
}

istream& operator>>(istream& is, DMatrix& dt ) {
	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	char junk; 
	int line,col;
	dt.~DMatrix();
	is >> line >> junk >> col >> junk >> junk;
	dt=DMatrix(line,col);
	for(register int m=0;m<line;m++){
		is >> junk;
		for(register int n=0;n<col;n++){
			is >> dt.mat[m][n];
		}
		is >> junk;
	}
	is >> junk;
	return is;
}

ostream& operator<< ( ostream& os, const WMatrix& dt )  {
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	os << dt.line << "," << dt.col  << " : [" << endl;
	for(register int i=0;i<dt.line;i++) {
		os << "[ ";
		for(register int j=0;j<dt.col;j++) {
			os << dt.mat[i][j]<< ' ';
		}
		os << " ]" << endl;
	}
	os << "]" << endl;
	return os;
}

istream& operator>>(istream& is, WMatrix& dt ) {
	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	char junk; 
	int line,col;
	dt.~WMatrix();
	is >> line >> junk >> col >> junk >> junk;
	dt=WMatrix(line,col);
	for(register int m=0;m<line;m++){
		is >> junk;
		for(register int n=0;n<col;n++){
			is >> dt.mat[m][n];
		}
		is >> junk;
	}
	is >> junk;
	return is;
}

ostream& operator<< (ostream& os, const ZMatrix& dt )  {
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	os << dt.line << "," << dt.col  << " : [" << endl;
	for(register int i=0;i<dt.line;i++) {
		os << "[ ";
		for(register int j=0;j<dt.col;j++) {
			os << dt.mat[i][j]<< ' ';
		}
		os << " ]" << endl;
	}
	os << "]" << endl;
	return os;
}

istream& operator>>(istream& is, ZMatrix& dt ) {
	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	char junk; 
	int line,col;
	dt.~ZMatrix();
	is >> line >> junk >> col >> junk >> junk;
	dt=ZMatrix(line,col);
	for(register int m=0;m<line;m++){
		is >> junk;
		for(register int n=0;n<col;n++){
			is >> dt.mat[m][n];
		}
		is >> junk;
	}
	is >> junk;
	return is;
}

void ZMatrix::inv() {
	if(mat==NULL || line <=0 || col<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(line!=col) throw CMatrixInversionException(INVALID_MATRIX_DIMENSION);;
	ZMatrix dummy(line, col);
	gaussj(*this, dummy);
}

ZMatrix ZMatrix::inv_nip() {
	if(mat==NULL || line <=0 || col<=0) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(line!=col) throw CMatrixInversionException(INVALID_MATRIX_DIMENSION);;
	ZMatrix inv(*this);
	ZMatrix dummy(line, col);
	gaussj(inv, dummy);
	return inv;
}

DCplx ZMatrix::sum() {
	if (col==0 || line==0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DCplx sum(0.0,0.0);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			sum+=mat[i][j];
		}
	}
	return sum;
}

ZMatrix ZMatrix::linesum(){
	if (col==0 || line==0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	ZMatrix sum(line,1);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			sum.mat[i][1]+=mat[i][j];
		}
	}
	return sum;
}

ZMatrix ZMatrix::colsum(){
	if (col==0 || line==0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	ZMatrix sum(1,col);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			sum.mat[0][j]+=mat[i][j];
		}
	}
	return sum;
}


/* 
	Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a.mat[1..n][1..n]
	is the input matrix. b.mat[1..n][1..m] is input containing the m right-hand side vectors. On
	output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of solution
	vectors.
*/
void gaussj(const DMatrix &a, const DMatrix &b){

	if(a.col!=a.line || a.col<=0 || b.col<=0) throw CMatrixInversionException(INVALID_MATRIX_DIMENSION);

	int m=b.col, n=a.line;

	//int *indxc,*indxr,*ipiv;
	WVector indxc(n),indxr(n),ipiv(n);
	int i,icol=0,irow=0,j,k,l,ll;
	double big,dum,pivinv;
	//indxc=ivector(1,n); //The integer arrays ipiv, indxr, and indxc are used for bookkeeping on the pivoting. 
	//indxr=ivector(1,n);
	//ipiv=ivector(1,n);
	for (j=0;j<n;j++) 
		ipiv.vect[j]=0;
	for (i=0;i<n;i++) { //This is the main loop over the columns to be reduced. 
		big=0.0;
		for (j=0;j<n;j++) //This is the outer loop of the search for a pivot element. 
			if (ipiv.vect[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv.vect[k] == 0) {
						if (fabs(a.mat[j][k]) >= big) {
							big=fabs(a.mat[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv.vect[k] > 1) throw CMatrixInversionException(SINGULAR_MATRIX_1_ERROR);
				}
			++(ipiv.vect[icol]);
			/*We now have the pivot element, so we interchange rows, if needed, to put the pivot
			element on the diagonal. The columns are not physically interchanged, only relabeled:
			indxc.vect[i], the column of the ith pivot element, is the ith column that is reduced, while
			indxr.vect[i] is the row in which that pivot element was originally located. If indxr.vect[i] 6=
			indxc.vect[i] there is an implied column interchange. With this form of bookkeeping, the
			solution b's will end up in the correct order, and the inverse matrix will be scrambled
			by columns.*/
			if (irow != icol) {
				for (l=0;l<n;l++) swap(a.mat[irow][l],a.mat[icol][l]);
					for (l=0;l<m;l++) swap(b.mat[irow][l],b.mat[icol][l]);
			}
			indxr.vect[i]=irow; //We are now ready to divide the pivot row by the pivot element, located at irow and icol. 
			indxc.vect[i]=icol;
			if (a.mat[icol][icol] == 0.0) throw CMatrixInversionException(SINGULAR_MATRIX_2_ERROR);
			pivinv=1.0/a.mat[icol][icol];
			a.mat[icol][icol]=1.0;
			for (l=0;l<n;l++) a.mat[icol][l] *= pivinv;
			for (l=0;l<m;l++) b.mat[icol][l] *= pivinv;

			for (ll=0;ll<n;ll++) //Next, we reduce the rows...
				if (ll != icol) { //...except for the pivot one, of course.
					dum=a.mat[ll][icol];
					a.mat[ll][icol]=0.0;
					for (l=0;l<n;l++) a.mat[ll][l] -= a.mat[icol][l]*dum;
					for (l=0;l<m;l++) b.mat[ll][l] -= b.mat[icol][l]*dum;
				}
	}
	/*This is the end of the main loop over columns of the reduction. It only remains to unscramble
	the solution in view of the column interchanges. We do this by interchanging pairs of
	columns in the reverse order that the permutation was built up.*/
	for (l=n-1;l>=0;l--) {
		if (indxr.vect[l] != indxc.vect[l])
			for (k=0;k<n;k++)
				swap(a.mat[k][indxr.vect[l]],a.mat[k][indxc.vect[l]]);
	} //And we are done.
//	free_ivector(ipiv,1,n);
//	free_ivector(indxr,1,n);
//	free_ivector(indxc,1,n);
}

void gaussj(const ZMatrix &a, const ZMatrix &b){

	if(a.col!=a.line || a.col<=0 || b.col<=0) throw CMatrixInversionException(INVALID_MATRIX_DIMENSION);

	int m=b.col, n=a.line;

	//int *indxc,*indxr,*ipiv;
	WVector indxc(n),indxr(n),ipiv(n);
	int i,icol=0,irow=0,j,k,l,ll;
	double big; DCplx pivinv ,dum;
	//indxc=ivector(1,n); //The integer arrays ipiv, indxr, and indxc are used for bookkeeping on the pivoting. 
	//indxr=ivector(1,n);
	//ipiv=ivector(1,n);
	for (j=0;j<n;j++) 
		ipiv.vect[j]=0;
	for (i=0;i<n;i++) { //This is the main loop over the columns to be reduced. 
		big=0.0;
		for (j=0;j<n;j++) //This is the outer loop of the search for a pivot element. 
			if (ipiv.vect[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv.vect[k] == 0) {
						if (a.mat[j][k].mag() >= big) {
							big=a.mat[j][k].mag();
							irow=j;
							icol=k;
						}
					} else if (ipiv.vect[k] > 1) throw CMatrixInversionException(SINGULAR_MATRIX_1_ERROR);
				}
			++(ipiv.vect[icol]);
			/*We now have the pivot element, so we interchange rows, if needed, to put the pivot
			element on the diagonal. The columns are not physically interchanged, only relabeled:
			indxc.vect[i], the column of the ith pivot element, is the ith column that is reduced, while
			indxr.vect[i] is the row in which that pivot element was originally located. If indxr.vect[i] 6=
			indxc.vect[i] there is an implied column interchange. With this form of bookkeeping, the
			solution b's will end up in the correct order, and the inverse matrix will be scrambled
			by columns.*/
			if (irow != icol) {
				for (l=0;l<n;l++) swap(a.mat[irow][l],a.mat[icol][l]);
					for (l=0;l<m;l++) swap(b.mat[irow][l],b.mat[icol][l]);
			}
			indxr.vect[i]=irow; //We are now ready to divide the pivot row by the pivot element, located at irow and icol. 
			indxc.vect[i]=icol;
			//if (a.mat[icol][icol] == 0.0) throw CMatrixInversionException(SINGULAR_MATRIX_2_ERROR);
			if (a.mat[icol][icol].re <= 1e-20 && a.mat[icol][icol].im <= 1e-20 && 
				a.mat[icol][icol].re >= -1e-20 && a.mat[icol][icol].im >= -1e-20) 
				throw CMatrixInversionException(SINGULAR_MATRIX_2_ERROR);
			pivinv=1.0/a.mat[icol][icol];
			a.mat[icol][icol]=1.0;
			for (l=0;l<n;l++) a.mat[icol][l] *= pivinv;
			for (l=0;l<m;l++) b.mat[icol][l] *= pivinv;

			for (ll=0;ll<n;ll++) //Next, we reduce the rows...
				if (ll != icol) { //...except for the pivot one, of course.
					dum=a.mat[ll][icol];
					a.mat[ll][icol]=0.0;
					for (l=0;l<n;l++) a.mat[ll][l] -= a.mat[icol][l]*dum;
					for (l=0;l<m;l++) b.mat[ll][l] -= b.mat[icol][l]*dum;
				}
	}
	/*This is the end of the main loop over columns of the reduction. It only remains to unscramble
	the solution in view of the column interchanges. We do this by interchanging pairs of
	columns in the reverse order that the permutation was built up.*/
	for (l=n-1;l>=0;l--) {
		if (indxr.vect[l] != indxc.vect[l])
			for (k=0;k<n;k++)
				swap(a.mat[k][indxr.vect[l]],a.mat[k][indxc.vect[l]]);
	} //And we are done.
//	free_ivector(ipiv,1,n);
//	free_ivector(indxr,1,n);
//	free_ivector(indxc,1,n);
}


/*
	Given a matrix a.vect[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
	permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
	indx.vect[1..n] is an output vector that records the row permutation eected by the partial
	pivoting; d is output as .1 depending on whether the number of row interchanges was even
	or odd, respectively. This routine is used in combination with lubksb to solve linear equations
	or invert a matrix.
*/

void ludcmp(const DMatrix &a, const WVector &indx, double  *d){
	if (a.col!=a.line || a.col<=0) throw CMatrixInversionException(INVALID_MATRIX_DIMENSION);
	int n = a.col;
	int i,imax=0,j,k;
	double big,dum,sum,temp;
	DVector vv(n); //vv stores the implicit scaling of each row.
	//vv=vector(1,n);
	*d=1.0; //No row interchanges yet.
	for (i=0;i<n;i++) { //Loop over rows to get the implicit scaling information.
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a.mat[i][j])) > big) big=temp;
			if (big == 0.0) throw CMatrixInversionException(SINGULAR_MATRIX_ERROR);
			//No nonzero largest element.
			vv.vect[i]=1.0/big; //Save the scaling.
	}
	for (j=0;j<n;j++) { //This is the loop over columns of Crout's method.
		for (i=0;i<j;i++) { //This is equation (2.3.12) except for i = j.
			sum=a.mat[i][j];
			for (k=0;k<i;k++) sum -= a.mat[i][k]*a.mat[k][j];
			a.mat[i][j]=sum;
		}
		big=0.0; //Initialize for the search for largest pivot element.
		for (i=j;i<n;i++) { //This is i = j of equation (2.3.12) and i = j+1 : ::N of equation (2.3.13). 
			sum=a.mat[i][j];
			for (k=0;k<j;k++)
				sum -= a.mat[i][k]*a.mat[k][j];
			a.mat[i][j]=sum;
			if ( (dum=vv.vect[i]*fabs(sum)) >= big) {
				//Is the figure of merit for the pivot better than the best so far?
				big=dum;
				imax=i;
			}
		}
		if (j != imax) { //Do we need to interchange rows?
			for (k=0;k<n;k++) { //Yes, do so...
				dum=a.mat[imax][k];
				a.mat[imax][k]=a.mat[j][k];
				a.mat[j][k]=dum;
			}
			*d = -(*d); //...and change the parity of d.
			vv.vect[imax]=vv.vect[j]; //Also interchange the scale factor.
		}
		indx.vect[j]=imax;
		if (a.mat[j][j] == 0.0) a.mat[j][j]=TINY;
		/*
		If the pivot element is zero the matrix is singular (at least to the precision of the
		algorithm). For some applications on singular matrices, it is desirable to substitute
		TINY for zero.
		*/
		if (j != n-1) { //Now, finally, divide by the pivot element.
			dum=1.0/(a.mat[j][j]);
			for (i=j+1;i<n;i++) a.mat[i][j] *= dum;
		}
	} 
	//Go back for the next column in the reduction.
	//free_vector(vv,1,n);
}

void lubksb(const DMatrix &a, const WVector &indx, const DVector &b){
/*
	Solves the set of n linear equations A.X = B. Here a.mat[1..n][1..n] is input, not as the matrix
	A but rather as its LU decomposition, determined by the routine ludcmp. indx.vect[1..n] is input
	as the permutation vector returned by ludcmp. b.vect[1..n] is input as the right-hand side vector
	B, and returns with the solution vector X. a, n, and indx are not modified by this routine
	and can be left in place for successive calls with difierent right-hand sides b. This routine takes
	into account the possibility that b will begin with many zero elements, so it is e.cient for use
	in matrix inversion.
*/
	int n=a.line;
	int i,ii=-1,ip,j;
	double sum;
	for (i=0;i<n;i++) { 
	/*
		When ii is set to a positive value, it will become the index of the first nonvanishing element of b. We now
		do the forward substitution, equation (2.3.6). The only new wrinkle is to unscramble the permutation as we go.
	*/
		ip=indx.vect[i];
		sum=b.vect[ip];
		b.vect[ip]=b.vect[i];
		if (ii!=-1)
			for (j=ii;j<=i;j++) sum -= a.mat[i][j]*b.vect[j];
		else if (sum) ii=i; //A nonzero element was encountered, so from now on we will have to do the sums in the loop above. 
		b.vect[i]=sum;
	}
	for (i=n-1;i>=0;i--) { //Now we do the backsubstitution, equation (2.3.7).
		sum=b.vect[i];
		for (j=i+1;j<=n-1;j++) sum -= a.mat[i][j]*b.vect[j];
		b.vect[i]=sum/a.mat[i][i]; //Store a component of the solution vector X.
	}// All done!
}


double convertsnrtosigma2(double snr){
	return 0.5*pow(10, -snr*0.1);
}

double convertsigma2tosnr(double sigma2){
	return 10.0*log10(0.5/sigma2);
}


double convertsnrtosigma2(double Eb, double snr){
	return 0.5*Eb*pow(10, -snr*0.1);
}

double convertsigma2tosnr(double Eb, double sigma2){
	return 10.0*log10(0.5*Eb/sigma2);
}


DMatrix operator/(const double lhs, const DMatrix &rhs) {
	if(rhs.col<=0 || rhs.line<=0||rhs.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	DMatrix n(rhs.line, rhs.col);
	if(lhs ==0) return  n;
	for(register int i=0;i<rhs.line;i++) {
		for(register int j=0;j<rhs.line;j++) {
			 if(rhs.mat[i][j]==0) throw CUtilitisException(__FILE__, __LINE__, DIVISION_BY_ZERO);
			 n.mat[i][j]=lhs/rhs.mat[i][j];
		}
	}
	return n;
}

WVector WVector::operator^(const WVector& v) {
	if(v.taille<=0 || taille<=0 || v.vect==NULL || vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(v.taille!=taille) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	WVector n(taille);
	for(register int i=0;i<taille;i++) {
		n.vect[i]=v.vect[i]^vect[i];
	}
	return n;
}

WVector WVector::operator~() {
	if(taille<=0 || vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	WVector n(taille);
	for(register int i=0;i<taille;i++) {
		n.vect[i]=~vect[i];
	}
	return n;
}

int WVector::sum() {
	if(taille<=0 || vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	int somme=0;
	for(register int i=0;i<taille;i++) {
		somme+=vect[i];
	}
	return somme;
}

int WVector::diff_count(const WVector& v) {
	if(v.taille<=0 || taille<=0 || v.vect==NULL || vect==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if(v.taille!=taille) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);
	
	int cnt=0;
	for(register int i=0;i<taille;i++) {
		cnt+=v.vect[i]==vect[i] ? 0:1;
	}
	return cnt;
}


DVector DVector::operator *(const DMatrix &d) const {
	if (taille<=0 || vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (d.line <=0 || d.col<=0 || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(taille != d.line) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	DVector n(d.col, 0.0);
	for(register int i=0; i< d.col;i++) {
		for(register int j=0;j<taille;j++)  {
			n.vect[i]+=vect[j]*d.mat[j][i];
		}
	}
	return n;
}

DVector DVector::operator *(const WMatrix &d) const {
	if (taille<=0 || vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (d.line <=0 || d.col<=0 || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(taille != d.line) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	DVector n(d.col, 0.0);
	for(register int i=0; i< d.col;i++) {
		for(register int j=0;j<taille;j++)  {
			n.vect[i]+=vect[j]*d.mat[j][i];
		}
	}
	return n;
}

ZVector DVector::operator *(const ZMatrix &d) const {
	if (taille<=0 || vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (d.line <=0 || d.col<=0 || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(taille != d.line) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZVector n(d.col);
	
	for(register int i=0; i< d.col;i++) {
		n.vect[i]=DCplx(0.0,0.0);
		for(register int j=0;j<taille;j++)  {
			n.vect[i]+=vect[j]*d.mat[j][i];
		}
	}
	return n;
}


ZVector ZVector::operator *(const DMatrix &d) const {
	if (taille<=0 || vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (d.line <=0 || d.col<=0 || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(taille != d.line) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZVector n(d.col);
	for(register int i=0; i< d.col;i++) {
		n.vect[i]=0;;
		for(register int j=0;j<taille;j++)  {
			n.vect[i]+=vect[j]*d.mat[j][i];
		}
	}
	return n;
}

ZVector ZVector::operator *(const ZMatrix &d) const {
	if (taille<=0 || vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (d.line <=0 || d.col<=0 || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(taille != d.line) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZVector n(d.col);
	
	for(register int i=0; i< d.col;i++) {
		n.vect[i]=DCplx(0.0,0.0);
		for(register int j=0;j<taille;j++)  {
			n.vect[i]+=vect[j]*d.mat[j][i];
		}
	}
	return n;
}

DVector WVector::operator *(const DMatrix &d) const {
	if (taille<=0 || vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (d.line <=0 || d.col<=0 || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(taille != d.line) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	DVector n(d.col, 0.0);
	for(register int i=0; i< d.col;i++) {
		for(register int j=0;j<taille;j++)  {
			n.vect[i]+=vect[j]*d.mat[j][i];
		}
	}
	return n;
}

WVector WVector::operator *(const WMatrix &d) const {
	if (taille<=0 || vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (d.line <=0 || d.col<=0 || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(taille != d.line) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	WVector n(d.col, 0);
	for(register int i=0; i< d.col;i++) {
		for(register int j=0;j<taille;j++)  {
			n.vect[i]+=vect[j]*d.mat[j][i];
		}
	}
	return n;
}

ZVector WVector::operator *(const ZMatrix &d) const {
	if (taille<=0 || vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (d.line <=0 || d.col<=0 || d.mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(taille != d.line) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_MATRIX_DIMENSION);
	ZVector n(d.col);
	
	for(register int i=0; i< d.col;i++) {
		n.vect[i]=DCplx(0.0,0.0);
		for(register int j=0;j<taille;j++)  {
			n.vect[i]+=vect[j]*d.mat[j][i];
		}
	}
	return n;
}

DVector DMatrix::operator *(const DVector &d) const {
	if (d.taille<=0 || d.vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (line <=0 || col<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (d.taille != col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	DVector n(line, 0);

	for(register int i= 0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			n.vect[i]+=mat[i][j]*d.vect[j];
		}
	}
	return n;
}

ZVector DMatrix::operator *(const ZVector &d) const {
	if (d.taille<=0 || d.vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (line <=0 || col<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (d.taille != col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	ZVector n(line);

	for(register int i= 0;i<line;i++) {
		n.vect[i]=0;
		for(register int j=0;j<col;j++) {
			n.vect[i]+=mat[i][j]*d.vect[j];
		}
	}
	return n;
}

DVector DMatrix::operator *(const WVector &d) const {
	if (d.taille<=0 || d.vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (line <=0 || col<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (d.taille != col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	DVector n(line, 0);

	for(register int i= 0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			n.vect[i]+=mat[i][j]*d.vect[j];
		}
	}
	return n;
}

DVector WMatrix::operator *(const DVector &d) const {
	if (d.taille<=0 || d.vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (line <=0 || col<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (d.taille != col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	DVector n(line, 0);

	for(register int i= 0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			n.vect[i]+=mat[i][j]*d.vect[j];
		}
	}
	return n;
}

ZVector WMatrix::operator *(const ZVector &d) const {
	if (d.taille<=0 || d.vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (line <=0 || col<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (d.taille != col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	ZVector n(line);

	for(register int i= 0;i<line;i++) {
		n.vect[i]=0;
		for(register int j=0;j<col;j++) {
			n.vect[i]+=mat[i][j]*d.vect[j];
		}
	}
	return n;
}

WVector WMatrix::operator *(const WVector &d) const {
	if (d.taille<=0 || d.vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (line <=0 || col<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (d.taille != col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	WVector n(line, 0);

	for(register int i= 0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			n.vect[i]+=mat[i][j]*d.vect[j];
		}
	}
	return n;
}

ZVector ZMatrix::operator *(const DVector &d) const {
	if (d.taille<=0 || d.vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (line <=0 || col<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (d.taille != col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	ZVector n(line);

	for(register int i= 0;i<line;i++) {
		n.vect[i]=0;
		for(register int j=0;j<col;j++) {
			n.vect[i]+=mat[i][j]*d.vect[j];
		}
	}
	return n;
}

ZVector ZMatrix::operator *(const ZVector &d) const {
	if (d.taille<=0 || d.vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (line <=0 || col<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (d.taille != col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	ZVector n(line);

	for(register int i= 0;i<line;i++) {
		n.vect[i]=0;
		for(register int j=0;j<col;j++) {
			n.vect[i]+=mat[i][j]*d.vect[j];
		}
	}
	return n;
}

ZVector ZMatrix::operator *(const WVector &d) const {
	if (d.taille<=0 || d.vect == NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	if (line <=0 || col<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if (d.taille != col) throw CUtilitisException(__FILE__, __LINE__, INCOMPATIBLE_VECTOR_DIMENSION);

	ZVector n(line);

	for(register int i= 0;i<line;i++) {
		n.vect[i]=0;
		for(register int j=0;j<col;j++) {
			n.vect[i]+=mat[i][j]*d.vect[j];
		}
	}
	return n;
}

// all the functions below are left as is for backward compatibility 
// with C version of utilitis developped by V. Meghdadi prior to 
// this version (http://www.ensil.unilim.fr/~meghdadi)

void Err_Message(char Msg[]){
	cerr << Msg << endl;
	exit (2);
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
	if (x.taille*2<data.taille) {
		cerr << endl << "Incompatible sizes in QPSK_MOD function" << endl;

		return;

	}
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
	 cerr << "File ReadData: Cannot open file "<< Fname<<endl;
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
	 cerr << "File ReadData: Cannot open file "<< Fname<<endl;
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
	 cerr << "File ReadData: Cannot open file "<< Fname<<endl;
	 exit (1);
  }
/* Read the values */
  while (!feof(fp)){
	fscanf(fp,"%ld\n",&tmp);
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
	 cerr << "File ReadData: Cannot open file "<< Fname<<endl;
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

	for (i=0; i<a.taille; i++)
		fprintf(fp,"%lf  %lf\n",a.vect[i].re, a.vect[i].im);
	fclose(fp);
	return 0;
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
	int min = a.taille> b.taille ? b.taille : a.taille; 
	for (i=0 ; i<min; i++)
		if (*p_a++ != *p_b++) cmpt++;
	return cmpt;
}


DCplx ZMatrix::trace() {
	if(line<=0 || col<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(line!=col) throw CUtilitisException(__FILE__, __LINE__, MATRIX_NOT_SQUARE);
	DCplx c(0,0);
	for(register int i=0;i<line;i++) {
		c+=mat[i][i];
	}
	return c;
}

int WMatrix::trace() {
	if(line<=0 || col<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(line!=col) throw CUtilitisException(__FILE__, __LINE__, MATRIX_NOT_SQUARE);
	int c = 0;
	for(register int i=0;i<line;i++) {
		c+=mat[i][i];
	}
	return c;
}

double DMatrix::trace() {
	if(line<=0 || col<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	if(line!=col) throw CUtilitisException(__FILE__, __LINE__, MATRIX_NOT_SQUARE);
	double  c = 0.0;
	for(register int i=0;i<line;i++) {
		c+=mat[i][i];
	}
	return c;
}



void WMatrix::reset(int val) {
	if(line<=0 || col<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			mat[i][j]=val;
		}
	}


}

int WMatrix::max() {
	if(line<=0 || col<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	int m = mat[0][0];
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			if(mat[i][j]>m) m = mat[i][j];
		}
	}
	return m;
}

int WMatrix::min() {
	if(line<=0 || col<=0 || mat==NULL) throw CUtilitisException(__FILE__, __LINE__, EMPTY_MATRIX);
	int m = mat[0][0];
	for(register int i=0;i<line;i++) {
		for(register int j=0;j<col;j++) {
			if(mat[i][j]<m) m = mat[i][j];
		}
	}
	return m;
}

double DVector::variance(){
    if(vect==NULL || taille <=0) throw CUtilitisException(EMPTY_MATRIX);
    register double sum=0.0, sum2=0.0;
    double isize=1.0/taille, mean;
    for(register int i=0;i<taille;i++) {
        sum+=vect[i];
        sum2+=vect[i]*vect[i];
    }
    mean=(sum*isize);
    return (sum2*isize-mean*mean);
}

double WVector::variance(){
    if(vect==NULL || taille <=0) throw CUtilitisException(EMPTY_MATRIX);
    register double sum=0.0, sum2=0.0;
    double isize=1.0/taille, mean;
    for(register int i=0;i<taille;i++) {
        sum+=vect[i];
        sum2+=vect[i]*vect[i];
    }
    mean=(sum*isize);

    return (sum2*isize-mean*mean);
}

