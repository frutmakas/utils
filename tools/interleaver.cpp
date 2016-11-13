/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/07/30 14:19:00 $
 * $Revision: 1.1.2.17 $
 * $Id: interleaver.cpp,v 1.1.2.17 2006/07/30 14:19:00 syed Exp $
 ********************************************************************/
#include "tools/interleaver.h"
#include "tools/tools.h"
#include <iostream>
#include<time.h>
using namespace std;
#include "tools/interleaver.h"
//dRandUniStatePtr interleave_seed;

CInterleaver::CInterleaver(int size, int seed) {
    padsize=0;
    if (seed<=0) seed=time(NULL);
    dRandUniInit(&interleaver_seed, seed, 0.0, 1.0);
    int lsize=size;
	WVector flag(size ,0, 'i');
	map = WVector(size);
	int p, tmp;
	for(register int s=0;s<size;s++,lsize--) {
        map.vect[s]= flag.vect[p = (int)(dRandUni(&interleaver_seed)*lsize)];
        // swap flag
        tmp = flag.vect[p];
        flag.vect[p] = flag.vect[size-s-1];
        flag.vect[size-s-1] = tmp;
	}
}

CInterleaver::CInterleaver(const WVector &pmap) {
    if(pmap.vect==NULL || pmap.taille<=0) throw CInterleaverException(BAD_PREDEFINED_MAPPING);

    padsize=0;
    dRandUniInit(&interleaver_seed, time(NULL), 0.0, 1.0);

    int size=pmap.taille;

	WVector flag(size ,0);
	map = WVector(size);

    for(int i=0;i<size;i++) {
        if(pmap.vect[i]>size || pmap.vect[i]<0) throw CInterleaverException(BAD_PREDEFINED_MAPPING);
        if(++flag.vect[pmap.vect[i]]>1) throw CInterleaverException(BAD_PREDEFINED_MAPPING);
        map.vect[i]=pmap.vect[i];
    }
}

CInterleaver::CInterleaver(char type, int line, int col){
    padsize=0;
    dRandUniInit(&interleaver_seed, time(NULL), 0.0, 1.0);


    int a=-1,b=-1, c,  as,bs, ap,bp, is, js, vmax, i2, imin;
    int I=line, J=col;

    WMatrix m(I,J); WVector v(I,0,'i');

    //find b .. b must be relative prime to I and b = floor(I/J)
    vmax = I/J;
    for(c=vmax;c>1;c--) {
        if (isRelativePrime(c,I)) {b = c; break;}
    }
    if(b==-1) throw CInterleaverException(CREATE_OPTIMIZED_INTERLEAVER_FAILED);//, __LINE__, __FILE__);
    //have b ready ..  get a now
    vmax = I-1;
    for(c=vmax;c>1;c--) {
        if (isRelativePrime(c,I) && isRelativePrime(c,b)) {a = c; break;}
    }
    if(a==-1) throw CInterleaverException(CREATE_OPTIMIZED_INTERLEAVER_FAILED);//, __LINE__, __FILE__);
    // have a and b .. calculate ap and bp

    ap = a%I;
    bp = b%I;

    // get the first column
    for(int i=0;i<I;i++) {
        m.mat[i][0]=(v.vect[i]=((a*i)%I))*J;
    }

    imin=-1;

    for(i2=0;i2<I;i2++)
        if (((a*i2)%I)==b) {
            imin =i2;
            break;
        }

    if(imin==-1) throw CInterleaverException(CREATE_OPTIMIZED_INTERLEAVER_FAILED);//, __LINE__, __FILE__);

    int offset;

    for(int j=1;j<J;j++){
        offset=imin*j;
        for(i2=0;i2<I;i2++) {
            m.mat[i2][j] = v.vect[(J*I+i2-offset)%I]*J+j;
        }
    }

    map = WVector(I*J,0);

    for(int i=0,p=0;i<I;i++) {
        for(int j=0;j<J;j++) {
            map.vect[p++]=m.mat[i][j];
        }
    }
}

CInterleaver::CInterleaver(const CInterleaver &src){
   padsize=src.padsize;
    interleaver_seed=src.interleaver_seed;
    map = src.map;
}

CInterleaver::~CInterleaver(){
}


ZVector CInterleaver::Apply(const ZVector &src) {
	int finalsize=src.taille%map.taille;
	padsize=0;
	if(finalsize) {
		padsize=map.taille-finalsize;
	}
	ZVector dest(src.taille+padsize);
	DCplx c(0.0,0.0);
	for(int s=0,d=0;d<dest.taille;d+=map.taille){
        finalsize=src.taille-s>map.taille?map.taille:src.taille-s;
 		for(int m=0;m<finalsize;m++) {
			dest.vect[map.vect[m]+d] = s<src.taille ? src.vect[s++] : c;
		}
	}
	return dest;
}

void CInterleaver::Apply(char *dest, int *destlen, char *src, int srclen ) {
    if(srclen<=0 || src==NULL || dest == NULL || *destlen<srclen)
        throw CInterleaverException(INVALID_RANGE);
	int finalsize=srclen%map.taille;
	padsize=0;
	if(finalsize) {
		padsize=map.taille-finalsize;
	}
    finalsize = srclen+padsize;
    if (*destlen < finalsize) throw CInterleaverException(INVALID_RANGE);
//	dest = new char[*destlen];
    *destlen = finalsize;
	char c=0;
	for(int s=0,d=0;d<*destlen;d+=map.taille){
        finalsize=srclen-s>map.taille?map.taille:srclen-s;
 		for(int m=0;m<finalsize;m++) {
			dest[map.vect[m]+d] = s<srclen ? src[s++] : c;
		}
	}
}

DVector CInterleaver::Apply(const DVector &src) {
	int finalsize=src.taille%map.taille;
	padsize=0;
	if(finalsize) {
		padsize=map.taille-finalsize;
	}
	DVector dest(src.taille+padsize);
	for(int s=0,d=0;d<dest.taille;d+=map.taille){
        finalsize=src.taille-s>map.taille?map.taille:src.taille-s;
 		for(int m=0;m<finalsize;m++) {
			dest.vect[map.vect[m]+d] = s<src.taille ? src.vect[s++] : 0;
		}
	}
	return dest;
}

WVector CInterleaver::Apply(const WVector &src) {
	padsize=0;
	int finalsize=src.taille%map.taille;
	if(finalsize) {
		padsize=map.taille-finalsize;
	}
	WVector dest(src.taille+padsize);
    if(map.taille<=0) throw CInterleaverException(INVALID_INTERLEAVER_STATE);
	for(int s =0,d=0;d<dest.taille;d+=map.taille){
        finalsize=src.taille-s>map.taille?map.taille:src.taille-s;
 		for(int m=0;m<finalsize;m++) {
			dest.vect[map.vect[m]+d] = src.vect[s++];
		}
	}
	return dest;
}

ZVector CInterleaver::Extract(const ZVector &src){
	if (src.taille%map.taille){
		throw CInterleaverException(INVALID_INTERLEAVED_VECTOR_SIZE);
	}
	int unpaddedsize=src.taille-padsize, finalsize;
	ZVector dest(unpaddedsize);
	for(int s=0,d=0;s<src.taille&&d<unpaddedsize;s+=map.taille){
        finalsize=(unpaddedsize-d)> map.taille ? map.taille : unpaddedsize-d;
		for(int m=0;m<finalsize;m++) {
			dest.vect[d++]=src.vect[map.vect[m]+s];
		}
	}
	return dest;
}

DVector CInterleaver::Extract(const DVector &src){
	if (src.taille%map.taille){
		throw CInterleaverException(INVALID_INTERLEAVED_VECTOR_SIZE);
	}
	int unpaddedsize=src.taille-padsize, finalsize;
	DVector dest(unpaddedsize);
	for(int s=0,d=0;s<src.taille&&d<unpaddedsize;s+=map.taille){
        finalsize=(unpaddedsize-d)> map.taille ? map.taille : unpaddedsize-d;
		for(int m=0;m<finalsize;m++) {
			dest.vect[d++]=src.vect[map.vect[m]+s];
		}
	}
	return dest;
}

WVector CInterleaver::Extract(const WVector &src){
	if (src.taille%map.taille){
		throw CInterleaverException(INVALID_INTERLEAVED_VECTOR_SIZE);
	}
	int unpaddedsize=src.taille-padsize, finalsize;
    WVector dest(unpaddedsize);
	for(int s=0,d=0;s<src.taille&&d<unpaddedsize;s+=map.taille){
        finalsize=(unpaddedsize-d)> map.taille ? map.taille : unpaddedsize-d;
		for(int m=0;m<finalsize;m++) {
			dest.vect[d++]=src.vect[map.vect[m]+s];
		}
	}
	return dest;
}


void CInterleaver::Extract(char *dest, int *destlen, char *src, int srclen ){
    if(srclen<=0 || src==NULL || dest == NULL || *destlen<srclen)
        throw CInterleaverException(INVALID_RANGE);

	if (srclen%map.taille){
		throw CInterleaverException(INVALID_INTERLEAVED_VECTOR_SIZE);
	}
	int unpaddedsize=srclen-padsize, finalsize;

    if(*destlen<unpaddedsize)
        throw CInterleaverException(INVALID_RANGE);

//    WVector dest(unpaddedsize);
    *destlen = unpaddedsize;
	for(int s=0,d=0;s<srclen&&d<unpaddedsize;s+=map.taille){
        finalsize=(unpaddedsize-d)> map.taille ? map.taille : unpaddedsize-d;
		for(int m=0;m<finalsize;m++) {
			dest[d++]=src[map.vect[m]+s];
		}
	}
//	return dest;
}
//******** in place operation ********

void CInterleaver::ipApply(const ZVector &src) {
    if(src.taille<=0 || src.vect==NULL) throw CInterleaverException(EMPTY_VECTOR);
	int finalsize=src.taille%map.taille;
	padsize=0;
	if(finalsize) {
		padsize=map.taille-finalsize;
	}
    if(padsize) throw CInterleaverException(CANNOT_PAD_IN_INPLACE_MODE);
    if(map.taille<=0) throw CInterleaverException(INVALID_INTERLEAVER_STATE);

	ZVector tmp(src);
	for(int s=0,d=0;d<src.taille;d+=map.taille){
 		for(int m=0;m<map.taille;m++) {
			src.vect[map.vect[m]+d] = tmp.vect[s++];
		}
	}
}

void CInterleaver::ipApply(const DVector &src) {
    if(src.taille<=0 || src.vect==NULL) throw CInterleaverException(EMPTY_VECTOR);
	int finalsize=src.taille%map.taille;
	padsize=0;
	if(finalsize) {
		padsize=map.taille-finalsize;
	}
    if(padsize) throw CInterleaverException(CANNOT_PAD_IN_INPLACE_MODE);
    if(map.taille<=0) throw CInterleaverException(INVALID_INTERLEAVER_STATE);
	DVector tmp(src);
	for(int s=0,d=0;d<src.taille;d+=map.taille){
 		for(int m=0;m<map.taille;m++) {
			src.vect[map.vect[m]+d] = tmp.vect[s++];
		}
	}
}

void CInterleaver::ipApply(const WVector &src) {
    if(src.taille<=0 || src.vect==NULL) throw CInterleaverException(EMPTY_VECTOR);
	padsize=0;
	int finalsize=src.taille%map.taille;
	if(finalsize) {
		padsize=map.taille-finalsize;
	}
    if(padsize) throw CInterleaverException(CANNOT_PAD_IN_INPLACE_MODE);
    if(map.taille<=0) throw CInterleaverException(INVALID_INTERLEAVER_STATE);
	WVector tmp(src);
	for(int s=0,d=0;d<src.taille;d+=map.taille){
 		for(int m=0;m<map.taille;m++) {
			src.vect[map.vect[m]+d] = tmp.vect[s++];
		}
	}
}

void CInterleaver::ipExtract(const ZVector &src){
    if(src.taille<=0 || src.vect==NULL) throw CInterleaverException(EMPTY_VECTOR);
    if(map.taille<=0 || map.vect==NULL) throw CInterleaverException(INVALID_INTERLEAVER_STATE);
	if (src.taille%map.taille){
		throw CInterleaverException(INVALID_INTERLEAVED_VECTOR_SIZE);
	}
    ZVector tmp(src);
	for(int s=0,d=0;s<src.taille;s+=map.taille){
		for(int m=0;m<map.taille;m++) {
			src.vect[d++]=tmp.vect[map.vect[m]+s];
		}
	}
}

void CInterleaver::ipExtract(const DVector &src){
    if(src.taille<=0 || src.vect==NULL) throw CInterleaverException(EMPTY_VECTOR);
    if(map.taille<=0 || map.vect==NULL) throw CInterleaverException(INVALID_INTERLEAVER_STATE);
	if (src.taille%map.taille){
		throw CInterleaverException(INVALID_INTERLEAVED_VECTOR_SIZE);
	}
    DVector tmp(src);
	for(int s=0,d=0;s<src.taille;s+=map.taille){
		for(int m=0;m<map.taille;m++) {
			src.vect[d++]=tmp.vect[map.vect[m]+s];
		}
	}
}

void CInterleaver::ipExtract(const WVector &src){
    if(src.taille<=0 || src.vect==NULL) throw CInterleaverException(EMPTY_VECTOR);
    if(map.taille<=0 || map.vect==NULL) throw CInterleaverException(INVALID_INTERLEAVER_STATE);
	if (src.taille%map.taille){
		throw CInterleaverException(INVALID_INTERLEAVED_VECTOR_SIZE);
	}
    WVector tmp(src);
	for(int s=0,d=0;s<src.taille;s+=map.taille){
		for(int m=0;m<map.taille;m++) {
			src.vect[d++]=tmp.vect[map.vect[m]+s];
		}
	}
}

//******** end in place operation ********

ostream &operator<<(ostream &os, CInterleaver &b){
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	os << b.map;
	os << b.padsize << endl;
	os << b.interleaver_seed << endl;
	return os;
}

istream &operator>>(istream &is, CInterleaver &b) {
	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	is >> b.map;
	is >> b.padsize;
	is >> b.interleaver_seed;
	return is;
}

int interleaver_selftest(){
    cout << endl << " Running interleaver self test and verification .. please wait " << endl;
    const int VSIZE=32452841, PNSIZE=33;
    WVector v(VSIZE, 0, 'd');
    unsigned long int starttime=time(NULL), stoptime;
    WVector pn(PNSIZE,0,'p'), out, check;
    int diff, cumul=0;
    cout << " Testing interleaver in not in place mode " << endl;

    for(int s=3;s<PNSIZE;s++) {
        CInterleaver i(pn.vect[s], starttime);
        check = i.Extract(i.Apply(v));
        diff = check.diff_count(v)  ;

        if (diff) cout << endl << "With interleaver size " << s << ", the test vector has " << diff  << " differences with interleaver extracted output" << endl;
        cout << ".";
    }
    cout << endl;
    stoptime=time(NULL);
    cout <<endl << "Elapsed time ... " << stoptime-starttime << " seconds" << endl;

    const int VSIZE2 = 39916800; //pn.vect[10]*pn.vect[4]*pn.vect[9]*pn.vect[16]*pn.vect[25];
    v = WVector(VSIZE2,0,'i');
    WVector checkv2 = v;
    cout << " Testing inplace mode interleaver  " << endl;
    cout << "New vector test size : " << VSIZE2 << endl;

    for(int s=1;s<13;s++) {
        int vv = s!= 10 ? s*10 : s*12;
        CInterleaver i(vv, starttime);
        i.ipApply(v) ;
        i.ipExtract(v);
        checkv2 = i.Extract(i.Apply(v));
        cumul += (diff = checkv2.diff_count(v));
        if (diff) cout << endl << "With interleaver size " << 12 << ", the test vector has " << diff  << " differences with interleaver extracted output" << endl;
        cout << ".";
    }
    stoptime=time(NULL);
    cout <<endl << "Elapsed time ... " << stoptime-starttime << " seconds" << endl;


    cout << "done ";
    return cumul ? 0 : 1;
}

bool CInterleaver::operator==(const CInterleaver &c)  const {
    if(map.taille<=0 || map.vect==NULL || c.map.taille<=0 || c.map.vect==NULL) throw CInterleaverException(INVALID_INTERLEAVER_STATE);
    if(map.taille!=c.map.taille) {cerr << "Woi apek" << endl;  return false; }
    return map==c.map;
}

bool CInterleaver::operator!=(const CInterleaver &c)  const {
    if(map.taille<=0 || map.vect==NULL || c.map.taille<=0 || c.map.vect==NULL) throw CInterleaverException(INVALID_INTERLEAVER_STATE);
    if(map.taille!=c.map.taille) return true;
    return map!=c.map;
}

CInterleaver &CInterleaver::operator=(const CInterleaver &c){
   padsize=c.padsize;
    interleaver_seed=c.interleaver_seed;
    map = c.map;
    return *this;
}
