/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/07/21 10:58:00 $
 * $Revision: 1.1.2.22 $
 * $Id: randgen.cpp,v 1.1.2.22 2006/07/21 10:58:00 syed Exp $
 ********************************************************************/
 
#include "globaldef.h"
/*#include <math.h>
#include <stdlib.h>*/
#include <iostream>
#include <fstream>
using namespace std;
#include "rand/randgen.h"

#ifdef WIN_DOS
#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 
#pragma intrinsic(sqrt, pow)
#endif
/*
Returns a uniform random deviate between 0:0 and 1:0. Set idum to any negative value to
initialize or reinitialize the sequence.
*/
void dRandUniInit(dRandUniStatePtr *statePtr, int seed, double min, double max){
	statePtr->idum=seed>0 ? -seed : seed;
	statePtr->iff=0;
	statePtr->min=min>max ? max : min;
	statePtr->max=min>max ? min : max;
	statePtr->dist=max-min;
	statePtr->MSEED=MSEED_MAX-(abs(seed%16384));
	dRandUni(statePtr);
}	

double dRandUni(dRandUniStatePtr *statePtr){
/*	static int statePtr->inext,statePtr->inextp;
	static long statePtr->ma[56]; // The value 56 (range statePtr->ma[1..55]) is special and should not be modied; see Knuth. 
	static int iff=0; */
	long mj,mk;
	register int i,ii,k;
	if (statePtr->idum < 0 || statePtr->iff == 0) { // Initialization.
		statePtr->iff=1;
		mj=labs(statePtr->MSEED-labs(statePtr->idum)); //Initialize statePtr->ma[55] using the seed idum and the	large number MSEED. 
		mj %= MBIG;
		statePtr->ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) { // Now initialize the rest of the table,in a slightly random order,with numbers that are not especially random.
			ii=(21*i) % 55; 
			statePtr->ma[ii]=mk; 
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=statePtr->ma[ii];
		}

		for (k=1;k<=4;k++) {//We randomize them by warming up the generator.
			for (i=1;i<=55;i++) {
				statePtr->ma[i] -= statePtr->ma[1+(i+30) % 55];
				if (statePtr->ma[i] < MZ) statePtr->ma[i] += MBIG;
			}
		}
		statePtr->inext=0; //Prepare indices for our rst generated number.
		statePtr->inextp=31; //The constant 31 is special; see Knuth.
		statePtr->idum=1;
	}

	if (++(statePtr->inext) == 56) statePtr->inext=1; // Increment statePtr->inext and statePtr->inextp, wrapping around 56 to 1. 
	if (++(statePtr->inextp) == 56) statePtr->inextp=1;
	mj=statePtr->ma[statePtr->inext]-statePtr->ma[statePtr->inextp]; // Generate a new random number subtractively.
	if (mj < MZ) mj += MBIG; // Be sure that it is in range.
	statePtr->ma[statePtr->inext]=mj; //Store it,
	return (mj*FAC)*(statePtr->dist)+statePtr->min; //and output the derived uniform deviate.
}



void dbRandUni(dRandUniStatePtr *statePtr, const DVector &dv){
	long mj,mk;
	register int i,ii,k;

	//dv.vect[i]=dRandUni(statePtr);
	if (statePtr->idum < 0 || statePtr->iff == 0) { // Initialization.
		statePtr->iff=1;
		mj=labs(statePtr->MSEED-labs(statePtr->idum)); //Initialize statePtr->ma[55] using the seed idum and the	large number MSEED. 
		mj %= MBIG;
		statePtr->ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) { // Now initialize the rest of the table,in a slightly random order,with numbers that are not especially random.
			ii=(21*i) % 55; 
			statePtr->ma[ii]=mk; 
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=statePtr->ma[ii];
		}

		for (k=1;k<=4;k++) {//We randomize them by warming up the generator.
			for (i=1;i<=55;i++) {
				statePtr->ma[i] -= statePtr->ma[1+(i+30) % 55];
				if (statePtr->ma[i] < MZ) statePtr->ma[i] += MBIG;
			}
		}
		statePtr->inext=0; //Prepare indices for our rst generated number.
		statePtr->inextp=31; //The constant 31 is special; see Knuth.
		statePtr->idum=1;
	}

	for (int vi =0; vi < dv.taille;vi++) {
		if (++(statePtr->inext) == 56) statePtr->inext=1; // Increment statePtr->inext and statePtr->inextp, wrapping around 56 to 1. 
		if (++(statePtr->inextp) == 56) statePtr->inextp=1;
		mj=statePtr->ma[statePtr->inext]-statePtr->ma[statePtr->inextp]; // Generate a new random number subtractively.
		if (mj < MZ) mj += MBIG; // Be sure that it is in range.
		statePtr->ma[statePtr->inext]=mj; //Store it,
		dv.vect[vi]=(mj*FAC)*(statePtr->dist)+statePtr->min; //and output the derived uniform deviate.
	}
}

void dRandGausInit(dRandGausStatePtr *statePtr, double mean, double variance){
	if (statePtr->unistate.idum == 0) statePtr->unistate.idum=-(int)ceil(fabs((mean-variance)*(mean+variance)));
	statePtr->iset=0;
	statePtr->mean=mean;
	statePtr->variance=variance;
	statePtr->sigma=sqrt(variance);//*variance;

	dRandUniInit(&(statePtr->unistate));
}	

void dRandGausInit(dRandGausStatePtr *statePtr, int seed, double mean, double variance){
	if (statePtr->unistate.idum == 0) statePtr->unistate.idum=seed<0?seed:-seed;
	statePtr->iset=0;
	statePtr->mean=mean;
	statePtr->variance=variance;
	statePtr->sigma=sqrt(variance);//*variance;

	dRandUniInit(&(statePtr->unistate), seed);
}

/*
Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
as the source of uniform deviates.
*/
double dRandGaus(dRandGausStatePtr *statePtr) {
	double fac,rsq,v1,v2;
	if (statePtr->unistate.idum < 0) statePtr->iset=0; // Reinitialize.
	if (statePtr->iset == 0) { // We don't have an extra deviate handy, so pick two uniform numbers in the square extending
		do {
			v1=2.0*dRandUni(&(statePtr->unistate))-1.0; //from -1 to +1 in each direction, 
			v2=2.0*dRandUni(&(statePtr->unistate))-1.0; 
			rsq=v1*v1+v2*v2; // see if they are in the unit circle,
		} while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try again.
		fac=statePtr->sigma*sqrt(-2.0*log(rsq)/rsq);
		/* Now make the Box-Muller transformation to get two normal deviates. Return one and save the other for next time.	*/
		statePtr->gset=v1*fac+statePtr->mean;
		statePtr->iset=1; // Set flag.
		return v2*fac+statePtr->mean;
	} else { // We have an extra deviate handy,
		statePtr->iset=0; //so unset the flag,
		return statePtr->gset; //and return it.
	}
}

//#define RAND_ANALYSIS
int dbRandGaus(dRandGausStatePtr *statePtr, const DVector &dv) {
	if (dv.taille==0 || dv.vect==NULL) return 0;
	double fac,rsq,v1,v2;
	statePtr->iset=0;
#ifdef RAND_ANALYSIS
    int miss = 0, aim=0;
#endif
	for(register int i=0;i<dv.taille;i+=2) {
#ifdef RAND_ANALYSIS
        aim=-1;
#endif
		do {
#ifdef RAND_ANALYSIS
            aim++;
#endif
			v1=2.0*dRandUni(&(statePtr->unistate))-1.0; //from -1 to +1 in each direction,
			v2=2.0*dRandUni(&(statePtr->unistate))-1.0;
			rsq=v1*v1+v2*v2; // see if they are in the unit circle,
		} while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try again.
#ifdef RAND_ANALYSIS
        miss+=aim;
#endif
		fac=statePtr->sigma*sqrt(-2.0*log(rsq)/rsq);
//		dv.vect[i]=statePtr->sigma*v2*fac+statePtr->mean;
		dv.vect[i]=v2*fac+statePtr->mean;
		if(i+1<dv.taille)
			dv.vect[i+1]=v1*fac+statePtr->mean;
		else {
			statePtr->iset=1;
			statePtr->gset=v1*fac+statePtr->mean;
		}
	}
#ifdef RAND_ANALYSIS
    cerr << endl << "randgen.cpp > DEBUG : Missed = " << miss << " out of " << dv.taille << endl;
#endif

	return 1;
}

int dmRandGaus(dRandGausStatePtr *statePtr, const DMatrix &dm) {
	if(dm.line<=0 || dm.col <=0 || dm.mat==0) throw CUtilitisException(EMPTY_MATRIX);
	
	for(register int i=0;i<dm.line;i++) {
		for(register int j=0;j<dm.col;j++) {
			dm.mat[i][j]=dRandGaus(statePtr);
		}
	}
	return 1;
}

int zmRandGaus(zRandGausStatePtr *statePtr, const ZMatrix &dm) {
	if(dm.line<=0 || dm.col <=0 || dm.mat==0) throw CUtilitisException(EMPTY_MATRIX);
	
	for(register int i=0;i<dm.line;i++) {
		for(register int j=0;j<dm.col;j++) {
			dm.mat[i][j]=zRandGaus(statePtr);
		}
	}
	return 1;
}

int dbApplyRandGaus(dRandGausStatePtr *statePtr, const DVector &dv) {
	if (dv.taille==0 || dv.vect==NULL) return 0;
	for(register int i = 0; i<dv.taille;i++) {
		dv.vect[i]+=dRandGaus(statePtr);
	}
	return 1;

}


void zRandGausInit(zRandGausStatePtr *statePtr, double mean, double variance){
	if (statePtr->unistate.idum == 0) statePtr->unistate.idum=-(int)ceil(fabs((mean-variance)*(mean+variance)));
	statePtr->mean=mean;
	statePtr->variance=fabs(0.5*variance);
	statePtr->sigma=sqrt(0.5*variance);//*variance;

	dRandUniInit(&(statePtr->unistate));
}	

void zRandGausInit(zRandGausStatePtr *statePtr, double mean, double variance, int seed){
	if (statePtr->unistate.idum == 0) statePtr->unistate.idum=-(int)ceil(fabs((mean-variance)*(mean+variance)));
	statePtr->mean=mean;
	statePtr->variance=fabs(0.5*variance);
	statePtr->sigma=sqrt(0.5*variance);//*variance;

	dRandUniInit(&(statePtr->unistate), seed);
}	
/*
Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
as the source of uniform deviates.
*/
DCplx zRandGaus(zRandGausStatePtr *statePtr) {
	double fac,rsq,v1,v2;
	do {
		v1=2.0*dRandUni(&(statePtr->unistate))-1.0; //from -1 to +1 in each direction, 
		v2=2.0*dRandUni(&(statePtr->unistate))-1.0; 
		rsq=v1*v1+v2*v2; // see if they are in the unit circle,
	} while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try again.
	fac=statePtr->sigma*sqrt(-2.0*log(rsq)/rsq);
	/* Now make the Box-Muller transformation to get two normal deviates. Return one and save the other for next time.	*/
	return DCplx(v2*fac+statePtr->mean, v1*fac+statePtr->mean);
/*
	SQRT( T1 )*EXP( CMPLX( ZERO, TWOPI*T2 ) )
	DCplx n;
	double sqrtv1 = sqrt(v1*fac);
	double angle = M_2_PI*v2*fac;
	n.re = sin(angle);
	n.im = cos(angle);
	n*=sqrtv1*statePtr->sigma+statePtr->mean;
	return n;
*/
}

int zbRandGaus(zRandGausStatePtr *statePtr, const ZVector &zv){
	
	if (zv.taille==0 || zv.vect==NULL) return 0;
	double fac,rsq,v1,v2;
	for(register int i=0;i<zv.taille;i++) {
		do {
			v1=2.0*dRandUni(&(statePtr->unistate))-1.0; //from -1 to +1 in each direction, 
			v2=2.0*dRandUni(&(statePtr->unistate))-1.0; 
			rsq=v1*v1+v2*v2; // see if they are in the unit circle,
		} while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try again.
		fac=statePtr->sigma*sqrt(-2.0*log(rsq)/rsq);
		zv.vect[i].re=v2*fac+statePtr->mean;
		zv.vect[i].im=v1*fac+statePtr->mean;
		//zv.vect[i]=zRandGaus(statePtr);
	}
	return 1;
}

const ZVector &zbRandGausApply(zRandGausStatePtr *statePtr, const ZVector &zv){
	if (zv.taille==0 || zv.vect==NULL)  throw CUtilitisException(__FILE__, __LINE__, EMPTY_VECTOR);
	double fac,rsq,v1,v2;
	for(register int i=0;i<zv.taille;i++) {
		do {
			v1=2.0*dRandUni(&(statePtr->unistate))-1.0; //from -1 to +1 in each direction, 
			v2=2.0*dRandUni(&(statePtr->unistate))-1.0; 
			rsq=v1*v1+v2*v2; // see if they are in the unit circle,
		} while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try again.
		fac=statePtr->sigma*sqrt(-2.0*log(rsq)/rsq);
		zv.vect[i].re+=v2*fac+statePtr->mean;
		zv.vect[i].im+=v1*fac+statePtr->mean;
		//zv.vect[i]=zRandGaus(statePtr);
	}
	return zv;
}


int zbRandGaus(zRandGausStatePtr *statePtr, DCplx *zv, int zv_size){
	
	if (zv_size==0 || zv==NULL) return 0;
	double fac,rsq,v1,v2;
	for(register int i=0;i<zv_size;i++) {
		do {
			v1=2.0*dRandUni(&(statePtr->unistate))-1.0; //from -1 to +1 in each direction, 
			v2=2.0*dRandUni(&(statePtr->unistate))-1.0; 
			rsq=v1*v1+v2*v2; // see if they are in the unit circle,
		} while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try again.
		fac=statePtr->sigma*sqrt(-2.0*log(rsq)/rsq);
		zv[i].re=v2*fac+statePtr->mean;
		zv[i].im=v1*fac+statePtr->mean;
	}
	return 1;
}

int dbRandGaus(dRandGausStatePtr *statePtr, double *dv, int dv_size) {
	if (dv_size==0 || dv==NULL) return 0;
	double fac,rsq,v1,v2;
	statePtr->iset=0;

	for(register int i=0;i<dv_size;i+=2) {
		do {
			v1=2.0*dRandUni(&(statePtr->unistate))-1.0; //from -1 to +1 in each direction, 
			v2=2.0*dRandUni(&(statePtr->unistate))-1.0; 
			rsq=v1*v1+v2*v2; // see if they are in the unit circle,
		} while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try again.
		fac=statePtr->sigma*sqrt(-2.0*log(rsq)/rsq);
		dv[i]=statePtr->sigma*v2*fac+statePtr->mean;
		if(i+1<dv_size) 
			dv[i+1]=v1*fac+statePtr->mean;
		else {
			statePtr->iset=1;
			statePtr->gset=v1*fac+statePtr->mean;
		}
	}
	return 1;
}

int RandBit1(TRandBitStatePtr *statePtr){
	unsigned int newbit;
	newbit = (statePtr->iseed >> 17) & 1 ^ (statePtr->iseed >> 4) &1 ^ (statePtr->iseed >> 4) & 1 ^ (statePtr->iseed &1);
	statePtr->iseed=(statePtr->iseed<<1)| newbit;
	return (int) newbit;
}

int RandBit2(TRandBitStatePtr  *statePtr){
	if (statePtr->iseed & IB18){
		statePtr->iseed=((statePtr->iseed ^MASK)<<1) | IB1;
		return 1;
	} else {
		statePtr->iseed <<=1;
		return 0;
	}
}

void bRandBit1(TRandBitStatePtr *statePtr, const WVector &dst) {
	if(dst.taille==0 ||dst.vect==NULL) return;
	for(register int i=0; i < dst.taille; i++) {
		dst.vect[i]=RandBit1(statePtr);
	}
}

void bRandBit2(TRandBitStatePtr *stateptr, const WVector &dst) {
	if(dst.taille==0 ||dst.vect==NULL) return;
	for(register int i=0; i < dst.taille; i++) {
		dst.vect[i]=RandBit2(stateptr);
	}
}

void bRandBit1(TRandBitStatePtr *statePtr, const WMatrix &dst) {
	if(dst.line<=0 ||dst.col<=0 || dst.mat==NULL) throw CUtilitisException(INVALID_MATRIX);
	for(register int i=0; i < dst.line; i++) {
		for(register int j=0; j<dst.col; j++) {
			dst.mat[i][j]=RandBit1(statePtr);
		}
	}
}

void bRandBit2(TRandBitStatePtr *stateptr, const WMatrix &dst) {
	if(dst.line<=0 ||dst.col<=0 || dst.mat==NULL) throw CUtilitisException(INVALID_MATRIX);
	for(register int i=0; i < dst.line; i++) {
		for(register int j=0; j<dst.col; j++) {
			dst.mat[i][j]=RandBit2(stateptr);
		}
	}
}


void bRandBit1NRZ(TRandBitStatePtr *statePtr, const WVector &dst) {
	if(dst.taille==0 ||dst.vect==NULL) return;
	for(register int i=0; i < dst.taille; i++) {
		dst.vect[i]=RandBit1(statePtr)==0?-1:1;
	}
}

void bRandBit2NRZ(TRandBitStatePtr *stateptr, const WVector &dst) {
	if(dst.taille==0 ||dst.vect==NULL) return;
	for(register int i=0; i < dst.taille; i++) {
		dst.vect[i]=RandBit2(stateptr)==0?-1:1;
	}
}

void bRandBit1NRZi(TRandBitStatePtr *statePtr, const WVector &dst) {
	if(dst.taille==0 ||dst.vect==NULL) return;
	for(register int i=0; i < dst.taille; i++) {
		dst.vect[i]=RandBit1(statePtr)==0?1:-1;
	}
}

void bRandBit2NRZi(TRandBitStatePtr *stateptr, const WVector &dst) {
	if(dst.taille==0 ||dst.vect==NULL) return;
	for(register int i=0; i < dst.taille; i++) {
		dst.vect[i]=RandBit2(stateptr)==0?1:-1;
	}
}

void bRandBit1NRZi(TRandBitStatePtr *statePtr, const WMatrix &dst) {
	if(dst.line<=0 ||dst.col<=0 || dst.mat==NULL) throw CUtilitisException(INVALID_MATRIX);
	for(register int i=0; i < dst.line; i++) {
		for(register int j=0; j<dst.col; j++) {
			dst.mat[i][j]=RandBit1(statePtr)?-1:1;
		}
	}
}

void bRandBit2NRZi(TRandBitStatePtr *stateptr, const WMatrix &dst) {
	if(dst.line<=0 ||dst.col<=0 || dst.mat==NULL) throw CUtilitisException(INVALID_MATRIX);
	for(register int i=0; i < dst.line; i++) {
		for(register int j=0; j<dst.col; j++) {
			dst.mat[i][j]=RandBit2(stateptr)?-1:1;
		}
	}
}

void bRandBit1NRZ(TRandBitStatePtr *statePtr, const WMatrix &dst) {
	if(dst.line<=0 ||dst.col<=0 || dst.mat==NULL) throw CUtilitisException(INVALID_MATRIX);
	for(register int i=0; i < dst.line; i++) {
		for(register int j=0; j<dst.col; j++) {
			dst.mat[i][j]=RandBit1(statePtr)?1:-1;
		}
	}
}

void bRandBit2NRZ(TRandBitStatePtr *stateptr, const WMatrix &dst) {
	if(dst.line<=0 ||dst.col<=0 || dst.mat==NULL) throw CUtilitisException(INVALID_MATRIX);
	for(register int i=0; i < dst.line; i++) {
		for(register int j=0; j<dst.col; j++) {
			dst.mat[i][j]=RandBit2(stateptr)?1:-1;
		}
	}
}

/*struct dRandUniStatePtr {
	double min, max, dist;
	int inext,inextp;
	int ma[56]; // The value 56 (range ma[1..55]) is special and should not be modied; see Knuth. 
	int iff;
	int idum;
	int  MSEED;
};*/
void dRandUniStatePtrCopy(dRandUniStatePtr &left, const dRandUniStatePtr &right) {
    left.min= right.min;
    left.max= right.max;
    left.dist= right.dist;
    left.inext= right.inext;
    left.inextp= right.inextp;
    for(register int i=0;i<56;i++) left.ma[i]=right.ma[i];
    left.iff= right.iff;
    left.idum= right.idum;
    left.MSEED= right.MSEED;
    //return left;
}

ostream& operator<< (ostream &os, dRandUniStatePtr &r){
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	os << r.min << " " << r.max << " " << r.dist << " " << r.inext << " " << r.inextp << endl;
	for(int i=0;i<56;i++) {
		os << r.ma[i] << " " ;
	}
	os << endl;
	os << r.iff << " " << r.idum << " " << r.MSEED << endl;
	return os;
}

istream& operator>> (istream &is, dRandUniStatePtr &r){
	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	is >> r.min  >> r.max  >> r.dist  >> r.inext  >> r.inextp ;
	for(int i=0;i<56;i++) {
		is >> r.ma[i];
	}
	is >> r.iff >> r.idum >> r.MSEED;
	return is;
}

/*struct dRandGausStatePtr {
	double mean, variance, sigma; 
	// int idum;
	dRandUniStatePtr unistate;
	int iset;
	double gset;
};*/

ostream& operator<< (ostream &os, dRandGausStatePtr &r){
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	os << r.mean << " " << r.variance << " " << r.sigma << endl;
	os << r.unistate;
	os << r.iset<<" "<<r.gset<< endl;
	return os;
}

istream& operator>> (istream &is, dRandGausStatePtr &r){
	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	is >> r.mean  >> r.variance >> r.sigma;
	is >> r.unistate;
	is >> r.iset>>r.gset;
	return is;
}

/*struct zRandGausStatePtr {
	double mean, variance, sigma; 
	// int idum;
	dRandUniStatePtr unistate;
};*/

ostream& operator<< (ostream &os, zRandGausStatePtr &r){
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	os << r.mean << " " << r.variance << " " << r.sigma << endl;
	os << r.unistate;
	return os;
}

istream& operator>> (istream &is, zRandGausStatePtr &r){
	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	is >> r.mean  >> r.variance >> r.sigma;
	is >> r.unistate;
	return is;
}
