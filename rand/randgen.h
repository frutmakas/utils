/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/07/21 10:58:00 $
 * $Revision: 1.1.2.13 $
 * $Id: randgen.h,v 1.1.2.13 2006/07/21 10:58:00 syed Exp $
 ********************************************************************/
 
#ifndef _RAND_GEN_H_
#define _RAND_GEN_H_
#include <iostream>
using namespace std;
#include "globaldef.h"

#include "tools/utilitis.h"

#define MBIG 1000000000
#define MSEED_MAX 161803398
#define MZ 0
#define FAC (1.0/MBIG)
// According to Knuth, any large MBIG, and any smaller (but still large) MSEED can be substituted for the above values.


struct dRandUniStatePtr {
	double min, max, dist;
	int inext,inextp;
	int ma[56]; // The value 56 (range ma[1..55]) is special and should not be modied; see Knuth. 
	int iff;
	int idum;
	int  MSEED;
};

ostream& operator<< (ostream &os, dRandUniStatePtr &r);
istream& operator>> (istream &is, dRandUniStatePtr &r);
//dRandUniStatePtr &operator=(dRandUniStatePtr &left, const dRandUniStatePtr &right);
void dRandUniStatePtrCopy(dRandUniStatePtr &left, const dRandUniStatePtr &right);

struct dRandGausStatePtr {
	double mean, variance, sigma; 
	// int idum;
	dRandUniStatePtr unistate;
	int iset;
	double gset;
};

ostream& operator<< (ostream &os, dRandGausStatePtr &r);
istream& operator>> (istream &is, dRandGausStatePtr &r);

struct zRandGausStatePtr {
	double mean, variance, sigma; 
	// int idum;
	dRandUniStatePtr unistate;
};

ostream& operator<< (ostream &os, zRandGausStatePtr &r);
istream& operator>> (istream &is, zRandGausStatePtr &r);



void dRandUniInit(dRandUniStatePtr *statePtr, int seed=0, double min=0.0, double max=1.0);
double dRandUni(dRandUniStatePtr *statePtr);
void dbRandUni(dRandUniStatePtr *statePtr, const DVector &dv);


void dRandGausInit(dRandGausStatePtr *statePtr, double mean=0.0, double variance=1.0);
void dRandGausInit(dRandGausStatePtr *statePtr, int seed, double mean, double variance);
double dRandGaus(dRandGausStatePtr *statePtr) ;

void zRandGausInit(zRandGausStatePtr *statePtr, double mean=0.0, double variance=1.0);
void zRandGausInit(zRandGausStatePtr *statePtr, double mean, double variance, int seed);

DCplx zRandGaus(zRandGausStatePtr *statePtr) ;

int dbRandGaus(dRandGausStatePtr *statePtr, const DVector &dv);
int zbRandGaus(zRandGausStatePtr *statePtr, double *dv, int dv_size);
int zbRandGaus(zRandGausStatePtr *statePtr, const ZVector &zv);
int zbRandGaus(zRandGausStatePtr *statePtr, DCplx *zv, int zv_size);

const ZVector& zbRandGausApply(zRandGausStatePtr *statePtr, const ZVector &zv);
int dbApplyRandGaus(dRandGausStatePtr *statePtr, const DVector &dv) ;

int dmRandGaus(dRandGausStatePtr *statePtr, const DMatrix &dm);
int zmRandGaus(zRandGausStatePtr *statePtr, const ZMatrix &dm);


#define IB1 1
#define IB2 2
#define IB5 16
#define IB18 131072
#define MASK (IB1+IB2+IB5)

struct TRandBitStatePtr {
	unsigned long iseed;
};

typedef TRandBitStatePtr wRandBitStatePtr;

int RandBit1(TRandBitStatePtr *statePtr);
int RandBit2(TRandBitStatePtr *statePtr);

void bRandBit1(TRandBitStatePtr *statePtr, const WVector &dst);
void bRandBit2(TRandBitStatePtr *statePtr, const WVector &dst);
void bRandBit1(TRandBitStatePtr *stateptr, const WVector &dst);
void bRandBit2(TRandBitStatePtr *stateptr, const WVector &dst);
void bRandBit1NRZi(TRandBitStatePtr *statePtr, const WMatrix &dst);
void bRandBit2NRZi(TRandBitStatePtr *stateptr, const WMatrix &dst);
void bRandBit1NRZ(TRandBitStatePtr *statePtr, const WMatrix &dst);
void bRandBit2NRZ(TRandBitStatePtr *stateptr, const WMatrix &dst);
void bRandBit1NRZ(TRandBitStatePtr *statePtr, const WVector &dst);
void bRandBit2NRZ(TRandBitStatePtr *stateptr, const WVector &dst);
void bRandBit1NRZi(TRandBitStatePtr *statePtr, const WVector &dst);
void bRandBit2NRZi(TRandBitStatePtr *stateptr, const WVector &dst);

#endif
