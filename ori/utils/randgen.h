#ifndef _RAND_GEN_H_
#define _RAND_GEN_H_
#include "globaldef.h"

#ifndef _NSP_INTEL_
#include "utilitis.h"

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

struct dRandGausStatePtr {
	double mean, variance, sigma; 
	// int idum;
	dRandUniStatePtr unistate;
	int iset;
	double gset;
};

struct zRandGausStatePtr {
	double mean, variance, sigma; 
	// int idum;
	dRandUniStatePtr unistate;
};



void dRandUniInit(dRandUniStatePtr *statePtr, int seed=0, double min=0.0, double max=1.0);
double dRandUni(dRandUniStatePtr *statePtr);

void dRandGausInit(dRandGausStatePtr *statePtr, double mean=0.0, double variance=1.0);
double dRandGaus(dRandGausStatePtr *statePtr) ;

void zRandGausInit(zRandGausStatePtr *statePtr, double mean=0.0, double variance=1.0);
DCplx zRandGaus(zRandGausStatePtr *statePtr) ;

int dbRandGaus(dRandGausStatePtr *statePtr, const DVector &dv);
int zbRandGaus(zRandGausStatePtr *statePtr, double *dv, int dv_size);
int zbRandGaus(zRandGausStatePtr *statePtr, const ZVector &zv);
int zbRandGaus(zRandGausStatePtr *statePtr, DCplx *zv, int zv_size);



#define IB1 1
#define IB2 2
#define IB5 16
#define IB18 131072
#define MASK (IB1+IB2+IB5)

struct TRandBitStatePtr {
	unsigned long iseed;
};

int RandBit1(TRandBitStatePtr *statePtr);
int RandBit2(TRandBitStatePtr *statePtr);

void bRandBit1(TRandBitStatePtr *statePtr, const WVector &dst);
void bRandBit2(TRandBitStatePtr *statePtr, const WVector &dst);


#endif //_NSP_INTEL_
#endif