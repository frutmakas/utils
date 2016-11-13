#include "globaldef.h"

#ifndef _NSP_INTEL_

#include <math.h> 
#include "randgen.h"



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
	statePtr->MSEED=MSEED_MAX-(seed%1024);
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


void dRandGausInit(dRandGausStatePtr *statePtr, double mean, double variance){
	statePtr->unistate.idum=-(int)ceil(fabs((mean-variance)*(mean+variance)));
	statePtr->iset=0;
	statePtr->mean=mean;
	statePtr->variance=variance;
	statePtr->sigma=variance*variance;

	dRandUniInit(&(statePtr->unistate));
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

int dbRandGaus(dRandGausStatePtr *statePtr, const DVector &dv) {
	if (dv.taille==0 || dv.vect==NULL) return 0;
	double fac,rsq,v1,v2;
	statePtr->iset=0;

	for(register int i=0;i<dv.taille;i+=2) {
		do {
			v1=2.0*dRandUni(&(statePtr->unistate))-1.0; //from -1 to +1 in each direction, 
			v2=2.0*dRandUni(&(statePtr->unistate))-1.0; 
			rsq=v1*v1+v2*v2; // see if they are in the unit circle,
		} while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try again.
		fac=statePtr->sigma*sqrt(-2.0*log(rsq)/rsq);
		dv.vect[i]=statePtr->sigma*v2*fac+statePtr->mean;
		if(i+1<dv.taille) 
			dv.vect[i+1]=v1*fac+statePtr->mean;
		else {
			statePtr->iset=1;
			statePtr->gset=v1*fac+statePtr->mean;
		}
	}
	return 1;
}

void zRandGausInit(zRandGausStatePtr *statePtr, double mean, double variance){
	statePtr->unistate.idum=-3;
	statePtr->mean=mean;
	statePtr->variance=fabs(variance);
	statePtr->sigma=variance*variance;

	dRandUniInit(&(statePtr->unistate));
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
	fac=sqrt(-2.0*log(rsq)/rsq);
	/* Now make the Box-Muller transformation to get two normal deviates. Return one and save the other for next time.	*/
	return DCplx(statePtr->sigma*v2*fac+statePtr->mean, statePtr->sigma*v1*fac+statePtr->mean);
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
	}
	return 1;
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
		dst.vect[i]=RandBit1(stateptr);
	}
}

#endif // _NSP_INTEL_
