/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/06/06 09:22:53 $
 * $Revision: 1.1.2.12 $
 * $Id: fft.cpp,v 1.1.2.12 2006/06/06 09:22:53 syed Exp $
 ********************************************************************/
 
#include <iostream>
using namespace std;

#include "tools/utilitis.h"
#include "fft/fft.h"
#include "tools/tools.h"
#ifdef WIN_DOS
#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 
#pragma intrinsic(atan, exp, log10, sqrt, atan2, log, sin, tan, cos, fabs, abs ,pow)
#endif

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
/*
This is a good place to remind you that you can also use a routine like four1
without modification even if your input data array is zero-offset, that is has the range
data[0..2*nn-1]. In this case, simply decrement the pointer to data by one when
four1 is invoked, e.g., four1(data-1,1024,1);. The real part of f0 will now be
returned in data[0], the imaginary part in data[1], and so on. See x1.2.
*/
/*
Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; or replaces
data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as -1.
data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn MUST
be an integer power of 2 (this is not checked for!).
*/
void four1(double *data, unsigned long nn, int isign) {
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta; // Double precision for the trigonometric recurrences. 
	double tempr,tempi;
	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) { //This is the bit-reversal section of the routine. 
		if (j > i) {
			SWAP(data[j],data[i]); // Exchange the two complex numbers.
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	/*Here begins the Danielson-Lanczos section of the routine. */
	mmax=2;
	
	while (n > mmax) { // Outer loop executed log2 nn times.
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax); // Initialize the trigonometric recurrence.
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) { // Here are the two nested inner loops.
			for (i=m;i<=n;i+=istep) {
				j=i+mmax; //This is the Danielson-Lanczos formula:
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr; //Trigonometric recurrence.
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

ZVector fft(const ZVector &td) {
	if (td.taille==0 ||td.vect==NULL) return ZVector();
	double x = log_2(td.taille);
	if (x!=ceil(x)) return ZVector(); // control taille= 2^N
	int limit=td.taille*2+1;
	double *dtd = new double[limit];
	
	for(register int i = 1, j=0; i < limit; j++, i+=2) {
		dtd[i]=td.vect[j].re;
		dtd[i+1]=td.vect[j].im;
	}
	DVector dbg;
	dbg.taille=limit-1;
	dbg.vect=dtd+1;
	//cout << "dbg1="<< dbg << endl;
	dbg.vect=NULL;
	dbg.taille=0;
	four1(dtd, td.taille, 1);
	dbg.taille=limit-1;
	dbg.vect=dtd+1;
	//cout << "dbg2="<< dbg << endl;
	dbg.vect=NULL;
	dbg.taille=0;

	ZVector fd(td.taille);
/*	for(register int f=1, g=0;f<limit;f+=2,g++) {
		fd.vect[g].re = dtd[f];
		fd.vect[g].im = dtd[f+1];
	}

	*/

	fd.vect[0].re=dtd[1];

	fd.vect[0].im=dtd[2];

	for(register int f=3, g=td.taille-1;f<limit;f+=2,g--) {

		fd.vect[g].re = dtd[f];

		fd.vect[g].im = dtd[f+1];

	}

#ifdef SHARED_NORMALIZATION

#ifdef WIN_DOS
#pragma message("INFORMATION : " __FILE__" ---->>> : using shared normalization")

#endif
	double normalize = 1.0/sqrt(double(td.taille));

	return fd*normalize;

#else 

#ifdef WIN_DOS
#pragma message("INFORMATION : " __FILE__" ---->>> : normalization is on ifft")

#endif
	return fd;

#endif


	delete[] dtd;
	return fd;
}

ZVector ifft(const ZVector &fd) {
	if (fd.taille==0 ||fd.vect==NULL) return ZVector();
	double x = log_2(fd.taille);
	if (x!=ceil(x)) return ZVector(); // control taille= 2^N
	int limit=fd.taille*2+1;
	double *dfd = new double[limit];
	
	for(register int i = 1, j=0; i < limit; j++, i+=2) {
		dfd[i]=fd.vect[j].re;
		dfd[i+1]=fd.vect[j].im;
	}
	four1(dfd, fd.taille, -1);
	ZVector td(fd.taille);
/*	for(register int f=1, g=0;f<limit;f+=2,g++) {
		td.vect[g].re = dfd[f];
		td.vect[g].im = dfd[f+1];
	}

	*/

	td.vect[0].re=dfd[1];

	td.vect[0].im=dfd[2];

	for(register int f=3, g=fd.taille-1;f<limit;f+=2,g--) {

		td.vect[g].re = dfd[f];

		td.vect[g].im = dfd[f+1];

	}

#ifdef SHARED_NORMALIZATION
	double normalize = 1.0/sqrt(double(fd.taille));
	return td*normalize;

#else 

	double normalize = 1.0/double(fd.taille);

	return td*normalize;

#endif
}

ZVector zero_center(const ZVector &ffted) {
	if (ffted.taille==0 || ffted.vect==NULL) return ZVector();
	ZVector zc(ffted.taille);
	int mid = ffted.taille >> 1;
	register int i, j;

	for(i = mid-1, j=ffted.taille-1; i >= 0; i--,j--) {
		zc.vect[i]=ffted.vect[j];
	}
	for(i=mid,j=0;i<ffted.taille;i++,j++) {
		zc.vect[i]=ffted.vect[j];
	}
	return zc;
}

ZVector butterfly(const ZVector &zero_centered) {
	if (zero_centered.taille==0 || zero_centered.vect==NULL) return ZVector();
	ZVector b(zero_centered.taille);
	int mid= zero_centered.taille>>1;
	register int i, j;
	for(i=mid-1,j=zero_centered.taille-1;i>=0;i--,j--) {
		b.vect[j]=zero_centered.vect[i];
	}
	for(i=mid,j=0;i<zero_centered.taille;i++,j++) {
		b.vect[j]=zero_centered.vect[i];
	}
	return b;
}

