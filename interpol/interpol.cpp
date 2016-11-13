/* 
	NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0- 43108-5)
	Copyright (C) 1988-1992 by Cambridge University Press.
	Programs Copyright (C) 1988-1992 by Numerical Recipes Software.
	Permission is granted for internet users to make one paper copy for their
	own personal use. Further reproduction, or any copying of machinereadable
	files (including this one) to any server computer, is strictly prohibited.
	To order Numerical Recipes books, diskettes, or CDROMs visit website 
	http://www.nr.com or call 1-800-872-7423 (North America only),
	or send email to trade@cup.cam.ac.uk (outside North America).
*/

#include <stdlib.h>
#include <math.h>
#include "nr/nrutil.h"
#include "interpol/interpol.h"
#ifdef WIN_DOS
#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 
#endif
/*
	Here is a routine for polynomial interpolation or extrapolation from N input
	points. Note that the input arrays are assumed to be unit-offset. If you have
	zero-offset arrays, remember to subtract 1
  
	Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns 
	a value y, and an error estimate dy. If P(x) is the polynomial of degree N - 1 
	such that P(xai) = yai; i = 1; : : : ;n, then the returned value y = P(x).
*/

void polint(double xa[], double ya[], int n, double x, double *y, double *dy) {

	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;
	dif=fabs(x-xa[1]);
	c=nr_vector(1,n);
	d=nr_vector(1,n);
	for (i=1;i<=n;i++) { //Here we find the index ns of the closest table entry,
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i]; //and initialize the tableau of c's and d's.
		d[i]=ya[i];
	}

	*y=ya[ns--]; //This is the initial approximation to y.
	
	for (m=1;m<n;m++) { //For each column of the tableau,
		for (i=1;i<=n-m;i++) { //we loop over the current c's and d's and update them.
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
			//This error can occur only if two input xa's are (to within roundoff) identical.
			den=w/den;
			d[i]=hp*den; ///Here the c's and d's are updated.
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
		/*  After each column in the tableau is completed, we decide which correction, c or d,
			we want to add to our accumulating value of y, i.e., which path to take through the
			tableau|forking up or down. We do this in such a way as to take the most \straight
			line" route through the tableau to its apex, updating ns accordingly to keep track of
			where we are. This route keeps the partial approximations centered (insofar as possible)
			on the target x. The last dy added is thus the error indication.
		*/
	}
	nr_free_vector(d,1,n);
	nr_free_vector(c,1,n);
}

/* 
	Rational Function Interpolation and Extrapolation 
	
	Given arrays xa[1..n] and ya[1..n], and given a value of x, this routine returns a value of
	y and an accuracy estimate dy. The value returned is that of the diagonal rational function,
	evaluated at x, which passes through the n points (xai; yai), i = 1:::n.

*/

#define TINY 1.0e-25 //A small number.
//#define FREERETURN {free_vector(d,1,n);free_vector(c,1,n);return;}

void ratint(double xa[], double ya[], int n, double x, double *y, double *dy) {
	int m,i,ns=1;
	double w,t,hh,h,dd,*c,*d;

	c=nr_vector(1,n);
	d=nr_vector(1,n);
	hh=fabs(x-xa[1]);
	for (i=1;i<=n;i++) {
		h=fabs(x-xa[i]);
		if (h == 0.0) {
			*y=ya[i];
			*dy=0.0;
			nr_free_vector(d,1,n);
			nr_free_vector(c,1,n);
			return;
		} else if (h < hh) {
			ns=i;
			hh=h;
		}
		c[i]=ya[i];
		d[i]=ya[i]+TINY; //The TINY part is needed to prevent a rare zero-over-zero condition. 
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			w=c[i+1]-d[i];
			h=xa[i+m]-x; // h will never be zero, since this was tested in the initializing loop. 
			t=(xa[i]-x)*d[i]/h;
			dd=t-c[i+1];
			if (dd == 0.0) nrerror("Error in routine ratint");
			/* This error condition indicates that the interpolating function 
			has a pole at the requested value of x. */
			dd=w/dd;
			d[i]=c[i+1]*dd;
			c[i]=t*dd;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	nr_free_vector(d,1,n);
	nr_free_vector(c,1,n);
	return;
}

/*
	Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., 
	yi = f(xi), with x1 <x2 < :: : < xN, and given values yp1 and ypn for the
	first derivative of the interpolating function at points 1 and n, 
	respectively, this routine returns an array y2[1..n] that contains the
	second derivatives of the interpolating function at the tabulated points xi. 
	If yp1 and/or ypn are equal to 1 . 1030 or larger, the routine is signaled 
	to set the corresponding boundary 	condition for a natural spline, with 
	zero second derivative on that boundary.
*/
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]){
	int i,k;
	double p,qn,sig,un,*u;
	u=nr_vector(1,n-1);
	if (yp1 > 0.99e30) //The lower boundary condition is set either to be \natural"
		y2[1]=u[1]=0.0;
	else { //or else to have a specified first derivative.
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	/*
		This is the decomposition loop of the tridiagonal algorithm.
		y2 and u are used for temporary	storage of the decomposed factors.
	*/
	for (i=2;i<=n-1;i++) { 
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30) //The upper boundary condition is set either to be \natural" 
		qn=un=0.0;
	else { // or else to have a specified first derivative.
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--) // This is the backsubstitution loop of the tridiagonal algorithm. 
		y2[k]=y2[k]*y2[k+1]+u[k];
	nr_free_vector(u,1,n-1);
}


/*
	Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with 
	the xai's in order), and given the array y2a[1..n], which is the output 
	from spline above, and given a value of x, this routine returns a 
	cubic-spline interpolated value y.
*/

void splint(double xa[], double ya[], double y2a[], int n, double x, double *y){
	void nrerror(char error_text[]);
	int klo,khi,k;
	double h,b,a;

	/*
		We will find the right place in the table by means of bisection. This 
		is optimal if sequential calls to this routine are at random values 
		of x. If sequential calls are in order, and closely spaced, one would 
		do better to store previous values of klo and khi and test if they 
		remain appropriate on the next call.
	*/
	khi=n;
	klo=1; 
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	} //klo and khi now bracket the input value of x.
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad xa input to routine splint"); //The xa's must be distinct.
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h; //Cubic spline polynomial is now evaluated.
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}


