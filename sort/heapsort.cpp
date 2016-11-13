/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/04/05 18:49:18 $
 * $Revision: 1.1.2.5 $
 * $Id: heapsort.cpp,v 1.1.2.5 2004/04/05 18:49:18 syed Exp $
 ********************************************************************/
 
#include "tools/utilitis.h"

#include "sort/exception.h"

#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 

void hpsort(unsigned long n, double ra[]) { 
	/*
	Sorts an array ra[1..n] into ascending numerical order using the Heapsort algorithm. n is
	input; ra is replaced on output by its sorted rearrangement.
	*/
	unsigned long i,ir,j,l;
	double rra;
	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	/*The index l will be decremented from its initial value down to 1 during the \hiring" (heap
	creation) phase. Once it reaches 1, the index ir will be decremented from its initial value
	down to 1 during the \retirement-and-promotion" (heap selection) phase.*/
	for (;;) {
		if (l > 1) { //Still in hiring phase.
			rra=ra[--l];
		} else { //In retirement-and-promotion phase.
			rra=ra[ir]; //Clear a space at end of array.
			ra[ir]=ra[1]; //Retire the top of the heap into it.
			if (--ir == 1) { //Done with the last promotion.
				ra[1]=rra; //The least competent worker of all!
				break;
			}
		}
		i=l; //Whether in the hiring phase or promotion phase, we here set up to sift down element rra to its proper level.
		j=l+l;
		while (j <= ir) {
			if (j < ir && ra[j] < ra[j+1]) j++; //Compare to the better underling.

			if (rra < ra[j]) { //Demote rra.
				ra[i]=ra[j];
				i=j;
				j <<= 1;
			} else break; // Found rra's level. Terminate the sift-down.
		}
		ra[i]=rra; // Put rra into its slot.
	}
}

void hpsort_index(unsigned long n, double ra[], int rb[]) { 
	/*
	Sorts an array ra[1..n] into ascending numerical order using the Heapsort algorithm. n is
	input; ra is replaced on output by its sorted rearrangement.
	*/

	unsigned long i,ir,j,l;
	double rra; int rrb;
	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	/*The index l will be decremented from its initial value down to 1 during the \hiring" (heap
	creation) phase. Once it reaches 1, the index ir will be decremented from its initial value
	down to 1 during the \retirement-and-promotion" (heap selection) phase.*/

	for (;;) {
		if (l > 1) { //Still in hiring phase.
			rra=ra[--l];
			rrb=rb[l];
		} else { //In retirement-and-promotion phase.
			rra=ra[ir]; //Clear a space at end of array.
			rrb=rb[ir];
			ra[ir]=ra[1]; //Retire the top of the heap into it.
			rb[ir]=rb[1];
			if (--ir == 1) { //Done with the last promotion.
				ra[1]=rra; //The least competent worker of all!
				rb[1]=rrb;
				break;
			}
		}
		i=l; //Whether in the hiring phase or promotion phase, we here set up to sift down element rra to its proper level.
		j=l+l;
		while (j <= ir) {
			if (j < ir && ra[j] < ra[j+1]) j++; //Compare to the better underling.

			if (rra < ra[j]) { //Demote rra.
				ra[i]=ra[j];
				rb[i]=rb[j];
				i=j;
				j <<= 1;
			} else break; // Found rra's level. Terminate the sift-down.
		}
		ra[i]=rra; // Put rra into its slot.
		rb[j]=rrb;
	}
}

void heapsort(const DVector &ra) {

	double *rx = new double[ra.taille+1];
	
	rx[0]=0.0;
	register int r=0;
	for(r=1;r<=ra.taille;r++) rx[r]=ra.vect[r-1];
	hpsort(ra.taille, rx);
	for(r=1;r<=ra.taille;r++) ra.vect[r-1]=rx[r];
}

WVector heapsort_with_index(const DVector &ra) {
	WVector windex(ra.taille);
	int *index = new int[ra.taille+1];
	double *rx = new double[ra.taille+1];

	register int idx;

	for(idx=1;idx<=ra.taille;idx++){
		index[idx]=idx;
		rx[idx]=ra.vect[idx-1];
	}

	hpsort_index(ra.taille, rx, index);

	for(idx=1;idx<=ra.taille;idx++){
		windex.vect[idx-1]=index[idx];
		ra.vect[idx-1]=rx[idx];
	}
	return windex;
}



