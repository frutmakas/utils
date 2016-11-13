/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/04/05 18:49:18 $
 * $Revision: 1.1.2.4 $
 * $Id: quicksort.cpp,v 1.1.2.4 2004/04/05 18:49:18 syed Exp $
 ********************************************************************/
 
#include "sort/quicksort.h"
#include "sort/exception.h"

//#include "nr/nrutil.h"
#include "tools/tools.h"
#ifdef WIN_DOS
#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 
#endif
#define M 7
#define NSTACK 500
/*Here M is the size of subarrays sorted by straight insertion and NSTACK is the required auxiliary
storage.*/
//void sort(unsigned long n, double arr[]){
void quicksort(const DVector &arr) {
	/*
		Sorts an array arr.vect[1..n] into ascending numerical order using the Quicksort algorithm. n is
		input; arr is replaced on output by its sorted rearrangement.
	*/
	unsigned long i,ir=arr.taille-1/*n*/,j,k,l=0;//1,*istack;
	int jstack=0;
	double a;//,temp;
	//istack=lvector(1,NSTACK);
	WVector istack=WVector(NSTACK+1);
	for (;;) { //Insertion sort when subarray small enough.
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr.vect[j];
				for (i=j-1;i>=l;i--) {
					if (arr.vect[i] <= a) break;
					arr.vect[i+1]=arr.vect[i];
				}
				arr.vect[i+1]=a;
			}
			if (jstack == 0) break;
			ir=istack.vect[jstack--]; //Pop stack and begin a new round of partitioning.
			l=istack.vect[jstack--];
		} else {
			k=(l+ir) >> 1; 
			//Choose median of left, center, and right elements as partitioning element a. Also rearrange so that a[l] . a[l+1] . a[ir].
			swap(arr.vect[k],arr.vect[l+1]);
			if (arr.vect[l] > arr.vect[ir]) {
				swap(arr.vect[l],arr.vect[ir]);
			}
			if (arr.vect[l+1] > arr.vect[ir]) {
				swap(arr.vect[l+1],arr.vect[ir]);
			}
			if (arr.vect[l] > arr.vect[l+1]) {
				swap(arr.vect[l],arr.vect[l+1]);
			}
			i=l+1; //Initialize pointers for partitioning.
			j=ir;
			a=arr.vect[l+1]; //Partitioning element.
			for (;;) { //Beginning of innermost loop.
				do i++; while (arr.vect[i] < a); //Scan up to nd element > a.
				do j--; while (arr.vect[j] > a); //Scan down to nd element < a.
				if (j < i) break; //Pointers crossed. Partitioning complete.
				swap(arr.vect[i],arr.vect[j]); //Exchange elements.
			} // End of innermost loop.
			arr.vect[l+1]=arr.vect[j]; //Insert partitioning element.
			arr.vect[j]=a;
			jstack += 2;
			//Push pointers to larger subarray on stack, process smaller subarray immediately.
			if (jstack > NSTACK) {
				throw CSortException("NSTACK too small in sort.", TYPE_QUICKSORT);
			}
			if (ir-i+1 >= j-l) {
				istack.vect[jstack]=ir;
				istack.vect[jstack-1]=i;
				ir=j-1;
			} else {
				istack.vect[jstack]=j-1;
				istack.vect[jstack-1]=l;
				l=i;
			}
		}
	}
//	free_lvector(istack,1,NSTACK);
}


//void sort2(unsigned long n, double arr[], double brr[]) {
WVector quicksort_with_index(const DVector &arr) {
/*Sorts an array arr.vect[1..n] into ascending order using Quicksort, while making the corresponding
	rearrangement of the array brr.vect[1..n].*/

	WVector brr(arr.taille);
	for(int idx=0;idx<arr.taille;idx++) 
		brr.vect[idx]=idx;
	unsigned long i,ir=arr.taille-1/*n*/,j,k,l=0;//1,*istack;
	int jstack=0;
	double a;//,temp;
	int b;
	
	WVector istack(NSTACK+1);
	for (;;) { //Insertion sort when subarray small enough.
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr.vect[j];
				b=brr.vect[j];
				for (i=j-1;i>=l;i--) {
					if (arr.vect[i] <= a) break;
					arr.vect[i+1]=arr.vect[i];
					brr.vect[i+1]=brr.vect[i];
				}
				arr.vect[i+1]=a;
				brr.vect[i+1]=b;
			}
			if (!jstack) {
				//free_lvector(istack,1,NSTACK);
				return brr;
			}
			ir=istack.vect[jstack]; // Pop stack and begin a new round of partitioning.
			l=istack.vect[jstack-1];
			jstack -= 2;
		} else {
			k=(l+ir) >> 1; //Choose median of left, center and right elements as partitioning element a. Also rearrange so that a[l] . a[l+1] . a[ir].
			swap(arr.vect[k],arr.vect[l+1]);
			swap(brr.vect[k],brr.vect[l+1]);
			if (arr.vect[l] > arr.vect[ir]) {
				swap(arr.vect[l],arr.vect[ir]);
				swap(brr.vect[l],brr.vect[ir]);
			}
			if (arr.vect[l+1] > arr.vect[ir]) {
				swap(arr.vect[l+1],arr.vect[ir]);
				swap(brr.vect[l+1],brr.vect[ir]);
			}
			if (arr.vect[l] > arr.vect[l+1]) {
				swap(arr.vect[l],arr.vect[l+1]);
				swap(brr.vect[l],brr.vect[l+1]);
			}
			i=l+1; // Initialize pointers for partitioning.
			j=ir;
			a=arr.vect[l+1]; // Partitioning element.
			b=brr.vect[l+1];
			for (;;) { // Beginning of innermost loop.
				do i++; while (arr.vect[i] < a); //Scan up to nd element > a.
				do j--; while (arr.vect[j] > a); //Scan down to nd element < a.
				if (j < i) break; // Pointers crossed. Partitioning complete.
				swap(arr.vect[i],arr.vect[j]); // Exchange elements of both arrays.
				swap(brr.vect[i],brr.vect[j]);
			}// End of innermost loop.
			arr.vect[l+1]=arr.vect[j]; // Insert partitioning element in both arrays.
			arr.vect[j]=a;
			brr.vect[l+1]=brr.vect[j];
			brr.vect[j]=b;
			jstack += 2;
			//Push pointers to larger subarray on stack, process smaller subarray immediately.
			if (jstack > NSTACK) throw CSortException("NSTACK too small in sort.", TYPE_QUICKSORT);
			if (ir-i+1 >= j-l) {
				istack.vect[jstack]=ir;
				istack.vect[jstack-1]=i;
				ir=j-1;
			} else {
				istack.vect[jstack]=j-1;
				istack.vect[jstack-1]=l;
				l=i;
			}
		}
	}
}
#ifdef M
#undef M 
#endif

