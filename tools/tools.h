/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/07/30 14:18:31 $
 * $Revision: 1.1.2.16 $
 * $Id: tools.h,v 1.1.2.16 2006/07/30 14:18:31 syed Exp $
 ********************************************************************/

#ifndef _TOOLS_H_
#define _TOOLS_H_

#include <iostream>
using namespace std;
#include <math.h>


#include "tools/utilitis.h"

double log_2(double x);

int bin2oct(int *data, int data_left);
void oct2bin(int pskvalue, int *data_out);
int bin2quad(int *data, int data_left);
void quad2bin(int pskvalue, int *data_out);

double round(double d) ;

double calc_doppler(double freq, double vitesse, double angle =0);

void vector_increment(const WVector &v, int m_nary) ;

double pythag(double a, double b);

inline double nr_sign(double a, double b) {
	return (b >= 0.0) ? fabs(a) : -fabs(a);
}

template<class T> inline T nr_max(const T &a, const T &b) {
	return (a > b) ? a : b;
}

template<class T> inline T nr_min(const T &a, const T &b) {
	return (a < b) ? a : b;
}

template<class T> inline void utils_swap(T &a, T &b) {
	T temp=a;a=b;b=temp;
}

inline int sign(int val) {
	return val<0?-1:0;
}

inline double sign(double val) {
	return val<0.0?-1.0:0.0;
}

int gcd(int a, int b);

int isRelativePrime(int a, int b);
int genfilename(char *dest, int destlen, char *ext, int extlen, int seed=0);
#endif
