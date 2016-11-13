/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/07/28 11:03:32 $
 * $Revision: 1.1.2.16 $
 * $Id: tools.cpp,v 1.1.2.16 2006/07/28 11:03:32 syed Exp $
 ********************************************************************/

#include "tools/tools.h"
#ifdef WIN_DOS
#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 
#pragma intrinsic(log)
#endif
double log_2(double x) {
	if (x<=0) return 0; 
	return log(x)/log(2.0);
}

int bin2oct(int *data, int data_left) {
	int x = data_left>=3 ?3 : data_left;
	int oct=0;
	for(int i=0;i<3;i++) {
		if (x>i) oct = (oct << 1) | data[i];
		else oct <<= 1; 
	}
	return oct;
}

void oct2bin(int pskvalue, int *data_out) {
	data_out[0]=(pskvalue & 4)==4 ? 1 : 0;
	data_out[1]=(pskvalue & 2)==2 ? 1:0;
	data_out[2]=(pskvalue & 1)==1 ? 1:0;
}


int bin2quad(int *data, int data_left) {
	int x = data_left>=2 ?2 : data_left;
	int quad=0;
	for(int i=0;i<2;i++) {
		if (x>i) quad = (quad<< 1) | data[i];
		else quad <<= 1; 
	}
	return quad;
}

void quad2bin(int pskvalue, int *data_out) {
	data_out[0]=(pskvalue & 2)==2 ? 1:0;
	data_out[1]=(pskvalue & 1)==1 ? 1:0;
}

double round(double d) {
	double t = floor(d);
	return ((d-t)>=0.5) ? ceil(d) : t;
}

void vector_increment(const WVector &v, int m_nary)  {

	v.vect[0]++;

	for(int i=0;i<v.taille;i++) {

		if(v.vect[i]==m_nary) {

			if (i<v.taille-1) v.vect[i+1]++;

			v.vect[i]=0;

		} else { break; }

	}

}

double calc_doppler(double freq /* Hz */, double vitesse /* km/h */, double angle /* radian */) {
	double vitesse_ms = vitesse*1000.0/3600.0;
	double vitesse_c = 299792458.0;
	return (vitesse_ms*freq)/vitesse_c*cos(angle);
}

//Computes (a2 + b2)1=2 without destructive underflow or overflow.
double pythag(double a, double b) {
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	//if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	if (absa > absb) {
		if(absb == 0) return absa;
		double abss = absb/absa;
		return absa*sqrt(1.0+abss*abss);
	}
	else {
		double abss = absa/absb;
		return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+abss*abss));
	}
}

int gcd(int a, int b) {
    if(a==b) return a;
    if(a<=0 || b<=0) return -1;

    for(int i=a<b ? a : b;i>1;i--)
        if(b%i==0 && a%i==0) return i;
    return 1;
}

int isRelativePrime(int a, int b) {
    if(a==b) return 0;
    if(a<=0 || b<=0)
        return -1;
    if ( (a&1)==0 && (b&1)==0)
        return 0; // if both multiple of 2 .. not relatively prime ..
    if(gcd(a,b)!=1)
        return 0;
    else
        return 1;
}

int genfilename(char *dest, int destlen, char *ext, int extlen, int seed) {
    if(dest==NULL || destlen <=0 || extlen>=destlen) return 0;
    int t;
    if(seed == 0) seed = time(NULL);
    int i=0;
    for(i=0;i<destlen-extlen-1;i++) {
        do {
            t= rand()%80+48;
        } while(t<48 || t >122 ||(t>57 && t<65) || (t>90 && t<97));
        dest[i]= t;
    }
    if(ext!=NULL && extlen>0) {
        for(int j=0;j<extlen;j++, i++) {
            dest[i]= ext[j];
        }
    }
    dest[destlen-1]=0;
    return 1;
}
