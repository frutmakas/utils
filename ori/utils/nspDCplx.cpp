#ifndef _NSPDCPLX_H_

#include "globaldef.h"

#ifdef _NSP_INTEL_

#define nsp_UsesVector
#define nsp_UsesConversion

#include <nsp.h>
#include <math.h>

DCplx operator+(const DCplx &c1, const DCplx &c2){
	DCplx r = {c1.re+c2.re, c1.im+c2.im};
	return r;
}

DCplx operator+(const DCplx &c1, const double &c2){
	DCplx r = {c1.re+c2, c1.im};
	return r;
}

DCplx operator+(const double &c2, const DCplx &c1){
	DCplx r = {c1.re+c2, c1.im};
	return r;
}

DCplx operator-(const DCplx &c1, const DCplx &c2){
	DCplx r = {c1.re-c2.re, c1.im-c2.im};
	return r;
}

DCplx operator-(const DCplx &c1, const double &c2){
	DCplx r = {c1.re-c2, c1.im};
	return r;
}

DCplx operator-(const double &c2, const DCplx &c1){
	DCplx r = {c2-c1.re, c1.im};
	return r;
}

DCplx operator*(const DCplx &c1, const DCplx &c2){
	return nspzMpy(c1,c2);
}

DCplx operator*(const DCplx &c1, const double &c2){
	DCplx r = {c2, 0.0};
	return nspzMpy(c1,r);
}

DCplx operator*(const double &c2, const DCplx &c1){
	DCplx r = {c2, 0.0};
	return nspzMpy(c1,r);
}

DCplx operator/(const DCplx &c1, const DCplx &c2){
	return nspzDiv(c1,c2);
}

DCplx operator/(const DCplx &c1, const double &c2){
	DCplx r = {c2, 0.0};
	return nspzDiv(c1,r);
}

DCplx operator/(const double &c2, const DCplx &c1){
	DCplx r = {c2, 0.0};
	return nspzDiv(r,c1);
}

DCplx conj(const DCplx &c){
	DCplx r={c.re, -c.im};
	return r;
}

double dcplxmag(DCplx &c) {
	return sqrt(c.re*c.re+c.im*c.im);
}

double dcplxmod2(DCplx &c) {
	return c.re*c.re+c.im*c.im;
}

#endif
#endif