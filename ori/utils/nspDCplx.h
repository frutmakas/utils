#ifndef _NSPDCPLX_H_

#ifdef _NSP_INTEL_
#define nsp_UsesVector
#include <nsp.h>


DCplx operator+(const DCplx &c1, const DCplx &c2);

DCplx operator+(const DCplx &c1, const double &c2);

DCplx operator+(const double &c2, const DCplx &c1);

DCplx operator-(const DCplx &c1, const DCplx &c2);

DCplx operator-(const DCplx &c1, const double &c2);

DCplx operator-(const double &c2, const DCplx &c1);

DCplx operator*(const DCplx &c1, const DCplx &c2);
DCplx operator*(const DCplx &c1, const double &c2);

DCplx operator*(const double &c2, const DCplx &c1);

DCplx operator/(const DCplx &c1, const DCplx &c2);
DCplx operator/(const DCplx &c1, const double &c2);

DCplx operator/(const double &c2, const DCplx &c1);

DCplx conj(const DCplx &c);
double dcplxmag(DCplx &c) ;
double dcplxmod2(DCplx &c) ;

#endif
#endif