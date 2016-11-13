#ifndef _NSPMATRICE_H_
#define _NSPMATRICE_H_

#include "../globaldef.h"
#ifdef _NSP_INTEL_

#define nsp_UsesVector
#define nsp_UsesConversion

#include <nsp.h>
#include <stdlib.h>

class DMatrix;

class ZMatrix  {
	public :
		int col, line;
		
		DCplx **mat;

		ZMatrix (int col, int line);
		ZMatrix (int col, int line, DCplx x);
		ZMatrix(const ZMatrix& mat);
		ZMatrix(const DMatrix& mat);
		ZMatrix() { col=line=0;mat=NULL;}
		~ZMatrix();
		//void Free();
		void transpose();
		ZMatrix  transpose_nip();

		ZMatrix  operator+(const ZMatrix &rhs);
		ZMatrix  operator+(const DMatrix &rhs);
		ZMatrix  operator-(const ZMatrix &rhs);
		ZMatrix  operator-(const DMatrix &rhs);

		ZMatrix operator*(const ZMatrix &rhs);
		ZMatrix operator*(const DMatrix &rhs);
		ZMatrix operator*(const double &rhs);
		ZMatrix operator*(const DCplx &rhs);
		

		ZMatrix operator/(const double &rhs);
		ZMatrix operator/(const DCplx &rhs);


		ZMatrix &operator=(const ZMatrix &rhs); 


		ZMatrix &operator=(const DMatrix &rhs);

		bool operator==(const ZMatrix &rhs);

		bool operator!=(const ZMatrix &rhs);

		DMatrix real();
		DMatrix imag();
		DMatrix abs();
		DMatrix arg();
		ZMatrix conj();
		/* on fera ca quand on aura besoin ... 
			DMatrix norm(const ZMatrix x);
			ZMatrix polar(const T rho, const T theta = 0);
			ZMatrix cos(const ZMatrix x);
			ZMatrix cosh(const ZMatrix x);
			ZMatrix exp(const ZMatrix x);
			ZMatrix log(const ZMatrix x);
			ZMatrix log10(const ZMatrix x);
			ZMatrix pow(const ZMatrix x, int y);
			ZMatrix pow(const ZMatrix x, const double y);
			ZMatrix sin(const ZMatrix x);
			ZMatrix sinh(const ZMatrix x);
			ZMatrix sqrt(const ZMatrix x);
		*/
};

class DMatrix  {
	public :
		int col, line;
		double **mat;
		DMatrix(int col, int line, double initial=0.0);
		DMatrix(const DMatrix& m);
		DMatrix() { col=line=0;mat=NULL; }
		//void Free();
		~DMatrix();

		void transpose();
		DMatrix transpose_nip();
		ZMatrix ztranspose_nip();

		DMatrix operator+(const DMatrix &rhs);
		DMatrix operator-(const DMatrix &rhs);
		DMatrix operator*(const DMatrix &rhs);

		ZMatrix operator+(const ZMatrix &rhs);
		ZMatrix operator-(const ZMatrix &rhs);
		ZMatrix operator*(const ZMatrix &rhs);
		
		DMatrix operator*(const double &rhs);
		ZMatrix operator*(const DCplx &rhs);
		

		DMatrix operator/(const double &rhs);
		ZMatrix operator/(const DCplx &rhs);

		DMatrix &operator=(const DMatrix &rhs);

		bool operator==(const DMatrix &rhs);

		bool operator!=(const DMatrix &rhs);
};

ZMatrix operator*(const DCplx &lhs, const ZMatrix &rhs);
ZMatrix operator*(const double &lhs, const ZMatrix &rhs);
DMatrix operator*(const double &lhs, const DMatrix &rhs);
ZMatrix operator*(const DCplx &lhs, const DMatrix &rhs);

#endif //_NSP_INTEL_

#endif