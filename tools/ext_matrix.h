/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/07/28 11:03:32 $
 * $Revision: 1.1.2.6 $
 * $Id: ext_matrix.h,v 1.1.2.6 2006/07/28 11:03:32 syed Exp $
 ********************************************************************/
#ifndef _EXT_MATRIX_H_
#define _EXT_MATRIX_H_

#include "tools/all.h"
//#include "clapack.h"

#define NO_CONVERGENCE 9979
#define UNHANDLED_EXCEPTION 9978

class CLapackSVDStatus {
	public: 
		long info;
		doublereal *rwork;
		long lwork;
		doublecomplex *work;
		long ldvt;
		doublecomplex *vt;
		long ldu;
		doublecomplex *u;
		doublereal *s;
		long lda;
		doublecomplex * a;
		long m,n;
		char jobu,jobvt;
		int worksize, rworksize, vtsize, usize, ssize, asize;
		CLapackSVDStatus(char JOBU=0,char JOBVT=0, long M=0, long N=0, 
							doublecomplex *A=NULL,  long LDA=0, int ASIZE=0, 
							doublereal *S=NULL, int SSIZE=0,
							doublecomplex *U=NULL, long LDU=0, int USIZE=0,
							doublecomplex *VT=NULL, long LDVT=0, int VTSIZE=0,
							doublecomplex *WORK=NULL, long LWORK=0, int WORKSIZE=0,
							doublereal *RWORK=NULL, int RWORKSIZE=0, long INFO=0) {

			if(rwork!=NULL) delete[] rwork;
			if(work!=NULL) delete[] work;
			if(vt!=NULL) delete[] vt;
			if(u!=NULL) delete[] u;
			if(s!=NULL) delete[] s;
			if(a!=NULL) delete[] a;

			info=INFO;
			lwork=LWORK;
			ldvt=LDVT;
			ldu=LDU;
			lda=LDA;
			m=M;n=N;
			jobu=JOBU,jobvt=JOBVT;
			worksize=WORKSIZE;
			rworksize = RWORKSIZE;
			vtsize = VTSIZE;
			usize = USIZE;
			ssize = SSIZE;
			asize=ASIZE;

			work=WORK==NULL ? NULL : (doublecomplex*)memcpy(work, WORK, sizeof(doublecomplex)*worksize);
			rwork=RWORK ==NULL ? NULL : (doublereal*)memcpy(rwork, RWORK, sizeof(doublereal)*rworksize);
			vt=VT==NULL ? NULL : (doublecomplex*)memcpy(vt,VT,sizeof(doublecomplex)*vtsize); 
			u=U==NULL ? NULL : (doublecomplex*)memcpy(u,U, sizeof(doublecomplex)*usize);
			s=S==NULL ? NULL : (doublereal*)memcpy(s,S, sizeof(doublereal)*ssize);
			a=A==NULL ? NULL : (doublecomplex*)memcpy(a,A, sizeof(doublecomplex)*asize);
		}

		~CLapackSVDStatus() {
			if(rwork!=NULL) delete[] rwork;
			if(work!=NULL) delete[] work;
			if(vt!=NULL) delete[] vt;
			if(u!=NULL) delete[] u;
			if(s!=NULL) delete[] s;
			if(a!=NULL) delete[] a;
		}
};


extern CLapackSVDStatus lastzgesvdstate; 

class CExtMatrixException {
public : 
	int err_mesg;
	CExtMatrixException (int mesg) : err_mesg(mesg) {}
};

int svd(const DMatrix &A, DMatrix & U, DMatrix &W, DMatrix &V) ;
int zsvd(const ZMatrix &A, ZMatrix & U, DMatrix &S, ZMatrix &V);

ZMatrix gen_toeplitz(const ZVector &col, const ZVector &row);
DMatrix gen_toeplitz(const DVector &col, const DVector &row);

int rank(const ZMatrix &A);
int rank(const DMatrix &A);

DMatrix pinv(const DMatrix &A);
ZMatrix pinv(const ZMatrix &A);


enum EMatrixNormType { MAX_COL_SUM = 1, MAX_ROW_SUM =2, FROBENIUS=3, MAX_SVD=4 };
#endif
