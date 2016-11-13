/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2003/06/11 13:10:54 $
 * $Revision: 1.1.2.9 $
 * $Id: ext_matrix.cpp,v 1.1.2.9 2003/06/11 13:10:54 syed Exp $
 ********************************************************************/

#include "tools/all.h"
#include "tools/ext_matrix.h"
#include "tools/utils2lapack.h"
#include "blaswrap.h"
#include "f2c.h"
#include "clapack.h"

#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 

#pragma intrinsic(atan, exp, log10, sqrt, atan2, log, sin, tan, cos, fabs, abs ,pow)

CLapackSVDStatus lastzgesvdstate; 
 
 /*Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A =
U.W.V'. The matrixU replaces a on output. The diagonal matrix of singular values W is output
as a vector w[1..n]. The matrix V (not the transpose V' ) is output as v[1..n][1..n].*/

void svdcmp(double **a, int m, int n, double w[], double  **v) {
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z; 
	DVector rv1(n+1);
	g=scale=anorm=0.0; // Householder reduction to bidiagonal form.
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1.vect[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -nr_sign(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -nr_sign(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1.vect[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1.vect[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=nr_max(anorm,(fabs(w[i])+fabs(rv1.vect[i])));
	}
	for (i=n;i>=1;i--) { // Accumulation of right-hand transformations.
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++) // Double division to avoid possible underflow.
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1.vect[i];
		l=i;
	}
	for (i=nr_min(m,n);i>=1;i--) { //Accumulation of left-hand transformations.
		l=i+1;
	g=w[i];
	for (j=l;j<=n;j++) a[i][j]=0.0;
	if (g) {
		g=1.0/g;
		for (j=l;j<=n;j++) {
			for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
			f=(s/a[i][i])*g;
			for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
		}
		for (j=i;j<=m;j++) a[j][i] *= g;
	} else 
		for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) { // Diagonalization of the bidiagonal form: Loop over singular values, and over allowed iterations. 
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) { // Test for splitting.
				nm=l-1; // Note that rv1.vect[1] is always zero.
				if ((double)(fabs(rv1.vect[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;  // Cancellation of rv1.vect[l], if l > 1.
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1.vect[i];
					rv1.vect[i]=c*rv1.vect[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) { // Convergence.
				if (z < 0.0) { // Singular value is made nonnegative.
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) throw CExtMatrixException(NO_CONVERGENCE); //nrerror("no convergence in 30 svdcmp iterations");
			x=w[l]; // Shift from bottom 2-by-2 minor.
			nm=k-1;
			y=w[nm];
			g=rv1.vect[nm];
			h=rv1.vect[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+nr_sign(g,f)))-h))/x;
			c=s=1.0; //Next QR transformation:
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1.vect[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1.vect[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z; // Rotation can be arbitrary if z = 0.
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1.vect[l]=0.0;
			rv1.vect[k]=f;
			w[k]=x;
		}
	}
//	free_vector(rv1,1,n);
}


/*Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A =
U.W.V'. The matrix U replaces a on output. The diagonal matrix of singular values W is output
as a vector w[1..n]. The matrix V (not the transpose V' ) is output as v[1..n][1..n].

  void svdcmp(double **a, int m, int n, double w[], double  **v) {
*/

int nr_dsvd(const DMatrix &A, DMatrix & U, DMatrix &W, DMatrix &V) {
	register int i, j;
	int m = A.line+1, n=A.col+1;

	double **a = new double*[m];

	for(i=0;i<m;i++) a[i]=new double[n];
	for(i=0;i<A.line;i++) {
		for(j=0;j<A.col;j++) {
			a[i+1][j+1]=A.mat[i][j];
		}
	}

	double *w = new double[n];

	double **v = new double*[n];
	for(i=0;i<n;i++) v[i]=new double[n];

	// preparation for use with nr done ..

	try {
		svdcmp(a, m-1, n-1, w, v);
	} catch (CExtMatrixException e)  {
		delete[] w;
		for(i=0;i<m;i++) delete[] a[i];
		delete[] a;

		for(i=0;i<n;i++) delete[] v[i];
		delete[] v;
		return -1;
	} catch (...) {
		throw CExtMatrixException(UNHANDLED_EXCEPTION);
	}

	// done .. convert svdcmp result to utilitis format ..

	// a ==> U
	U = DMatrix(A.line, A.col);
	for(i=0;i<A.line; i++) {
		for(j=0;j<A.col;j++) {
			U.mat[i][j]=a[i+1][j+1];
		}
	}

	// w ==> W

	W = DMatrix(A.col, A.col);
	for(i=0;i<A.col;i++) {
		W.mat[i][i]=w[i+1];
	}


	// v ==> v
	V = DMatrix(A.col, A.col);
	for(i=0;i<A.col; i++) {
		for(j=0;j<A.col;j++) {
			V.mat[i][j]=v[i+1][j+1];
		}
	}
	return 0;
	// done ...
}

// output A = U*S*V
int dsvd(const DMatrix &A, DMatrix & U, DMatrix &S, DMatrix &V) {
	// using lapack dgesvd

	if(A.line <=0 || A.col<=0 || A.mat == NULL) throw CUtilitisException(INVALID_MATRIX);
	char JOBU = 'A', JOBVT='A';
	long M = A.line, N=A.col;
	doublereal *a = DMatrixToDoubleRealVector(A);
	long lda = M;
	int ssize = nr_min<long>(M,N); //,nr_min(M,N)

	doublereal *s = new doublereal[ssize];
	doublereal *u = new doublereal[M*M];
        long ldu = M;
	doublereal *vt = new doublereal[N*N]; long ldvt = N;
	long lwork = 5*(2*nr_min(M,N)+nr_max<long>(M,N));
	doublereal *work = new doublereal[lwork];
//	doublereal *rwork = new doublereal[(nr_max<long>(3*nr_min<long>(M,N),5*nr_min<long>(M,N)-4))];
	integer info;

/* Subroutine int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n,
            doublereal *a, integer *lda, doublereal *s, doublereal *u,
            integer * ldu, doublereal *vt, integer *ldvt, doublereal *work,
            integer *lwork, integer *info)
*/

	dgesvd_(&JOBU, &JOBVT, &M, &N,
                a, &lda, s, u,
                &ldu, vt, &ldvt, work,
                &lwork, &info);

	/*lastzgesvdstate = CLapackSVDStatus(JOBU, JOBVT, M,N, a, lda, M*N,
		                               s, ssize,
									   u, ldu, M*M, vt,ldvt, N*N,
									   work, lwork, lwork,
									   rwork, (nr_max<long>(3*nr_min<long>(M,N),5*nr_min<long>(M,N)-4)),
									   info);
*/
	//U = ZMatrix(M,M);
	DoubleRealVectorToDMatrix(u, M, M, U);

	//V = ZMatrix(N,N);
	DoubleRealVectorToDMatrix(vt, N,N, V);

	V = V.transpose_nip();
	S = DMatrix(ssize, ssize);
	for(register int i=0; i < ssize;i++) {
		S.mat[i][i]=s[i];
	}
	return info;
}

// output A = U*S*V
int zsvd(const ZMatrix &A, ZMatrix & U, DMatrix &S, ZMatrix &V) {
	// using lapack zgesvd

	if(A.line <=0 || A.col<=0 || A.mat == NULL) throw CUtilitisException(INVALID_MATRIX);
	char JOBU = 'A', JOBVT='A';
	long M = A.line, N=A.col;
	doublecomplex *a = ZMatrixToDoubleComplexVector(A);
	long lda = M;
	int ssize = nr_min<long>(M,N); //,nr_min(M,N)

	doublereal *s = new doublereal[ssize];
	doublecomplex *u = new doublecomplex[M*M]; 	long ldu = M;
	doublecomplex *vt = new doublecomplex[N*N]; long ldvt = N;
	long lwork = 5*(2*nr_min(M,N)+nr_max<long>(M,N));
	doublecomplex *work = new doublecomplex[lwork];
	doublereal *rwork = new doublereal[(nr_max<long>(3*nr_min<long>(M,N),5*nr_min<long>(M,N)-4))];
	long info;

/* Subroutine int zgesvd_(char *jobu, char *jobvt, integer *m, integer *n,
	doublecomplex *a, integer *lda, doublereal *s, doublecomplex *u, integer *
	ldu, doublecomplex *vt, integer *ldvt, doublecomplex *work, integer *lwork,
	doublereal *rwork, integer *info)*/

	zgesvd_(&JOBU, &JOBVT, &M, &N, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, &info);

	/*lastzgesvdstate = CLapackSVDStatus(JOBU, JOBVT, M,N, a, lda, M*N,
		                               s, ssize,
									   u, ldu, M*M, vt,ldvt, N*N,
									   work, lwork, lwork,
									   rwork, (nr_max<long>(3*nr_min<long>(M,N),5*nr_min<long>(M,N)-4)),
									   info);
*/
	//U = ZMatrix(M,M);
	DoubleComplexVectorToZMatrix(u, M, M, U);

	//V = ZMatrix(N,N);
	DoubleComplexVectorToZMatrix(vt, N,N, V);

	V = V.hermit();
	S = DMatrix(ssize, ssize);
	for(register int i=0; i < ssize;i++) {
		S.mat[i][i]=s[i];
	}
	return info;
}

ZMatrix gen_toeplitz(const ZVector &col, const ZVector &row){
	if (col.vect==NULL ||row.vect==NULL || row.taille==0 || col.taille ==0) throw CUtilitisException(EMPTY_VECTOR);
	ZMatrix toep(col.taille, row.taille);
	if(row.vect[0]!=col.vect[0]) cerr << "gen_toeplitz : Warning ! Column wins diagonal conflict. " << endl;
	register int i,j;
	// upper diagonal construction
	for(i=0;i<min(row.taille, col.taille);i++) {
		for(j=0;j<row.taille-i;j++) {
			toep.mat[i][j+i]=row.vect[j];
		}
	}

	// lower diagonal construction
	for(i=0;i<min(row.taille, col.taille);i++)  { // column sweep
		for(j=0;j<col.taille-i;j++) {
			toep.mat[j+i][i]=col.vect[j];
		}
	}

	return toep;
}

DMatrix gen_toeplitz(const DVector &col, const DVector &row){
	if (col.vect==NULL ||row.vect==NULL || row.taille==0 || col.taille ==0) throw CUtilitisException(EMPTY_VECTOR);
	DMatrix toep(col.taille, row.taille);
	if(row.vect[0]!=col.vect[0]) cerr << "gen_toeplitz : Warning ! Column wins diagonal conflict. " << endl;
	register int i,j;
	// upper diagonal construction
	for(i=0;i<min(row.taille, col.taille);i++) {
		for(j=0;j<row.taille-i;j++) {
			toep.mat[i][j+i]=row.vect[j];
		}
	}

	// lower diagonal construction
	for(i=0;i<min(row.taille, col.taille);i++)  { // column sweep
		for(j=0;j<col.taille-i;j++) {
			toep.mat[j+i][i]=col.vect[j];
		}
	}

	return toep;
}


int rank(const ZMatrix &A) {
/*
	s = svd(A);
	tol = max(size(A))*s(1)*eps;
	r = sum(s > tol);
*/

	ZMatrix U,V;
	DMatrix S;
	zsvd(A, U, S, V);
	int maxsize=nr_max(A.line, A.col);
	double tol = maxsize*S.mat[0][0]*DOUBLE_EPS;
	int cnt=0;
	for(int i=0;i<S.line;i++) {
		cnt += S.mat[i][i] > tol ? 1 :0;
	}
	return cnt;
}

int rank(const DMatrix &A) {
/*
	s = svd(A);
	tol = max(size(A))*s(1)*eps;
	r = sum(s > tol);
*/

	DMatrix U,S,V;
	dsvd(A, U, S, V);
	int maxsize=nr_max(A.line, A.col);
	double tol = maxsize*S.mat[0][0]*DOUBLE_EPS;
	int cnt=0;
	for(int i=0;i<S.line;i++) {
		cnt += S.mat[i][i] > tol ? 1 :0;
	}
	return cnt;
}

/*doublereal zlange_(char *norm, integer *m, integer *n, doublecomplex *a,
                   integer *lda, doublereal *work)*/

double norm(const ZMatrix &A, EMatrixNormType normtype = MAX_SVD) {
	cerr << "Warning : Matrix representationin norm is unclear !" << endl;
	if(A.line<=0 || A.col<=0 || A.mat==NULL) throw CUtilitisException(EMPTY_MATRIX);
	doublecomplex *lapack_a = ZMatrixToDoubleComplexVector(A);
	doublereal *work = new doublereal[5*A.line];
	if(work==NULL) throw CExtMatrixException(OUT_OF_MEMORY);
	char OP='I';
	switch (normtype) {
		case MAX_COL_SUM : OP='1'; break;
		case MAX_ROW_SUM : OP='I'; break;
		case MAX_SVD : OP='M'; break;
		case  FROBENIUS : OP='F'; break;
		default : OP='M'; break;
	}

	integer m=A.line, n=A.col;

	doublereal result = zlange_(&OP, &m, &n, lapack_a, &m, work);
	delete[] lapack_a;
	delete[] work;
	return result;
}

/*doublereal dlange_(char *norm, integer *m, integer *n, doublereal *a,
                   integer *lda, doublereal *work)*/

double norm(const DMatrix &A, EMatrixNormType normtype = MAX_SVD) {
	cerr << "Warning : Matrix representationin norm is unclear !" << endl;
	if(A.line<=0 || A.col<=0 || A.mat==NULL) throw CUtilitisException(EMPTY_MATRIX);
	doublereal *lapack_a = DMatrixToDoubleRealVector(A);
	doublereal *work = new doublereal[5*A.line];
	if(work==NULL) throw CExtMatrixException(OUT_OF_MEMORY);
	char OP='M';
	switch (normtype) {
		case MAX_COL_SUM : OP='1'; break;
		case MAX_ROW_SUM : OP='I'; break;
		case MAX_SVD : OP='M'; break;
		case  FROBENIUS : OP='F'; break;
		default : OP='M'; break;
	}

	integer m=A.line, n=A.col;

	doublereal result = dlange_(&OP, &m, &n, lapack_a, &m, work);
	delete[] lapack_a;
	delete[] work;
	return result;
}


ZMatrix pinv(const ZMatrix &A) {
	ZMatrix U, V;
	DMatrix S;
	// A.dim=MxN ...
	zsvd(A,U,S,V);
//	double tol = nr_max(A.line, A.col) * norm(A) * DOUBLE_EPS;
	double tol = nr_max<int>(A.line, A.col) * S.mat[0][0] * DOUBLE_EPS;
	// U.dim=M*M, S.dim=min(M,N)*min(M,N), V.dim=N*N
	// pinv(A) = VD"U.hermit
	//http://www.ee.ic.ac.uk/hp/staff/dmb/matrix/property.html#pseudoinverse
	DMatrix pS(S.line,S.col);
	DMatrix pSf(V.col, U.col);
	for(int i=0;i<S.line;i++) {
		if(S.mat[i][i]!=0) {
			double tmp= 1.0/S.mat[i][i];
			pS.mat[i][i] = tmp < tol ? 0 : tmp;
		}
	}

	InsertSubMatrixIntoMatrix(pSf, pS,0,0);

	return V*pSf*U.hermit();
}

DMatrix pinv(const DMatrix &A) {
	DMatrix U, V;
	DMatrix S;
	// A.dim=MxN ...
	dsvd(A,U,S,V);
        ofstream matlab("e:\\sandbox\\utils\\dpinv.m");
        matlaboutput(matlab, "A", A, 15);
        matlaboutput(matlab, "U", U, 15);
        matlaboutput(matlab, "S", S, 15);
        matlaboutput(matlab, "V", V, 15);
        matlab.close();

	//double tol = nr_max(A.line, A.col) * norm(A) * DOUBLE_EPS;
	double tol = nr_max(A.line, A.col) * S.mat[0][0] * DOUBLE_EPS;
	// U.dim=M*M, S.dim=min(M,N)*min(M,N), V.dim=N*N
	// pinv(A) = VD"U.hermit
	//http://www.ee.ic.ac.uk/hp/staff/dmb/matrix/property.html#pseudoinverse
	DMatrix pS(S.line,S.col);
	DMatrix pSf(V.col, U.col);
	for(int i=0;i<S.line;i++) {
		if(S.mat[i][i]!=0) {
			double tmp= 1.0/S.mat[i][i];
			pS.mat[i][i] = tmp < tol ? 0 : tmp;
		}
	}
	InsertSubMatrixIntoMatrix(pSf, pS,0,0);

	return V*pSf*U.transpose_nip();
}