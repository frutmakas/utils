#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#define NR_END 1
#define FREE_ARG char*

/* Numerical Recipes standard error handler */
void nrerror(char *error_text) {
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

/* allocate a double vector with subscript range v[nl..nh] */
double *nr_vector(long nl, long nh) {
	double *v;
	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}


/* allocate an int vector with subscript range v[nl..nh] */
int *nr_ivector(long nl, long nh) {
	int *v;
	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}


/* allocate an unsigned char vector with subscript range v[nl..nh] */
unsigned char *nr_cvector(long nl, long nh) {
	unsigned char *v;
	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}


/* allocate an unsigned long vector with subscript range v[nl..nh] */
unsigned long *nr_lvector(long nl, long nh) {
	unsigned long *v;
	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}


/* allocate a double vector with subscript range v[nl..nh] */
double *nr_dvector(long nl, long nh) {
	double *v;
	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}


/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **nr_matrix(long nrl, long nrh, long ncl, long nch) {
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}


/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **nr_dmatrix(long nrl, long nrh, long ncl, long nch) {
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}


/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
int **nr_imatrix(long nrl, long nrh, long ncl, long nch) {
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;
	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
double **nr_submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl) {
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	double **m;
	/* allocate array of pointers to rows */
	m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;
	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

/* allocate a double matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
double **nr_convert_matrix(double *a, long nrl, long nrh, long ncl, long nch) {
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	/* allocate pointers to rows */
	m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;
	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) 
		m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}


/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
double ***nr_f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh) {
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;
	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;
	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;
	for(j=ncl+1;j<=nch;j++) 
		t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) 
			t[i][j]=t[i][j-1]+ndep;
	}
	/* return pointer to array of pointers to rows */
	return t;
}

/* free a double vector allocated with vector() */
void nr_free_vector(double *v, long nl, long nh) {
	free((FREE_ARG) (v+nl-NR_END));
}


/* free an int vector allocated with ivector() */
void nr_free_ivector(int *v, long nl, long nh) {
	free((FREE_ARG) (v+nl-NR_END));
}


/* free an unsigned char vector allocated with cvector() */
void nr_free_cvector(unsigned char *v, long nl, long nh) {
	free((FREE_ARG) (v+nl-NR_END));
}


/* free an unsigned long vector allocated with lvector() */
void nr_free_lvector(unsigned long *v, long nl, long nh){
	free((FREE_ARG) (v+nl-NR_END));
}


/* free a double vector allocated with dvector() */
void nr_free_dvector(double *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}

/* free a double matrix allocated by matrix() */
void nr_free_matrix(double **m, long nrl, long nrh, long ncl, long nch) {
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

/* free a double matrix allocated by dmatrix() */
void nr_free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch){
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

/* free an int matrix allocated by imatrix() */
void nr_free_imatrix(int **m, long nrl, long nrh, long ncl, long nch) {
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

/* free a submatrix allocated by submatrix() */
void nr_free_submatrix(double **b, long nrl, long nrh, long ncl, long nch) {
	free((FREE_ARG) (b+nrl-NR_END));
}


/* free a matrix allocated by convert_matrix() */
void nr_free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch) {
	free((FREE_ARG) (b+nrl-NR_END));
}

/* free a double f3tensor allocated by f3tensor() */
void nr_free_f3tensor(double ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh) {
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

