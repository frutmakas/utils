/* Start of Listing */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <clapack.h>
#include "tools/utilitis.h"
#include "tools/ext_matrix.h"
#include "canal/canal_ofdm.h"

#define SIZE 4

void MAIN_(){}
void MAIN__(){}
void _MAIN_(){}


void pinvtestmain() {
	ifstream pinvdeb("samdol1.dat");
	ZMatrix gn;
	pinvdeb >> gn;
	ZMatrix pGn=pinv(gn);
	pinvdeb.close();
}

void main() {
	DMatrix TEB, SER = MIMO1OFDM(TEB);
}

void pinvmain() {
	DMatrix d(7,9);
	dRandGausStatePtr dptr;
	dRandGausInit(&dptr,7, 50);
	dmRandGaus(&dptr, d);
	ZMatrix z(7,9);
	zRandGausStatePtr zptr;
	zRandGausInit(&zptr, 7, 50.1);
	zmRandGaus(&zptr, z);

	ofstream matlab("e:\\sandbox\\utils\\pinvtest.m");
	matlaboutput(matlab, "d", d, 15);
	matlaboutput(matlab, "z", z, 15);

	DMatrix pd=pinv(d);
	ZMatrix pz=pinv(z);
	matlaboutput(matlab, "pd", pd, 15);
	matlaboutput(matlab, "pz", pz, 15);
	cout << d << pd << endl << endl;
	cout << z << pz << endl << endl;
}

void test_main() {
/*	DVector a(5), b(12);
	int i=0;
	for(i=0;i<a.taille;i++) a.vect[i]=2*i+1;
	for(i=0;i<b.taille;i++) b.vect[i]=i+3;
	DMatrix t1 = gen_toeplitz(a, b);
	DMatrix t2 = gen_toeplitz(b,a);
	cout << a << b << t1 << t2 ;
*/
	ZVector a(5), b(7);
	zRandGausStatePtr ptr;
	zRandGausInit(&ptr, 15, 100);
	zbRandGaus(&ptr, a);
	zbRandGaus(&ptr, b);
	ZVector c = a.conv_nip(b);
	ofstream mmo("j:\\work\\cances.old\\conv_test.m");
	matlaboutput(mmo, "va", a,15);
	matlaboutput(mmo, "vb", b,15);
	matlaboutput(mmo, "vc", c,15);

	ZVector d=a.corr_nip(b);
	matlaboutput(mmo, "vd", d,15);
	ZVector e=a.xcorr_nip();
	matlaboutput(mmo, "ve", e,15);
	a.xcorr();
	matlaboutput(mmo, "vf", a,15);

	mmo.close();


	ZMatrix C, D;
	ifstream is("xC.txt");
	is >> C >> D;
	mmo.open("j:\\work\\cances.old\\initdebug.m");
	matlaboutput(mmo, "xC", C);
	matlaboutput(mmo, "xD", D);
	C=C.hermit();
	D=D.hermit(); 
	matlaboutput(mmo, "tC", C);
	matlaboutput(mmo, "tD", D);
	mmo << "C = tC; D=tD;";
	mmo.close();

	int Lf=5;
	double sigma =0.3972;
	int fd=100;
	int training_length=10;
//	ZMatrix CF = ofdm_time_frequency_channel_estimation(C,D, training_length, Lf, sigma, fd);


}

main_menu( )
{
     char JOBU;
     char JOBVT;

     int i;
    
     integer M = SIZE;
     integer N = SIZE;
    
     integer LDA = M;
     integer LDU = M;
     integer LDVT = N;
     integer LWORK;
     integer INFO;
   
     integer mn = min( M, N );
    
     integer MN = max( M, N );
     
     double a[SIZE*SIZE] = { 16.0, 5.0, 9.0 , 4.0, 2.0, 11.0, 7.0 , 14.0, 3.0, 10.0, 6.0, 15.0, 13.0, 8.0, 12.0, 1.0};
     double s[SIZE];
     double wk[201];
     double uu[SIZE*SIZE];
     double vt[SIZE*SIZE];
     
       JOBU = 'A';
     
       JOBVT = 'A';

    LWORK =  201;

    printf("\nbefore = [");
	 for(i=0;i<SIZE;i++) {
		 printf("\n[");

		 for(int j=0;j<SIZE;j++) {
			 printf("%lf\t", a[i*SIZE+j]);
		 }
		 printf("]");
	 }
	 printf("\n];");

/* Subroutine int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
        doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
        ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
        integer *info)
*/
   dgesvd_( &JOBU, &JOBVT, &M, &N, a, &LDA, s, uu, 
          &LDU, vt, &LDVT, wk, &LWORK, &INFO);

    printf("\n INFO=%d", INFO );          

    for ( i= 0; i< SIZE; i++ ) {
        printf("\n s[ %d ] = %f", i, s[ i ] );
    }

	 printf("\nAfter = [");
	 for(i=0;i<SIZE;i++) {
		 printf("\n[");

		 for(int j=0;j<SIZE;j++) {
			 printf("%lf\t", a[i*SIZE+j]);
		 }
		 printf("]");
	 }
	 printf("\n];");

    for ( i= 0; i< SIZE; i++ ) {
        printf("\n s[ %d ] = %f", i, s[ i ] );
    }
	printf("\n");

	 printf("\ncuu = [");
	 for(i=0;i<SIZE;i++) {
		 printf("\n[");

		 for(int j=0;j<SIZE;j++) {
			 printf("%lf\t", uu[i*SIZE+j]);
		 }
		 printf("]");
	 }
	 printf("\n];");

 	 printf("\ncvt = [");
	 for(i=0;i<SIZE;i++) {
		 printf("\n[");

		 for(int j=0;j<SIZE;j++) {
			 printf("%lf\t", vt[i*SIZE+j]);
		 }
		 printf("]");
	 }
	printf("\n];\n");


   doublecomplex ca[SIZE*SIZE];

	ca[0].r = 16.0;
	ca[1].r = 5.0; 
	ca[2].r = 9.0 ;
	ca[3].r = 4.0;
	ca[4].r = 2.0;
	ca[5].r = 11.0;
	ca[6].r = 7.0 ;
	ca[7].r = 14.0;
	ca[8].r = 3.0;
	ca[9].r = 10.0;
	ca[10].r = 6.0;
	ca[11].r = 15.0;
	ca[12].r = 13.0;
	ca[13].r = 8.0;
	ca[14].r = 12.0;
	ca[15].r = -81.0;

	ca[0].i = 6.0;
	ca[1].i = 3.0; 
	ca[2].i = 7.0 ;
	ca[3].i = 2.0;
	ca[4].i = -44.0;
	ca[5].i = 2.0;
	ca[6].i = 17.0 ;
	ca[7].i = 1.0;
	ca[8].i = 10.0;
	ca[9].i = 6.0;
	ca[10].i = 16.0;
	ca[11].i = 1.0;
	ca[12].i = -3.0;
	ca[13].i = 4.7;
	ca[14].i = 5.0;
	ca[15].i = 9.0;

    doublecomplex cuu[SIZE*SIZE];
    doublecomplex cvt[SIZE*SIZE];
	doublereal rwk[5*SIZE];
	doublecomplex cwk[201];

/* Subroutine int zgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
	doublecomplex *a, integer *lda, doublereal *s, doublecomplex *u, integer *
	ldu, doublecomplex *vt, integer *ldvt, doublecomplex *work, integer *lwork, 
	doublereal *rwork, integer *info)*/
	 printf("\nbefore = [");
	 for(i=0;i<SIZE;i++) {
		 printf("\n[");

		 for(int j=0;j<SIZE;j++) {
			 printf("%f+%fi\t", ca[i*SIZE+j].r, ca[i*SIZE+j].i);
		 }
		 printf("]");
	 }
	 printf("\n];");

    zgesvd_( &JOBU, &JOBVT, &M, &N, ca, &LDA, s, cuu, 
          &LDU, cvt, &LDVT, cwk, &LWORK, rwk, &INFO);

	 printf("\n INFO=%d ", INFO );

	 printf("\nAfter = [");
	 for(i=0;i<SIZE;i++) {
		 printf("\n[");

		 for(int j=0;j<SIZE;j++) {
			 printf("%lf+%lfi\t", ca[i*SIZE+j].r, ca[i*SIZE+j].i);
		 }
		 printf("]");
	 }
	 printf("\n];");

    for ( i= 0; i< SIZE; i++ ) {
        printf("\n s[ %d ] = %f", i, s[ i ] );
    }
	printf("\n");

	 printf("\ncuu = [");
	 for(i=0;i<SIZE;i++) {
		 printf("\n[");

		 for(int j=0;j<SIZE;j++) {
			 printf("%lf+%lfi\t", cuu[i*SIZE+j].r, cuu[i*SIZE+j].i);
		 }
		 printf("]");
	 }
	 printf("\n];");

 	 printf("\ncvt = [");
	 for(i=0;i<SIZE;i++) {
		 printf("\n[");

		 for(int j=0;j<SIZE;j++) {
			 printf("%lf+%lfi\t", cvt[i*SIZE+j].r, cvt[i*SIZE+j].i);
		 }
		 printf("]");
	 }
	printf("\n];\n");


    return 0;
}

/* End of Listing */

