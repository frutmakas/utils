#ifndef UTILITIS_H
#define UTILITIS_H

#include <stdio.h>
#include <iostream.h>

#ifndef M_PI
#define M_PI 3.14159265359
#endif

#define min(x,y)  (((x)<(y)) ? (x):(y))
#define max(x,y)  (((x)>(y)) ? (x):(y))
#define puiss2(a) ((a)*(a))

typedef struct vector{ 
        double * vect;
		long taille;
} vector;

typedef struct long_vector{ 
        long * vect;
		long taille;
} long_vector;

typedef struct cmplx {
			double Re, Im;
} cmplx;

typedef struct cmplx_vect{
		cmplx * vect;
		long taille;
} cmplx_vect;

typedef struct long_matrice{
		long ** mat;
		long lin, col;
} long_matrice;

typedef struct matrice{
		double ** mat;
		long lin, col;
} matrice;

typedef struct cmplx_mat {
			cmplx **mat;
			long lin, col;
} cmplx_mat;

cmplx        init_cmplx(double a, double b);
cmplx        conjug(cmplx a);
double       modul2(cmplx a);
cmplx        cmplx_add(cmplx a,cmplx b);
cmplx        cmplx_sub(cmplx a,cmplx b);
cmplx        cmplx_mul(cmplx a,cmplx b);
double       argum(cmplx x);
vector       init_vect(long ,double);
vector       init_vect(long);
long_vector  init_long_vect(long,long);
void         liberer(vector a);
void         liberer(long_vector);
cmplx_vect   init_cmplx_vect(long , cmplx);
cmplx_vect   init_cmplx_vect(long);
void         liberer(cmplx_vect a);
long_matrice init_long_mat(long LIN, long COL, long valeur);
matrice      init_matrice(long LIN, long COL, double valeur);
void         liberer(matrice a);
void         liberer(long_matrice a);
vector       col_mat2vect(matrice a, long col);
cmplx_vect   conjug(cmplx_vect a);
void         Err_Message(char Msg[]);
void	     cmplx_s_mul_vect (cmplx_vect a, cmplx_vect b, cmplx c);
void	     cmplx_mul_vect(cmplx_vect a, cmplx_vect b, cmplx_vect c,double );
void	     cmplx_mul_vect(cmplx_vect a, cmplx_vect b, cmplx_vect c);
void         sub_matrice(matrice a, matrice b, matrice c);
void	     print(cmplx a);
void         print(cmplx_vect a);
void         print(vector a);
void         print(long_vector a);
void         print(matrice a);
void         print(long_matrice a);
void		 print(cmplx_mat a);
double       n_gaussien( double moy , double var );
void         random_bit(long_vector a);
void         rand_M(long_vector a, long M);
void         QPSK_mod(long_vector data,cmplx_vect x);
void         QPSK_demod(cmplx_vect x, long_vector data_hat);
void         convol (vector x, vector h, vector y);
void         convol (cmplx_vect x, vector h, cmplx_vect y);
void         convol (cmplx_vect x, cmplx_vect h, cmplx_vect y);
cmplx_vect   correlation(cmplx_vect a, cmplx_vect b);
vector       File2Vect (char Fname[]);
long_vector  File2long_vect (char Fname[]);
void         File2long_matrice (char Fname[],long_matrice a);
cmplx_vect   File2Vect_cmplx (char Fname[]);
void         File_Print (char Fname[],long_vector x);
void         File_Print (char Fname[],vector x);
long          File_Print(char * nom,cmplx_vect a);
void         Fill_vect(vector a, double val);
void         zero_pad(vector x, vector y, long rate);
void         zero_pad(cmplx_vect x, cmplx_vect y, long rate);
void         decimate (vector x, vector y, long delai,long rate);
void         decimate (cmplx_vect x, cmplx_vect y, long delai,long rate);
long         indice_max(vector a);
long         indice_min(vector a);
void         shift_right(long_vector a, long decal);
void         shift_right(vector a, long decal);
void         shift_right(cmplx_vect a, long decal);
void         shift_left(long_vector a,long decal);
void         shift_left(vector a,long decal);
void         shift_left(cmplx_vect a,long decal);
void         lin_vect(vector a,double initial, double step);
cmplx_vect   cmplx_exp(vector a);
vector       copy(vector a);
cmplx_vect   copy(cmplx_vect a);
long_vector  copy(long_vector a);
matrice      copy(matrice a);
long_matrice copy(long_matrice a);
void         awgn(vector a, double SNR, double factor);
void         awgn(cmplx_vect a, double SNR, double factor);
vector       BPSK_MOD(long_vector a);
long_vector	 BPSK_DEM(vector a);
long         comparer(long_vector a, long_vector b);

#endif
