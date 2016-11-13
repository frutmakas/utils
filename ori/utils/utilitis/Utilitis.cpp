#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utilitis.h"


cmplx init_cmplx(double a, double b)
{
	cmplx c;
	c.Re=a;
	c.Im=b;
	return c;
}
/*-----------------------------------------------*/
cmplx conjug(cmplx a)
{
	a.Im=-a.Im;
	return a;
}
/*-----------------------------------------------*/
double modul2(cmplx a)
{
	return a.Re*a.Re+a.Im*a.Im;
}
/*-----------------------------------------------*/
cmplx cmplx_add(cmplx a,cmplx b)
{
	cmplx c;
	c.Re=a.Re+b.Re;
	c.Im=a.Im+b.Im;
	return c;
}
/*-----------------------------------------------*/
cmplx cmplx_sub(cmplx a,cmplx b)
{
	cmplx c;
	c.Re=a.Re-b.Re;
	c.Im=a.Im-b.Im;
	return c;
}
/*-----------------------------------------------*/
cmplx cmplx_mul(cmplx a,cmplx b)
{
	cmplx c;
	double temp;
	temp=a.Re*b.Re-a.Im*b.Im;
	c.Im=a.Re*b.Im+a.Im*b.Re;
	c.Re=temp;
	return c;
}
/*-----------------------------------------------*/
double argum(cmplx x)
{
	double phase;
	phase =((x.Re==0.0) && (x.Im ==0.0)) ?
		0.0 : atan2(x.Im,x.Re);
	return phase;
}
/*----------------------------------------------*/
vector init_vect(long taille ,double valeur)
{
	vector a;
	long i;
	a.vect=(double *)calloc(taille , sizeof(double));
	if (a.vect==NULL) Err_Message("Memory allocation Error in init_vect");
	a.taille = taille;
	for (i=0; i<taille ; i++) a.vect[i]=valeur;
	return a;
}
/*----------------------------------------------*/
vector init_vect(long taille)
{
	vector a;
	a.vect=(double *)calloc(taille , sizeof(double));
	if (a.vect==NULL) Err_Message("Memory allocation Error in init_vect");
	a.taille = taille;
	return a;
}


/*----------------------------------------------*/
long_vector init_long_vect(long taille ,long valeur=0)
{
	long_vector a;
	long i;
	a.vect=(long *)calloc(taille , sizeof(long));
	if (a.vect==NULL) Err_Message("Memory allocation Error in init_long_vect");
	a.taille = taille;
	for (i=0; i<taille ; i++) a.vect[i]=valeur;
	return a;
}

/*----------------------------------------------*/
void liberer(vector a)
{
	free(a.vect);
}

/*----------------------------------------------*/
void liberer(long_vector a)
{
	free(a.vect);
}
/*-------------------------------------------------*/
cmplx_vect init_cmplx_vect(long taille, cmplx val)
{
	cmplx_vect a;
	long i;
	a.vect=(cmplx *)calloc(taille, sizeof(cmplx));
	if (a.vect==NULL) Err_Message("Memory allocation Error in init_cmplx_vect");
	for (i=0; i<taille ; i++) 
		a.vect[i]=val;
	a.taille=taille;
	return a;
}
/*-------------------------------------------------*/
cmplx_vect init_cmplx_vect(long taille)
{
	cmplx_vect a;
	a.vect=(cmplx *)calloc(taille, sizeof(cmplx));
	if (a.vect==NULL) Err_Message("Memory allocation Error in init_cmplx_vect");
	a.taille=taille;
	return a;
}

/*------------------------------------------------*/
void liberer(cmplx_vect a)
{
	free(a.vect);
}

/*-----------------------------------------------*/
long_matrice init_long_mat(long LIN, long COL, long valeur)
{
	long_matrice a;
	long i,j;
	a.mat=(long * *)calloc(LIN,sizeof(long *));
	if (a.mat==NULL) Err_Message("Memory allocation Error in init_matrice");
	a.mat[0]=(long *)calloc(LIN*COL, sizeof(long));
	if (a.mat[0]==NULL) Err_Message("Memory allocation Error in init_matrice");
	for (i=1; i<LIN; i++)
		a.mat[i]=COL+a.mat[i-1];
	for (i=0; i<LIN; i++)
		for (j=0 ; j<COL ; j++)
			a.mat[i][j]=valeur;
	a.lin=LIN;
	a.col=COL;
	return a;
}

/*---------------------------------------------*/
void liberer(long_matrice a)
{
	free (a.mat[0]);
	free(a.mat);
}

/*-----------------------------------------------*/
matrice init_matrice(long LIN, long COL, double valeur)
{
	matrice a;
	long i,j;
	a.mat=(double * *)calloc(LIN,sizeof(double *));
	if (a.mat==NULL) Err_Message("Memory allocation Error in init_matrice");
	a.mat[0]=(double *)calloc(LIN*COL, sizeof(double));
	if (a.mat[0]==NULL) Err_Message("Memory allocation Error in init_matrice");
	for (i=1; i<LIN; i++)
		a.mat[i]=COL+a.mat[i-1];
	for (i=0; i<LIN; i++)
		for (j=0 ; j<COL ; j++)
			a.mat[i][j]=valeur;
	a.lin=LIN;
	a.col=COL;
	return a;
}
/*---------------------------------------------*/
void liberer(matrice a)
{
	free (a.mat[0]);
	free(a.mat);
}

/*-------------------------------------------------------------*/
vector col_mat2vect(matrice a, long col)
{
	vector b;
	long i;
	b=init_vect(a.lin,0);
	for (i=0; i<a.lin; i++)
		b.vect[i]=a.mat[i][col];
	return b;
}

/*-------------------------------------------------------------*/
cmplx_vect conjug(cmplx_vect a)
{
	long	i;
	cmplx_vect	c;
	c=init_cmplx_vect(a.taille,init_cmplx(0,0));
	for(i=0; i<a.taille; i++)
		c.vect[i]=conjug(a.vect[i]);
	return c;
}

/*-------------------------------------------------------------*/
void Err_Message(char Msg[])
{
	cout << Msg << endl;
	exit (2);
}
/*------------------------------------------------*/
void cmplx_s_mul_vect (cmplx_vect a, cmplx_vect b, cmplx c)
{
	long i;
	for (i=0 ; i<a.taille ; i++)
		b.vect[i]=cmplx_mul(a.vect[i],c);
}
/*----------------------------------------------*/
void cmplx_mul_vect(cmplx_vect a, cmplx_vect b, cmplx_vect c,double conj=1.0)
{
	long i,size;
	size=min(min(a.taille,b.taille),c.taille);
	if (conj==1)
		for (i=0 ; i<size ; i++)
			c.vect[i]=cmplx_mul(a.vect[i],b.vect[i]);
	else
		for (i=0 ; i<size ; i++)
			c.vect[i]=cmplx_mul(a.vect[i],conjug(b.vect[i]));
}
/*----------------------------------------------*/
void sub_matrice(matrice a, matrice b, matrice c)
{
	long i,j;
	if(a.lin!=b.lin || a.col!=b.col || a.lin!=c.lin || a.col != c.col)
		Err_Message("taille incompatible dans sub_matrice");
	for (i=0; i<a.lin; i++)
		for (j=0; j<a.col; j++)
			c.mat[i][j]=a.mat[i][j]-c.mat[i][j];
}
/*----------------------------------------------*/

void print (cmplx a)
{
	if (a.Im < 0)
		cout <<a.Re<<"-j"<<-a.Im<<endl;
	else
		cout <<a.Re<<"+j"<<a.Im<<endl;
}
/*----------------------------------------------*/
void print(cmplx_vect a)
{
	long i;
	cout << endl;
	for (i=0 ; i<a.taille ; i++)
		print(a.vect[i]);
	cout << endl;
}

void print(cmplx_mat a){
	long i,k;
	cout << endl;
	for(i=0;i<a.lin;i++) {
		for(k=0;k<a.col;k++)
			print(a.mat[k][i]);
		cout << endl;
	}
	cout << endl;

}
/*---------------------------------------------*/
void print(vector a)
{
	long i;
	cout << endl;
	for (i=0; i<a.taille ; i++) cout << "   " << a.vect[i] << endl;
	cout << endl;
}

/*---------------------------------------------*/
void print(long_vector a)
{
	long i;
	cout << endl;
	for (i=0; i<a.taille ; i++) cout << "  " << a.vect[i] << endl;
	cout << endl;
}
/*----------------------------------------------*/
void print(matrice a)
{
	long i,j;
	cout << endl;
	for (i=0; i<a.lin ; i++){
		for (j=0; j<a.col ; j++)
			printf("%.4f ",a.mat[i][j]);
		cout << endl;
	}
}
/*----------------------------------------------*/
void print(long_matrice a)
{
	long i,j;
	cout << endl;
	for (i=0; i<a.lin ; i++){
		for (j=0; j<a.col ; j++)
			printf("%5d ",a.mat[i][j]);
		cout << endl;
	}
}

/*----------------------------------------------*/
double n_gaussien( double moy , double var ) //___________LOI GAUSSIENNE
//	retourne un evenement suivant une loi gaussienne de moyenne moy
//	et de variance var
{
	static long iset;
	static double gset;
	double fac,rsq,v1,v2;
	if (iset==0){
		do{
			v1=rand()*2.0/RAND_MAX-1.0;
			v2=rand()*2.0/RAND_MAX-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq ==0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac*sqrt(var)+moy;
	} else {
		iset=0;
		return gset*sqrt(var)+moy;
	}
}
/*----------------------------------------------*/
/* Generation d'un vecteur de bits aleatoires  */
void random_bit(long_vector a)
{
	long i;
	for (i=0; i<a.taille; i++)
		a.vect[i]=rand() & 1;
}

/*-------------------------------------------
*   Génération d'un vecteur de symbole aléatoires entre 0 et M-1
*   Attention : M doit etre 2**k
*/
void rand_M(long_vector a, long M)
{
	long i;
	for (i=0; i<a.taille; i++)
		a.vect[i]=rand() & M-1;
}
/*-------------------------------------------------------*/
/*  Modulation QPSK shifted, codage Gray  */
void QPSK_mod(long_vector data,cmplx_vect x)
{
	long i;
	double amp;
	amp=sqrt(2.0)/2.0;
	if (x.taille*2<data.taille)
		Err_Message("incompatible sizes in QPSK_MOD function");
	for (i=0; i<data.taille/2; i++){
		x.vect[i].Re=(data.vect[2*i]==0) ? amp : -amp;
		x.vect[i].Im=(data.vect[2*i+1]==0) ? amp : -amp;
	}
}
/*-------------------------------------------------------*/
/*  Demodulation QPSK shifted, codage Gray  */
void QPSK_demod(cmplx_vect x, long_vector data_hat)
{
	long i;
	for (i=0; i<data_hat.taille/2; i++){
		data_hat.vect[2*i] = (x.vect[i].Re >= 0) ? 0 : 1;
		data_hat.vect[2*i+1] = (x.vect[i].Im >= 0) ? 0 : 1;
	}
}
/*-------------------------------------------------------*/
void convol (vector x, vector h, vector y)
{
  long m;
  long j;
  double sum;

  /* Loop over output points */
  for (m = y.taille-1; m >= 0 ; m--) {

	 /* Convolution */
	 sum = 0.0;
	 for (j = 0; j < h.taille; ++j) {
		if (m-j>=0 && m-j<x.taille) sum += h.vect[j] * x.vect[m-j];
	 }
	 y.vect[m] = sum;
  }
}
/*-------------------------------------------------------*/
/* convolution : entree et sortie complexes, filtre réel  */
/* La sortie peut etre le meme vector que l'entrée */
void convol (cmplx_vect x, vector h, cmplx_vect y)
{
  long m;
  long j;
  double sum_R, sum_I;

  /* Loop over output points */
  for (m = y.taille-1; m >= 0 ; m--) {

	 /* Convolution */
	 sum_R = sum_I = 0.0;
	 for (j = 0; j < h.taille; ++j) {
		if (m-j>=0 && m-j<x.taille){
			sum_R += h.vect[j] * x.vect[m-j].Re;
			sum_I += h.vect[j] * x.vect[m-j].Im;
		}
	 }
	 y.vect[m].Re = sum_R;
	 y.vect[m].Im = sum_I;
  }
}
/*-------------------------------------------------------*/
/* convolution : entree, sortie et filtre complexes  */
/* La sortie peut etre le meme vector que l'entrée */
void convol (cmplx_vect x, cmplx_vect h, cmplx_vect y)
{
  long m;
  long j;
  cmplx sum;

  /* Loop over output points */
  for (m = y.taille-1; m >= 0 ; m--) {

	 /* Convolution */
	 sum = init_cmplx(0,0);
	 for (j = 0; j < h.taille; ++j)
		 if (m-j>=0 && m-j<x.taille)
			sum = cmplx_add(sum,cmplx_mul(h.vect[j] , x.vect[m-j]) );
	 y.vect[m] = sum;
  }
}

/*------------------------------------------------------*/
cmplx_vect	correlation(cmplx_vect a, cmplx_vect b)
{
	cmplx_vect R;
	long	n,m;
	cmplx	sum;
	R=init_cmplx_vect(a.taille+b.taille-1, init_cmplx(0,0));
	for (m=b.taille-1 ; m>-a.taille ; m--){
		sum=init_cmplx(0,0);
		for (n=0;n<a.taille; n++)
			if (m+n>=0 && m+n<b.taille)
				sum=cmplx_add(sum,cmplx_mul(a.vect[n],conjug(b.vect[m+n])));
		R.vect[b.taille-1-m]=sum;
	}
	return R;
}
/*------------------------------------------------------*/
vector File2Vect (char Fname[])
{
  FILE *fp;
  long Nval=0;
  double tmp;
  vector x;
  double *xp;
/* Open the data file */
  fp = fopen (Fname, "r");
  if (fp == NULL) {
	 cout << "File ReadData: Cannot open file "<< Fname<<endl;
	 exit (1);
  }
/* Read the values */
  while (!feof(fp)){
	fscanf(fp,"%lf\n",&tmp);
	Nval++;
  }
  fseek(fp,0L,SEEK_SET);
  x=init_vect(Nval,0);
  xp=x.vect;
  while (!feof(fp))
		fscanf(fp,"%lf\n",xp++);
/* Close the file */
  fclose (fp);
  x.taille=Nval;
  return x;
}
/*------------------------------------------------------*/
long_vector File2long_vect (char Fname[])
{
  FILE *fp;
  long Nval=0;
  long tmp;
  long_vector x;
  long *xp;
/* Open the data file */
  fp = fopen (Fname, "r");
  if (fp == NULL) {
	 cout << "File ReadData: Cannot open file "<< Fname<<endl;
	 exit (1);
  }
/* Read the values */
  while (!feof(fp)){
	fscanf(fp,"%d\n",&tmp);
	Nval++;
  }
  fseek(fp,0L,SEEK_SET);
  x=init_long_vect(Nval,0);
  xp=x.vect;
  while (!feof(fp))
		fscanf(fp,"%d\n",xp++);
/* Close the file */
  fclose (fp);
  x.taille=Nval;
  return x;
}
/*------------------------------------------------------*/
void File2long_matrice (char Fname[],long_matrice a)
{
  FILE *fp;
  long Nval=0,i,j;
  long tmp;
/* Open the data file */
  fp = fopen (Fname, "r");
  if (fp == NULL) {
	 cout << "File ReadData: Cannot open file "<< Fname<<endl;
	 exit (1);
  }
/* Read the values */
  while (!feof(fp)){
	fscanf(fp,"%d\n",&tmp);
	Nval++;
  }
  if (Nval != a.col*a.lin) Err_Message("File size incompatible!");
  fseek(fp,0L,SEEK_SET);
  for (i=0; i<a.lin; i++)
	  for (j=0; j<a.col ; j++)
		  fscanf(fp,"%d\n",&a.mat[i][j]);
/* Close the file */
  fclose (fp);
}
/*------------------------------------------------------*/
cmplx_vect File2Vect_cmplx (char Fname[])
{
  FILE *fp;
  long Nval=0;
  double temp1, temp2;
  cmplx *p;
  cmplx_vect x;

/* Open the data file */
  fp = fopen (Fname, "r");
  if (fp == NULL) {
	 cout << "File ReadData: Cannot open file "<< Fname<<endl;
	 exit (1);
  }
/* Read the values */
  while (!feof(fp)){
	fscanf(fp,"%lf %lf\n",&temp1, &temp2);
	Nval++;
  }
  fseek(fp,0L,SEEK_SET);
  x=init_cmplx_vect(Nval,init_cmplx(0,0));
  p=x.vect;
  while (!feof(fp)){
	fscanf(fp,"%lf %lf\n",&p->Re, &p->Im);
	p++;
  }
/* Close the file */
  fclose (fp);
  x.taille=Nval;
  return x;
}
/*-----------------------------------------------------------*/
void File_Print (char Fname[],vector x)
{
	FILE *fp;
	long i;
	fp=fopen (Fname, "w");
	for (i=0; i<x.taille ; i++)
		fprintf(fp,"%lf\n",*x.vect++);
	fclose(fp);
}
/*-----------------------------------------------------------*/
void File_Print (char Fname[],long_vector x)
{
	FILE *fp;
	long i;
	fp=fopen (Fname, "w");
	for (i=0; i<x.taille ; i++)
		fprintf(fp,"%ld\n",*x.vect++);
	fclose(fp);
}
/* -------------------------------------------------------*/
/* ecrire des données complexes dans  un fichier  */
long File_Print(char nom[],cmplx_vect a)
{
	FILE *fp;
	long i;
	if ((fp = fopen(nom,"w")) ==NULL) return 1;
	else {
		for (i=0; i<a.taille; i++)
			fprintf(fp,"%lf  %lf\n",a.vect[i].Re, a.vect[i].Im);
		fclose(fp);
		return 0;
	}
}
/*--------------------------------------------------------*/
void Fill_vect(vector a, double val)
{
	long i;
	for (i=0; i<a.taille ; i++)
		*a.vect++ = val;
}
/*--------------------------------------------------------*/
/* Fill a comlex vector */
void Fill_cmplx_vect(cmplx_vect a, cmplx val)
{
	long i;
	for (i=0; i<a.taille ; i++)
		*a.vect++ = val;
}
/*--------------------------------------------------------*/
void zero_pad(vector x, vector y, long rate)
{
	long i;
	if (y.taille<x.taille*rate)
		Err_Message("Error in zero_pad function: incompatible output vector size");
	Fill_vect(y,0);
	for (i=0; i<y.taille ; i++)
		y.vect[i*rate]=x.vect[i];
}
/*----------------------------------------------------*/
/* zero padding for complex signals  */
void zero_pad(cmplx_vect x, cmplx_vect y, long rate)
{
	long i;
	if (y.taille<x.taille*rate)
		Err_Message("Error in zero_pad function: incompatible output vector size");
	Fill_cmplx_vect(y,init_cmplx(0,0));
	for (i=0; i<x.taille ; i++)
		y.vect[i*rate]=x.vect[i];
}
/*--------------------------------------------------------*/
void decimate (vector x, vector y, long delai,long rate)
{
	long i,ind;
	for (i=0 ; i<y.taille ; i ++){
		ind=i*rate+delai;
		if (ind<x.taille)
			*y.vect++ = x.vect[ind];
		else *y.vect++ = 0;
	}
}
/*--------------------------------------------------------*/
void decimate (cmplx_vect x, cmplx_vect y, long delai,long rate)
{
	long i,ind;
	for (i=0 ; i<y.taille ; i ++){
		ind=i*rate+delai;
		if (ind<x.taille)
			*y.vect++ = x.vect[i*rate+delai];
		else
			*y.vect++ = init_cmplx(0,0);
	}
}
/*---------------------------------------------------------*/
long indice_min(vector a)
{
	long i, indice;
	double mini;
	mini=a.vect[0];
	indice=0;
	for (i=1 ; i<a.taille ; i++)
		if (a.vect[i]<mini) {
			mini=a.vect[i];
			indice=i;
		}
	return indice;
}
/*---------------------------------------------------------*/
long indice_max(vector a)
{
	long i, indice;
	double maxi;
	maxi=a.vect[0];
	indice=0;
	for (i=1 ; i<a.taille ; i++)
		if (a.vect[i]>maxi) {
			maxi=a.vect[i];
			indice=i;
		}
	return indice;
}
/*-----------------------------------------------------------*/
void shift_right(long_vector a, long decal)
{
	long *des,*src;
	long i;
	des=&a.vect[a.taille-1];
	src=&a.vect[a.taille-decal-1];
	for (i=0;i<a.taille-decal;i++)
		*des-- = *src--;
	for (i=0; i<decal;i++)
		*des-- = 0;
}
/*-----------------------------------------------------------*/
void shift_right(vector a, long decal)
{
	double *des,*src;
	long i;
	des=&a.vect[a.taille-1];
	src=&a.vect[a.taille-decal-1];
	for (i=0;i<a.taille-decal;i++)
		*des-- = *src--;
	for (i=0; i<decal;i++)
		*des-- = 0;
}
/*-----------------------------------------------------------*/
void shift_right(cmplx_vect a, long decal)
{
	cmplx *des,*src;
	long i;
	des=&a.vect[a.taille-1];
	src=&a.vect[a.taille-decal-1];

	for (i=0;i<a.taille-decal;i++)
		*des-- = *src--;
	for (i=0; i<decal;i++)
		*des-- = init_cmplx(0,0);
}
/*-----------------------------------------------------------*/
void shift_left(long_vector a,long decal)
{
	long *des,*src,i;
	des=&a.vect[0];
	src=&a.vect[decal];
	for (i=0;i<a.taille-decal;i++)
		*des++ = *src++;
	for (i=0; i<decal;i++)
		*des++ = 0;
}
/*-----------------------------------------------------------*/
void shift_left(vector a,long decal)
{
	double *des,*src;
	long i;
	des=&a.vect[0];
	src=&a.vect[decal];
	for (i=0;i<a.taille-decal;i++)
		*des++ = *src++;
	for (i=0; i<decal;i++)
		*des++ = 0;
}
/*-----------------------------------------------------------*/
void shift_left(cmplx_vect a,long decal)
{
	cmplx *des,*src;
	long i;
	des=a.vect;
	src=&a.vect[decal];

	for (i=0;i<a.taille-decal;i++)
		*des++ = *src++;
	for (i=0; i<decal;i++)
		*des++ = init_cmplx(0,0);
}
/*------------------------------------*/
void lin_vect(vector a,double initial, double step)
{
	long i;
	a.vect[0]=initial;
	for (i=1; i<a.taille; i++)
		a.vect[i] = a.vect[i-1] + step;
}
/*-----------------------------------------*/
/* cette routine calcule exp(j*input)  */
cmplx_vect cmplx_exp(vector a)
{
	long i;
	cmplx_vect x;
	x=init_cmplx_vect(a.taille,init_cmplx(0,0));
	for (i=0 ; i<a.taille ; i++){
		x.vect[i].Re=cos(*a.vect);
		x.vect[i].Im=sin(*a.vect++);
	}
	x.taille=a.taille;
	return x;
}
/*-----------------------------------------*/
vector copy(vector a)
{
	vector b;
	long i;
	b=init_vect(a.taille,0);
	b.taille=a.taille;
	for (i=0; i<a.taille; i++)
		b.vect[i]=a.vect[i];
	return b;
}
/*-----------------------------------------*/
cmplx_vect copy(cmplx_vect a)
{
	cmplx_vect b;
	long i;
	b=init_cmplx_vect(a.taille,init_cmplx(0,0));
	b.taille=a.taille;
	for (i=0; i<a.taille; i++)
		b.vect[i]=a.vect[i];
	return b;
}
/*-----------------------------------------*/
long_vector copy(long_vector a)
{
	long_vector b;
	long i;
	b=init_long_vect(a.taille,0);
	b.taille=a.taille;
	for (i=0; i<a.taille; i++)
		b.vect[i]=a.vect[i];
	return b;
}
/*-----------------------------------------*/
matrice copy(matrice a)
{
	matrice b;
	double *p_a,*p_b;
	long i;
	p_a=&a.mat[0][0];
	b=init_matrice(a.lin,a.col,0);
	p_b=&b.mat[0][0];
	b.col=a.col;
	b.lin=a.lin;
	for (i=0; i<a.lin*a.col; i++)
		*p_b++ = *p_a++;
	return b;
}
/*-----------------------------------------*/
long_matrice copy(long_matrice a)
{
	long_matrice b;
	long *p_a,*p_b;
	long i;
	p_a=&a.mat[0][0];
	b=init_long_mat(a.lin,a.col,0);
	p_b=&b.mat[0][0];
	b.col=a.col;
	b.lin=a.lin;
	for (i=0; i<a.lin*a.col; i++)
		*p_b++ = *p_a++;
	return b;
}
/*-----------------------------------------*/
/* cette fonction présume que la puissance du signal en entrée est 
egale à un */
/*
/*This model adds discrete-time zero-mean white gaussian noise of variance
   FACTOR * 10**(-SNR/10) * (Signal_power)
to the input signal, where Signal_ power is 1.
Note:
The input signal of AWGN must be real,
	FACTOR = Ns/(2*K)
	SNR    = 10*log(Eb/N0)
where Ns   = number of samples per symbol
	K    = T/Tb = (transmitted symbol duration) / (information bit duration)
	Eb   = input signal energy per bit of the original
			information-carrying bit stream of rate 1/Tb
	N0/2 = power spectral density of the continuous-
			time Gaussian noise which is assumed to have
			a doublesided bandwidth equal to the sampling rate
*/
void awgn(vector a, double SNR, double factor)
{
	double var;
	long i;
	var=factor*pow(10.0,-SNR/10.0);
	for (i=0 ; i<a.taille ; i++)
		a.vect[i] += n_gaussien(0,var);
}

/*-----------------------------------------*/
/* cette fonction présume que la puissance du signal en entrée est 
egale à un */
/*
/*This model adds discrete-time zero-mean white gaussian noise of variance
   FACTOR * 10**(-SNR/10) * (Signal_power)
to the input signal, where Signal_ power is 1.
Note:
The input signal of AWGN must be real,
	FACTOR = Ns/2K
	SNR    = 10*log(Eb/N0)
where Ns   = number of samples per symbol
	K    = T/Tb = (transmitted symbol duration) / (information bit duration)
	Eb   = input signal energy per bit of the original
			information-carrying bit stream of rate 1/Tb
	N0/2 = power spectral density of the continuous-
			time Gaussian noise which is assumed to have
			a doublesided bandwidth equal to the sampling rate
*/

void awgn(cmplx_vect a, double SNR, double factor)
{
	double var;
	long i;
	var=factor*pow(10.0,-SNR/10.0);
	for (i=0 ; i<a.taille ; i++){
		a.vect[i].Re += n_gaussien(0,var);
		a.vect[i].Im += n_gaussien(0,var);
	}
}

/*----------------------------------------*/
vector BPSK_MOD(long_vector a)
{
	vector b;
	long i;
	b=init_vect(a.taille,0);
	for (i=0 ; i<a.taille ; i++)
		b.vect[i]=((a.vect[i]==1) ? 1.0 : -1.0);
	return b;
}
/*----------------------------------------*/
long_vector	BPSK_DEM(vector a)
{
	long_vector b;
   long i;
	b=init_long_vect(a.taille,0);
	for (i=0 ; i<a.taille ; i++)
		b.vect[i]=((a.vect[i]>0) ? 1 : 0);
	return b;
}
/*----------------------------------------*/
long	comparer(long_vector a, long_vector b)
{
	long cmpt=0,i;
	long *p_a, *p_b;
	p_a=a.vect;
	p_b=b.vect;
	for (i=0 ; i<min(a.taille, b.taille) ; i++)
		if (*p_a++ != *p_b++) cmpt++;
	return cmpt;
}
