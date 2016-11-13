#include "globaldef.h"

#ifdef _NSP_INTEL_H

#define nsp_UsesConversion
#define nsp_UsesVector

#include "nsp.h"
#include "nspmatrice.h"
#include "nspdcplx.h"

// input : 2 echantillons a transmettre (4 Ts)
// output : 2 echantillons sur 2 antennes (8 Ts) 
//			matrice echantillon * antenne

ZMatrix G2(DCplx *ech){
	ZMatrix n(2, 2);

	n.mat[0][0]=ech[0];
	n.mat[0][1]=ech[1];

	n.mat[1][0].re=-ech[1].re;
	n.mat[1][0].im=ech[1].im;
	n.mat[1][1].re=ech[0].re;
	n.mat[1][1].im=-ech[0].im;
	return n;
}

ZMatrix G2(DCplx *ech, int nb_ech){
	if (nb_ech%2!=0) return ZMatrix();

	ZMatrix n(nb_ech, 2);
	int j;
	for(int i = 0;i<nb_ech;i+=2) {
		j=i+1;
		n.mat[i][0]=ech[i];
		n.mat[i][1]=ech[j];
		
		n.mat[j][0].re=-ech[j].re;
		n.mat[j][0].im=ech[j].im;
		n.mat[j][1].re=ech[i].re;
		n.mat[j][1].im=-ech[i].im;
	}

	return n;
}

// input : 4 echantillons a transmettre (4 Ts)
// output : 8 echantillons sur 3 antennes (8 Ts) 
//			matrice echantillon * antenne

ZMatrix G3(DCplx *ech){
	ZMatrix n(8, 3);
	DCplx *e_c=nspzMalloc(4);
	DCplx *n_e=nspzMalloc(4);
	DCplx *n_e_c=nspzMalloc(4);
	nspzbConj2(ech, e_c, 4);
	for(int i =0;i<4;i++) { n_e[i].re=-ech[i].re; n_e[i].im=-ech[i].im; }
	for(i =0;i<4;i++) { n_e_c[i].re=-e_c[i].re; n_e_c[i].im=-e_c[i].im; }

	n.mat[0][0]=ech[0];
	n.mat[0][1]=ech[1];
	n.mat[0][2]=ech[2];

	n.mat[1][0]=n_e[1];
	n.mat[1][1]=ech[0];
	n.mat[1][2]=n_e[3];

	n.mat[2][0]=n_e[2];
	n.mat[2][1]=ech[3];
	n.mat[2][2]=ech[0];

	n.mat[3][0]=n_e[3];
	n.mat[3][1]=n_e[2];
	n.mat[3][2]=ech[1];

	n.mat[4][0]=e_c[0];
	n.mat[4][1]=e_c[1];
	n.mat[4][2]=e_c[2];

	n.mat[5][0]=n_e_c[1];
	n.mat[5][1]=e_c[0];
	n.mat[5][2]=n_e_c[3];

	n.mat[6][0]=n_e_c[2];
	n.mat[6][1]=e_c[3];
	n.mat[6][2]=e_c[0];

	n.mat[7][0]=n_e_c[3];
	n.mat[7][1]=n_e_c[2];
	n.mat[7][2]=e_c[1];

	nspFree(n_e);
	nspFree(n_e_c);
	nspFree(e_c);


	return n;
}


ZMatrix G3(DCplx *ech, int nb_ech){
	if(nb_ech%4!=0) return ZMatrix();

	int nb_ech_2=2*nb_ech;
	ZMatrix n(nb_ech_2, 3);
	DCplx *n_e=nspzMalloc(nb_ech);
	for(int i =0;i<nb_ech;i++) { n_e[i].re=-ech[i].re; n_e[i].im=-ech[i].im; }
	int j;
	
	for(i=0, j=0;i<nb_ech_2;i+=8, j+=4) {
		n.mat[i][0]=ech[j];
		n.mat[i][1]=ech[j+1];
		n.mat[i][2]=ech[j+2];

		n.mat[i+1][0]=n_e[j+1];
		n.mat[i+1][1]=ech[j];
		n.mat[i+1][2]=n_e[j+3];

		n.mat[i+2][0]=n_e[j+2];
		n.mat[i+2][1]=ech[j+3];
		n.mat[i+2][2]=ech[j];

		n.mat[i+3][0]=n_e[j+3];
		n.mat[i+3][1]=n_e[j+2];
		n.mat[i+3][2]=ech[j+1];

		nspzbConj2(n.mat[i], n.mat[i+4], 3);
		nspzbConj2(n.mat[i+1], n.mat[i+5], 3);
		nspzbConj2(n.mat[i+2], n.mat[i+6], 3);
		nspzbConj2(n.mat[i+3], n.mat[i+7], 3);
	}

	nspFree(n_e);

	return n;
}

ZMatrix G4(DCplx *ech){
	ZMatrix n(8, 4);
	DCplx *e_c=nspzMalloc(4);
	DCplx *n_e=nspzMalloc(4);
	DCplx *n_e_c=nspzMalloc(4);
	nspzbConj2(ech, e_c, 4);
	for(int i =0;i<4;i++) { n_e[i].re=-ech[i].re; n_e[i].im=-ech[i].im; }
	for(i =0;i<4;i++) { n_e_c[i].re=-e_c[i].re; n_e_c[i].im=-e_c[i].im; }

	n.mat[0][0]=ech[0];
	n.mat[0][1]=ech[1];
	n.mat[0][2]=ech[2];
	n.mat[0][3]=ech[3];

	n.mat[1][0]=n_e[1];
	n.mat[1][1]=ech[0];
	n.mat[1][2]=n_e[3];
	n.mat[1][3]=ech[2];

	n.mat[2][0]=n_e[2];
	n.mat[2][1]=ech[3];
	n.mat[2][2]=ech[0];
	n.mat[2][3]=n_e[1];

	n.mat[3][0]=n_e[3];
	n.mat[3][1]=n_e[2];
	n.mat[3][2]=ech[1];
	n.mat[3][3]=ech[0];

	n.mat[4][0]=e_c[0];
	n.mat[4][1]=e_c[1];
	n.mat[4][2]=e_c[2];
	n.mat[4][3]=e_c[3];

	n.mat[5][0]=n_e_c[1];
	n.mat[5][1]=e_c[0];
	n.mat[5][2]=n_e_c[3];
	n.mat[5][3]=e_c[2];

	n.mat[6][0]=n_e_c[2];
	n.mat[6][1]=e_c[3];
	n.mat[6][2]=e_c[0];
	n.mat[6][3]=n_e_c[1];

	n.mat[7][0]=n_e_c[3];
	n.mat[7][1]=n_e_c[2];
	n.mat[7][2]=e_c[1];
	n.mat[7][3]=e_c[0];

	nspFree(n_e);
	nspFree(n_e_c);
	nspFree(e_c);

	return n;
}


ZMatrix G4(DCplx *ech, int nb_ech){
	if(nb_ech%4!=0) return ZMatrix();

	int nb_ech_2=2*nb_ech;
	ZMatrix n(nb_ech_2, 4);
	DCplx *n_e=nspzMalloc(nb_ech);
	for(int i =0;i<nb_ech;i++) { n_e[i].re=-ech[i].re; n_e[i].im=-ech[i].im; }
	int j;
	
	for(i=0, j=0;i<nb_ech_2;i+=8, j+=4) {
		n.mat[i][0]=ech[j];
		n.mat[i][1]=ech[j+1];
		n.mat[i][2]=ech[j+2];
		n.mat[i][3]=ech[j+3];

		n.mat[i+1][0]=n_e[j+1];
		n.mat[i+1][1]=ech[j];
		n.mat[i+1][2]=n_e[j+3];
		n.mat[i+1][3]=ech[j+2];

		n.mat[i+2][0]=n_e[j+2];
		n.mat[i+2][1]=ech[j+3];
		n.mat[i+2][2]=ech[j];
		n.mat[i+2][3]=n_e[j+1];

		n.mat[i+3][0]=n_e[j+3];
		n.mat[i+3][1]=n_e[j+2];
		n.mat[i+3][2]=ech[j+1];
		n.mat[i+3][3]=ech[j];

		nspzbConj2(n.mat[i], n.mat[i+4], 4);
		nspzbConj2(n.mat[i+1], n.mat[i+5], 4);
		nspzbConj2(n.mat[i+2], n.mat[i+6], 4);
		nspzbConj2(n.mat[i+3], n.mat[i+7], 4);
	}

	nspFree(n_e);

	return n;
}
#else 
// input : 2 echantillons a transmettre (4 Ts)
// output : 2 echantillons sur 2 antennes (8 Ts) 
//			matrice echantillon * antenne
#include "utilitis.h"

ZMatrix G2(const ZVector &ech){
	ZMatrix n(2, 2);

	n.mat[0][0]=ech.vect[0];
	n.mat[0][1]=ech.vect[1];

	n.mat[1][0]=-ech.vect[1].conj();
	//n.mat[1][0].im=ech.vect[1].im;
	n.mat[1][1]=ech.vect[0].conj();
	//n.mat[1][1].im=-ech.vect[0].im;
	return n;
}

ZMatrix G2(const ZVector &ech, int nb_ech){
	if (nb_ech%2!=0) return ZMatrix();

	ZMatrix n(nb_ech, 2);
	int j;
	for(int i = 0;i<nb_ech;i+=2) {
		j=i+1;
		n.mat[i][0]=ech.vect[i];
		n.mat[i][1]=ech.vect[j];
		
		n.mat[j][0]=-ech.vect[j].conj();
		n.mat[j][1]=ech.vect[i].conj();
	}

	return n;
}

// input : 4 echantillons a transmettre (4 Ts)
// output : 8 echantillons sur 3 antennes (8 Ts) 
//			matrice echantillon * antenne

ZMatrix G3(const ZVector &ech){
	if (ech.taille!=4) return ZMatrix();
	ZMatrix n(8, 3);

	n.mat[0][0]=ech.vect[0];
	n.mat[0][1]=ech.vect[1];
	n.mat[0][2]=ech.vect[2];

	n.mat[1][0]=-ech.vect[1];
	n.mat[1][1]=ech.vect[0];
	n.mat[1][2]=-ech.vect[3];

	n.mat[2][0]=-ech.vect[2];
	n.mat[2][1]=ech.vect[3];
	n.mat[2][2]=ech.vect[0];

	n.mat[3][0]=-ech.vect[3];
	n.mat[3][1]=-ech.vect[2];
	n.mat[3][2]=ech.vect[1];

	n.mat[4][0]=ech.vect[0].conj();
	n.mat[4][1]=ech.vect[1].conj();
	n.mat[4][2]=ech.vect[2].conj();

	n.mat[5][0]=-ech.vect[1].conj();
	n.mat[5][1]=ech.vect[0].conj();
	n.mat[5][2]=-ech.vect[3].conj();

	n.mat[6][0]=-ech.vect[2].conj();
	n.mat[6][1]=ech.vect[3].conj();
	n.mat[6][2]=ech.vect[0].conj();

	n.mat[7][0]=-ech.vect[3].conj();
	n.mat[7][1]=-ech.vect[2].conj();
	n.mat[7][2]=ech.vect[1].conj();

	return n;
}


ZMatrix G3(const ZVector &ech, int nb_ech){
	if (ech.taille!=4) return ZMatrix();
	if(nb_ech%4!=0) return ZMatrix();

	int nb_ech_2=2*nb_ech;
	ZMatrix n(nb_ech_2, 3);
	
	for(register int i=0, j=0;i<nb_ech_2;i+=8, j+=4) {
		n.mat[i][0]=ech.vect[j];
		n.mat[i][1]=ech.vect[j+1];
		n.mat[i][2]=ech.vect[j+2];

		n.mat[i+1][0]=-ech.vect[j+1];
		n.mat[i+1][1]=ech.vect[j];
		n.mat[i+1][2]=-ech.vect[j+3];

		n.mat[i+2][0]=-ech.vect[j+2];
		n.mat[i+2][1]=ech.vect[j+3];
		n.mat[i+2][2]=ech.vect[j];

		n.mat[i+3][0]=-ech.vect[j+3];
		n.mat[i+3][1]=-ech.vect[j+2];
		n.mat[i+3][2]=ech.vect[j+1];

		ZVector tmp, res;
		tmp.taille=3;
		for(register int k=0, l=4; k<4;k++,l++) {
			tmp.vect=n.mat[k];
			res = tmp.conj();
			n.mat[l]=res.vect;
			tmp.vect=res.vect=NULL; 
		}

	}
	return n;
}

ZMatrix G4(const ZVector &ech){
	if (ech.taille!=4) return ZMatrix();

	ZMatrix n(8, 4);

	n.mat[0][0]=ech.vect[0];
	n.mat[0][1]=ech.vect[1];
	n.mat[0][2]=ech.vect[2];
	n.mat[0][3]=ech.vect[3];

	n.mat[1][0]=-ech.vect[1];
	n.mat[1][1]=ech.vect[0];
	n.mat[1][2]=-ech.vect[3];
	n.mat[1][3]=ech.vect[2];

	n.mat[2][0]=-ech.vect[2];
	n.mat[2][1]=ech.vect[3];
	n.mat[2][2]=ech.vect[0];
	n.mat[2][3]=-ech.vect[1];

	n.mat[3][0]=-ech.vect[3];
	n.mat[3][1]=-ech.vect[2];
	n.mat[3][2]=ech.vect[1];
	n.mat[3][3]=ech.vect[0];

	n.mat[4][0]=ech.vect[0].conj();
	n.mat[4][1]=ech.vect[1].conj();
	n.mat[4][2]=ech.vect[2].conj();
	n.mat[4][3]=ech.vect[3].conj();

	n.mat[5][0]=-ech.vect[1].conj();
	n.mat[5][1]=ech.vect[0].conj();
	n.mat[5][2]=-ech.vect[3].conj();
	n.mat[5][3]=ech.vect[2].conj();

	n.mat[6][0]=-ech.vect[2].conj();
	n.mat[6][1]=ech.vect[3].conj();
	n.mat[6][2]=ech.vect[0].conj();
	n.mat[6][3]=-ech.vect[1].conj();

	n.mat[7][0]=-ech.vect[3].conj();
	n.mat[7][1]=-ech.vect[2].conj();
	n.mat[7][2]=ech.vect[1].conj();
	n.mat[7][3]=ech.vect[0].conj();

	return n;
}


ZMatrix G4(const ZVector &ech, int nb_ech){
	if(ech.taille!=4) return ZMatrix();

	if(nb_ech%4!=0) return ZMatrix();

	int nb_ech_2=2*nb_ech;
	ZMatrix n(nb_ech_2, 4);
	
	for(register int i=0, j=0;i<nb_ech_2;i+=8, j+=4) {
		n.mat[i][0]=ech.vect[j];
		n.mat[i][1]=ech.vect[j+1];
		n.mat[i][2]=ech.vect[j+2];
		n.mat[i][3]=ech.vect[j+3];

		n.mat[i+1][0]=-ech.vect[j+1];
		n.mat[i+1][1]=ech.vect[j];
		n.mat[i+1][2]=-ech.vect[j+3];
		n.mat[i+1][3]=ech.vect[j+2];

		n.mat[i+2][0]=-ech.vect[j+2];
		n.mat[i+2][1]=ech.vect[j+3];
		n.mat[i+2][2]=ech.vect[j];
		n.mat[i+2][3]=-ech.vect[j+1];

		n.mat[i+3][0]=-ech.vect[j+3];
		n.mat[i+3][1]=-ech.vect[j+2];
		n.mat[i+3][2]=ech.vect[j+1];
		n.mat[i+3][3]=ech.vect[j];

		ZVector tmp, res;
		tmp.taille=4;
		for(register int k=0, l=4; k<4;k++,l++) {
			tmp.vect=n.mat[k];
			res = tmp.conj();
			n.mat[l]=res.vect;
			tmp.vect=res.vect=NULL; 
		}
	}

	return n;
}

#endif