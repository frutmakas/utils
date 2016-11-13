#include "globaldef.h"

// input : 2 echantillons a transmettre (4 Ts)
// output : 2 echantillons sur 2 antennes (8 Ts) 
//			matrice echantillon * antenne
#include "tools/utilitis.h"
#ifdef WIN_DOS
#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 
#endif

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

ZMatrix G2_type2(const ZVector &ech, int nb_ech){
	if (nb_ech%2!=0) return ZMatrix();

	ZMatrix n(nb_ech, 2);
	int j;
	for(int i = 0;i<nb_ech;i+=2) {
		j=i+1;
		n.mat[i][0]=ech.vect[i];
		n.mat[i][1]=ech.vect[j];
		
		n.mat[j][0]=ech.vect[j].conj();
		n.mat[j][1]=-ech.vect[i].conj();
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
