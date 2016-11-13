/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/02/11 21:48:38 $
 * $Revision: 1.1.2.40 $
 * $Id: ldpc.cpp,v 1.1.2.40 2004/02/11 21:48:38 syed Exp $
 ********************************************************************/
#include <iostream.h>
#include "rand/randgen.h"
#include "tools/utilitis.h"
#include "tools/tools.h"
#include "tools/binmat.h"
#include "queue/fifo.h"
#include "queue/list.h"
#include "coding/block/ldpc.h"
#include <stdlib.h>

#include "globaldef.h"
#ifdef WIN_DOS
#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 
#endif
#ifdef WIN_DOS
#include <conio.h>
#endif
/* CREATE A SPARSE PARITY-CHECK MATRIX.  Of size M by N, stored in H. */
int MAX_ITERATION=30;

WLDPCParityCheck::WLDPCParityCheck(int row, int col, int nb_bit_per_col, make_method method): WBinarySparseMatrix(row,col), G(), G_cols(), G_rows(){
	
	this->line=row; this->col=col; this->nb_bit_per_col=nb_bit_per_col;
	int cb=nb_bit_per_col;

	if (line <=0 || col <= 0 || cb<=0){
		fprintf(stderr, "Problem in matrix size...!");
		return;
	}
	if (cb>line) {
		fprintf(stderr,"Number of checks per bit (%d) is greater than total checks (%d)\n", cb, line);
		return;
	}

	int no4cycle  = (cb==line && col>1) ? 0 : 1;

	dRandUniInit(&seed, (int)floor(rand()*rand())%5413, 0.0, 1.0);

	// start ldpc generation

	int added, uneven, elim4;
	int i, j, k, t;
	BMatrix parity_check(line,col); // temporary matrix...

	WVector u; // int *u;

	/* Create the initial version of the parity check matrix. */

	switch (method) { 
		case Evencol: 
			{  
				for (j = 0; j<parity_check.col; j++) {
					for (k = 0; k<cb; k++) {
						do { 
							i = (int)floor(parity_check.line*dRandUni(&seed)); 
						} while (parity_check.mat[i][j]==1);
						parity_check.Set(i,j); 
						Set(i,j);
					}
				}
				break;
			}

		case Evenboth:
			{	
				u = WVector(cb*parity_check.col, 0); 
				for (k = cb*parity_check.col-1; k>=0; k--) { 
					u.vect[k] = k%parity_check.line; 
				}

				uneven = 0;
				t = 0;

				for (j = 0; j<parity_check.col; j++){ 
					for (k = 0; k<cb; k++){ 
						for (i = t; i<cb*parity_check.col && parity_check.mat[u.vect[i]][j]==1 ; i++) ;

						if (i==cb*parity_check.col){
							uneven += 1;
							do{
								i = (int)floor(parity_check.line*dRandUni(&seed)); 
							} while (parity_check.mat[i][j]==1);
							parity_check.Set(i,j);
							Set(i,j);
						}
						else { 
							do { 
								i = t + (int)floor((cb*parity_check.col-t)*dRandUni(&seed));
							} while (parity_check.mat[u.vect[i]][j]==1);
							parity_check.Set(u.vect[i],j);
							Set(u.vect[i],j);
							u.vect[i] = u.vect[t];
							t += 1;
						}
					} // L49
				} // L48

				if (uneven>0){
					fprintf(stderr,"Had to place %d checks in rows unevenly\n",uneven);
				}

				break;
			}

		default: 
			fprintf(stderr, "Unknown make method");
			break;;
	}

	/* Add extra bits to avoid rows with less than two checks. */
	/* PART 2 ----------------------- */
	// ajouter la redondance de la transmission

	added = 0;
	//int max_bit_per_row = (int)ceil(parity_check.col*cb/parity_check.line);
	for (i = 0; i<parity_check.line; i++) { 
		while (row_index[i].count<2) {
			do {
				j = (int)floor(parity_check.col*dRandUni(&seed));
			} while (parity_check.mat[i][j]==1);
			parity_check.Set(i,j);
			Set(i,j);
			added++;
		}
	}

	if (added>0) { 
		fprintf(stderr, "Added %d extra bit-checks to make row counts at least two\n",added);
	}

	/* PART 3 ------------------------------------------*/
	/* Add extra bits to try to avoid problems with even column counts. */
	if (cb%2==0 && cb<parity_check.line/*M*/ && parity_check.col/*N*/>1 && added<2) { 
		int a;
		for (a = 0; added+a<2; a++) { 
			do {
				i = (int)floor(parity_check.line*dRandUni(&seed));
				j = (int)floor(parity_check.col *dRandUni(&seed)); 
			} while (parity_check.mat[i][j]==1);
			parity_check.Set(i,j);
			Set(i,j);
		}

		fprintf(stderr,	"Added %d extra bit-checks to try to avoid problems from even column counts\n",	a);
	}
	
	//------------
	if (no4cycle){ 
		/*
			on compare les colonnes entre eux et on veut pas que 2 colonnes ait plus d'un "1" en commun ...*/
			
		elim4=0;
		int f, h, rowsumviolation=0;;

		for(t=0;t<10;t++){
			k=0;
			/* on traite toutes les colonnes*/
			for(j=0;j<parity_check.col;j++) {  // all cols  '--'
				/* on cherche les lignes ayant un "1" */
				for(int e=0;e<col_index[j].count/*colsum.vect[j]*/;e++){ // all filled rows in j-th col ' | '
					int rowpos1 = col_index[j].Get(e+1).position;  // e-th row in j-th col
					/* on cherche d'autres colonnes ayant un "1" sur la meme ligne*/
					for(int fx=0;fx</*rowsum.vect*/row_index[rowpos1].count;fx++){					// find ones in other col at same row
						f=row_index[rowpos1].Get(fx+1);
						if (f==j) continue;
						/* on essaie de voir si cette colonne et la colonne j a un autre "1" sur la meme ligne */
						/* on cherche d'abord d'autre "1" sur la deuxieme colonne*/
						for(int g=0;g<col_index[f].count/*colsum.vect[f]*/;g++) { // go through all 1s in f-th col ' | '
							int rowpos2 =col_index[f].Get(g+1).position;
							if (rowpos1==rowpos2) continue;
							/*  maintenant on verifie si les deux colonnes ait plus d'un "1" en commun..*/
							for(int hx=0;hx<row_index[rowpos2].count/*rowsum.vect[rowpos2]*/;hx++) {
								h=row_index[rowpos2].Get(hx+1);
								if (h==j) { 
									/* pas de bol.. ils ont au moins deux "1" en commun.
									   il faut deplacer un de ces quatre "1"
									*/
									do { 
										i = (int)floor(parity_check.line*dRandUni(&seed)); // on cherche la nouvelle position de 1 sur la colonne j
									} while (parity_check.mat[i][j]==1);
									parity_check.Reset(rowpos1,j);
									Unset(rowpos1,j);
									// added by myself
									if(row_index[rowpos1].count<2) {
										int ww=0;
										do {
											ww = (int)floor(parity_check.col*dRandUni(&seed));
										}while (parity_check.mat[rowpos1][ww]);
										parity_check.Set(rowpos1,ww);
										Set(rowpos1,ww);
										rowsumviolation++;
									}
									// end add
									// ver la ligne i sur la meme colonne
									parity_check.Set(i,j);
									Set(i,j);
									elim4++;
									k++;
									goto nextj; // pour ne pas se casser la tete de faire un if partout
								}
							}
						}
						
					}
				}
nextj:; // c'est pas bo mais c'est pas grave ;)
			}
			if(k==0) break;
		}
		if(elim4>0) {
			fprintf(stderr, "Eliminated %d cycles of length 4 by moving checks within colums \n", elim4);
		}
		if(rowsumviolation>0) {
			fprintf(stderr, "Corrected %d row<2 violation during cycles of 4 elimination\n", rowsumviolation);
		}

		if (t==10) { 
		fprintf(stderr, "Couldn't eliminate all cycles of length four in 10 passes\n");
		}
	}
}

WLDPCParityCheck::~WLDPCParityCheck(){ // c'est l'heure de nettoyage 
}

ZVector operator*(const ZVector &z, const WLDPCParityCheck &l){
	if (z.taille!=l.line) return ZVector();
	ZVector nz(l.col);
	register int i, j;
	for(j=0;j<l.col;j++){
		for(i=1;i<=l.col_index[j].count/*colsum.vect[j]*/;i++){
			nz.vect[j]+=z.vect[l.col_index[j].Get(i).position];
		}
	}
	return nz;
}

DVector operator*(const DVector &z, const WLDPCParityCheck &l){
	if (z.taille!=l.line) return DVector();
	DVector nz(l.col);
	register int i, j;
	for(j=1;j<l.col;j++){
		for(i=0;i<=l.col_index[j].count/*colsum.vect[j]*/;i++){
			nz.vect[j]+=z.vect[l.col_index[j].Get(i).position];
		}
	}
	return nz;
}

WVector operator*(const WVector &z, const WLDPCParityCheck &l){
	if (z.taille!=l.line) return WVector();
	WVector nz(l.col);
	register int i, j;
	for(j=0;j<l.col;j++){
		for(i=1;i<=l.col_index[j].count/*colsum.vect[j]*/;i++){
			nz.vect[j]+=z.vect[l.col_index[j].Get(i).position];
		}
	}
	return nz;							    
}

WVector WLDPCParityCheck::operator*(const WVector &c){
	if(c.taille!=col) return WVector();
	WVector check(line);
	register int m, mc,n;
	for(m=0;m<line;m++){
		check.vect[m]=0;
		for(mc=1;mc<=row_index[m].count;mc++){
			n=row_index[m].Get(mc);
			if (c.vect[n]==1) check.vect[m]^=1;
		}
	}
	return check;
}

BMatrix WLDPCParityCheck::make_dense_mixed(make_method method, int verbose){ 
	int *rows_inv;
	int M=line, N=col;
	G_type=method;


	BIMatrix DH, A, A2(M,N), AI(M,M), B(M,N-M);

	DH=get_clmatrix();

	int i, j, k, n, n2, w, k0=0, /*b0,*/ c, c2, R;
	n = DH.line;

	w=DH.mat[0].word_count;
	n2 = DH.col; 

	G_rows = WVector(n);
	G_cols = WVector(n2);
	for (i = 0; i<n; i++) {
		//rows[i] = i;
		G_rows.vect[i]=i;
		G_cols.vect[i]=i;
	}

	for (j = n; j<n2; j++){
		//cols[j] = j;
		G_cols.vect[j]=j;
	}

	R = 0;
	i = 0;

	for (;;)  { 
		while (i<n-R){
			k0 = G_rows.vect[i];
			
			for (j = i; j<n2; j++) { 
				if(DH.mat[G_cols.vect[j]][k0]) break;
			}

			if (j<n2) break;

			R += 1;
			c = G_rows.vect[i];
			swap(G_rows.vect[i],G_rows.vect[n-R]);
		}

		if (i==n-R) break;

		c = G_cols.vect[j];
		swap(G_cols.vect[i], G_cols.vect[j]);

		A2.Set(G_rows.vect[i],c);

		for (j = 0; j<n2; j++) { 
			if (j!=c && DH.mat[j][k0]){ 
				DH.mat[j]^=DH.mat[c];
				A2.mat[j]^=A2.mat[c];
			}
		}
		i++;
	}
	for (j = n-R; j<n; j++)  { 
		A2.mat[G_cols.vect[j]].Clear();
	}

	n=R;
	
	DH=get_clmatrix();

	if (n>0){ 
		fprintf(stderr,"Note: Parity check matrix has %d redundant checks\n",n);
	}

	rows_inv = new int[M];

	for (i = 0; i<M; i++){
		rows_inv[G_rows.vect[i]] = i;
	}

	A=BIMatrix(M,N);
	for (i = 0; i<A.line; i++)  { 
		for (j = 0; j<A2.col; j++){
			if(A2.mat[j][G_rows.vect[i]]) A.Set(i,j);// else A.Reset(i,j);
		}
    }

	for (j = 0; j<A2.col; j++)  {
		for (k = 0; k<A.mat[0].word_count; k++){ 
			A2.mat[j].vect[k]=A.mat[G_cols.vect[j]].vect[k];
		}
		for ( ; k<A2.mat[0].word_count; k++){ 
			A2.mat[j].vect[k] = 0;
		}
	}
	for (j = 0; j<AI.col; j++)  {
		for (k = 0; k<A2.mat[0].word_count; k++){ 
			AI.mat[j].vect[k]=A2.mat[rows_inv[j]].vect[k];
		}
		for ( ; k<AI.mat[0].word_count; k++){ 
			AI.mat[j].vect[k] = 0;
		}
	}

	for (j = 0; j<B.col; j++)  {
		for (k = 0; k<DH.mat[0].word_count; k++){ 
			B.mat[j].vect[k]=DH.mat[G_cols.vect[j+M]].vect[k];
		}
		for ( ; k<B.mat[0].word_count; k++){ 
			B.mat[j].vect[k] = 0;
		}
	}

	/* Form final generator matrix. */
	if (method==Dense) { 
		BIMatrix Gt=BIMatrix(M,N-M);

		for (j = 0; j<Gt.col; j++){ 
			for (i = 0; i<B.line; i++){
  				if (B.mat[j][i]){ 
					for (k = 0; k<Gt.mat[0].word_count; k++){ 
						Gt.mat[j].vect[k]^=AI.mat[i].vect[k];
					}
				}
			}
		}
		G=Gt;
	}
	else if (method==Mixed)	{ 
		G = AI;
	}

  /* Compute and print number of 1s. */

	if(verbose) {
		if (method==Dense){ 
			c = 0;
			for (i = 0; i<M; i++) { 
				for (j = 0; j<N-M; j++){ 
					c += G.mat[i][j];
				}
			}
			fprintf(stderr, "Number of 1s per check in Inv(A) X B is %.1f\n", (double)c/M);
		}

		if (method==Mixed)  { 
			c = 0;
			for (i = 0; i<M; i++)    { 
				for (j = 0; j<M; j++){ 
					c += G.mat[i][j];
				}
			}
			c2 = 0;
			for (i = M; i<N; i++) { 
				c2 += col_index[G_cols.vect[i]].count;
			}
			fprintf(stderr, "Number of 1s per check in Inv(A) is %.1f, in B is %.1f, total is %.1f\n", (double)c/M, (double)c2/M, (double)(c+c2)/M);
		}
	}
	
	return G;
}

WVector WLDPCParityCheck::encode(const WVector &bitsrc){
	WMatrix u, v;

	int N = col;
	int M = line;
	if(G_type==Dense) {
		u = WMatrix(N-M,1); 
		v = WMatrix(M,1); 
	}

	if(G_type==Mixed) {
		u = WMatrix(M,1);
		v = WMatrix(M,1);
	}

	/* Open source file. */

	if (bitsrc.vect==NULL || bitsrc.taille==0) {
		throw CLDPCParityCheckException(INVALID_SOURCE);
	}

	/* Create encoded output */


	WVector sblk(N-M), cblk(N), chks(M,0);

	/* Encode successive blocks. */

	int output_size=N*(int)ceil(bitsrc.taille/(N-M));
	WVector output(output_size,0);

	for(int srcptr=0, outptr=0;srcptr<bitsrc.taille;srcptr+=N-M, outptr+=N) 	{
		register int i;
		for(i=0;i<sblk.taille;i++) sblk.vect[i]=0;
		sblk.insert(bitsrc.copy(srcptr,srcptr+N-M), 0);
		/* Compute encoded block. */
		switch (G_type){ 
			case Dense: { 
				/* Copy source bits to the systematic part of the coded block. */
				register int j;
				for (j = M; j<N; j++) { 
					cblk.vect[G_cols.vect[j]] = sblk.vect[j-M];
				}

				/* Multiply by Inv(A) X B to produce check bits. */

				for (j = M; j<N; j++){ 
					u.mat[j-M][0]=sblk.vect[j-M];
				}
				
				for(i=0;i<v.line;i++) {
					for(j=0;j<v.col;j++){
						v.mat[i][j]=0;
						for(register int k=0;k<u.line;k++){
							v.mat[i][j]^=G.mat[i][k]*u.mat[k][j];
						}
					}
				}


				for (j = 0; j<M; j++){ 
					cblk.vect[G_cols.vect[j]] = v.mat[j][0]; 
				}
				break;
			}

			case Mixed: { 
				register int j;
				for (j = M; j<N; j++){ 
					cblk.vect[G_cols.vect[j]] = sblk.vect[j-M];

					if (sblk.vect[j-M]==1){ 
						for (int e = col_index[G_cols.vect[j]].Get(1).position, ct=1; ct</*colsum.vect*/col_index[G_cols.vect[j]].count;e = col_index[G_cols.vect[j]].Get(++ct).position) { 
							u.mat[e][0]=1^u.mat[e][0];
						}
					}
				}

				/* Multiply by Inv(A) to produce check bits. */
				for(i=0;i<v.line;i++) {
					for(j=0;j<v.col;j++){
						for(register int k=0;k<u.line;k++){
							v.mat[i][j]^=G.mat[i][k]*u.mat[k][j];
						}
					}
				}
				/* Copy check bits to the right places in the coded block. */
				for (j = 0; j<M; j++){ 
					cblk.vect[G_cols.vect[j]] = v.mat[j][0]; //mod2dense_get(v,j,0);
				}

				break;
			}
			default: 
				cerr << "Unhandled type in " << __LINE__ << " in file " << __FILE__ << endl;
		}

		for(i=0;i<line;i++){
			int sum =0;
			for(register int j=1;j<=row_index[i].count/*rowsum.vect[i]*/;j++) {
				int one = row_index[i].Get(j);
				sum+=cblk.vect[one];
			}
			chks.vect[i]=sum%2;
		}
		
		for (i = 0; i<M; i++) 	{ 
			if (chks.vect[i]==1){ 
				fprintf(stderr,"Output block %d is not a code word!  (Fails check %d)\n",srcptr,i);
				throw CLDPCParityCheckException(ENCODING_PARITY_CHECK_FAILED);
			}
		}

		/* Write encoded block to encoded output file. */
		output.insert(cblk, outptr);
	}

	//fprintf(stderr,"Encoded %d blocks, source block size %d, encoded block size %d\n",(int)ceil(bitsrc.taille/(N-M)),N-M,N);

	return output;

}

void WLDPCParityCheck::special(){
	line=10; col = 20;
	col_index = new TSortedList<CSparseMatrixColInformation>[col];
	row_index = new TSortedList<int>[line];
/*
 seed=0 

   0 1 2 3 4  5 6 7 8 9  0 1 2 3 4  5 6 7 8 9 
   ---------- ---------- ---------- ---------

 0|0 1 0 0 0  0 0 0 0 1  0 0 1 1 0  0 0 0 0 0

 1|1 0 0 0 0  0 0 0 1 0  1 0 0 0 0  1 1 0 0 1

 2|0 0 0 0 0  1 0 0 0 0  0 0 0 0 0  0 0 0 0 0

 3|0 0 1 0 0  0 0 1 0 0  0 0 0 0 1  0 0 0 1 0

 4|0 0 0 1 0  0 0 0 0 0  0 0 0 0 0  0 0 0 0 0



 5|0 0 0 0 1  0 0 0 0 0  0 0 0 0 0  0 0 1 0 1

 6|0 0 0 0 0  0 1 0 0 1  0 0 1 0 0  0 0 0 0 0

 7|0 0 1 0 1  0 1 0 0 1  1 1 0 1 1  0 0 0 1 0

 8|1 1 1 1 0  1 1 1 1 0  1 1 0 0 0  1 1 1 0 0

 9|1 1 0 1 1  1 0 1 1 0  0 1 1 1 1  1 1 1 1 1

*/
	int i=0;
	     Set(1,0);Set(8,0);Set(9,0); //0
	i++; Set(0,i);Set(8,i);Set(9,i);
	i++; Set(3,i);Set(7,i);Set(8,i);
	i++; Set(4,i);Set(8,i);Set(9,i);
	i++; Set(5,i);Set(7,i);Set(9,i);
	

	i++; Set(2,i);Set(8,i);Set(9,i); //5
	i++; Set(6,i);Set(7,i);Set(8,i); //6
	i++; Set(3,i);Set(8,i);Set(9,i); //7
	i++; Set(1,i);Set(8,i);Set(9,i); //8
	i++; Set(0,i);Set(6,i);Set(7,i); //9
	

	i++; Set(1,i);Set(7,i);Set(8,i); //10
	i++; Set(7,i);Set(8,i);Set(9,i); //11
	i++; Set(0,i);Set(6,i);Set(9,i); //12
	i++; Set(0,i);Set(7,i);Set(9,i); //13
	i++; Set(3,i);Set(7,i);Set(9,i); //14


	i++; Set(1,i);Set(8,i);Set(9,i); //15
	i++; Set(1,i);Set(8,i);Set(9,i); //16
	i++; Set(5,i);Set(8,i);Set(9,i); //17
	i++; Set(3,i);Set(7,i);Set(9,i); //18
	i++; Set(1,i);Set(5,i);Set(9,i); //5
}

WLDPCParityCheck::WLDPCParityCheck(){
}

#ifdef LDPC_LLR_INCLUDED
#ifdef WIN_DOS
#pragma message("Compiling LDPC decoding with llr return")
#endif
WVector WLDPCParityCheck::decode(const DVector& input, DVector &llr){

#else 
#ifdef WIN_DOS
#pragma message("Compiling LDPC decoding WITHOUT llr return")
#endif

WVector WLDPCParityCheck::decode(const DVector& input){

#endif

	cout <<endl << "input : " << endl << input << endl;
	cout << " ldpc status " << endl << (*this) << endl;
//	ofstream matlab("ldpc.m");
//	matlaboutput(matlab, "input", input, 15);
	int M=line, N=col ,K=N-M;
	WVector output((int)(floor(input.taille/N))*K);
	prev_decode_status = WVector((int)(floor(input.taille/N)));
	WVector c(N);
	register int m, n;
	int mc, nc;
	DVector z(N);


#ifdef LDPC_LLR_INCLUDED
	//DVector z(N);
	llr = DVector(input.taille, 0);

#endif
	for(int b=0, o=0, status_idx=0 ;b<input.taille;b+=N, status_idx++) { // on ignore le dernier bloc incomplet.... 
		DVector y;
		y=input.copy(b,b+N);
		cout << y;
//		matlaboutput(matlab, "y", y, 15);

		n=0;
		DVector zmn_min1(M,99999), zmn_min2(M,99999); 

		WVector  zmn_min1_idx(M,99999), zmn_min2_idx(M,99999);
#ifdef WIN_DOS 
		bool printstate=false;
		if(kbhit()) {
			printstate=true;
		}
#endif
		//zmn_min1.reset(99999);

/* Step  0: initialisation ******************************************/
		for(n=0;n<N;n++){
			for(int nc=1;nc<=col_index[n].count;nc++){ // initialisation de Zmn
				CSparseMatrixColInformation ci = col_index[n].Get(nc);
				double zmintmp =fabs(ci.Zmn=y.vect[n]);
				if (zmn_min1.vect[ci.position]>=zmintmp) {
					if(n!=zmn_min1_idx.vect[ci.position]){  
						zmn_min2.vect[ci.position]=zmn_min1.vect[ci.position];
						zmn_min2_idx.vect[ci.position]=zmn_min1_idx.vect[ci.position]; //downgrade
						zmn_min1.vect[ci.position]=zmintmp;
						zmn_min1_idx.vect[ci.position]=n; //update
					} else {
						zmn_min2.vect[ci.position]=zmintmp;
						zmn_min2_idx.vect[ci.position]=n; //update
					}

				} else if (zmn_min2.vect[ci.position]>zmintmp) {
					zmn_min2.vect[ci.position]=zmintmp;
					zmn_min2_idx.vect[ci.position]=n; //update
					cout << " update !!! " << endl;
				}
				col_index[n].Update(nc,ci);
			}
		}

		cout << zmn_min1 << zmn_min2 << zmn_min1_idx << zmn_min2_idx;
//		matlaboutput(matlab, "zmn_min1", zmn_min1, 15);
//		matlaboutput(matlab, "zmn_min2", zmn_min2, 15);



		WVector Sm(M); // sigma(m);
		CSparseMatrixColInformation tmp;

/* Start iteration process .... *********************************/
		int t;
		for(t=0;t<MAX_ITERATION;t++){ // doing MAX_ITERATION or less
/* Step 1 : Horizontal step (processing in check nodes ..... **********************/	
			Sm.reset();

			// processing sigma(m,n) and sigma(m)
			for(n=0;n<N;n++){ // all col
				for(nc=1;nc<=col_index[n].count;nc++) { // all line in n-th column
					tmp=col_index[n].Get(nc);
					m=tmp.position; // we are dealing with m-th line of n-th col
					if(/*smn=*/tmp.Zmn>0?1:0){ // sigma(m,n)=1 if Zmn>0, or else 0; [STEP 1.a]
						Sm.vect[m]^=1; // Sigma(m)=Sum(sigma(m,n) on col n) modulo 2 [STEP 1.b]
					}
				}
			}
			cout << "Sm : " << Sm << endl;

			DVector SumLm(N);
			SumLm.reset(0);

			for(m=0;m<M;m++){ // [STEP 1.c] : update Lmn
				for(mc=1;mc<=row_index[m].count;mc++){
					int nctmp;
					n=row_index[m].Get(mc);
					tmp=col_index[n].Get(nctmp=col_index[n].isExist(m));
					double zmnp_min=n!=zmn_min1_idx.vect[m]?zmn_min1.vect[m]:zmn_min2.vect[m];
					int SmSmnBar=!(tmp.Zmn>0?Sm.vect[m]^1:Sm.vect[m]);
					tmp.Lmn=SmSmnBar==0?zmnp_min:-zmnp_min; // Lmn=(-1)^(SmSmBar)*min|Zmn'|
					col_index[n].Update(nctmp, tmp);
					SumLm.vect[n]+=tmp.Lmn; // for STEP 2
				}
			}

			cout << "SumLm :"<< SumLm <<endl;

/* STEP 2 : Vertical  step (processing in bit nodes .... ******************/

			zmn_min1.reset(99999.0);
			zmn_min2.reset(99999.0);
			zmn_min1_idx.reset(99999);
			zmn_min2_idx.reset(99999);

			cout << zmn_min1 << zmn_min2 << zmn_min1_idx << zmn_min2_idx;
			
/*
#ifndef LDPC_LLR_INCLUDED

			DVector z(N);
#endif*/

			for(n=0;n<N;n++){
				for(nc=1;nc<=col_index[n].count;nc++){
					tmp=col_index[n].Get(nc);
					int m=tmp.position;
					tmp.Zmn=y.vect[n]+SumLm.vect[n]-tmp.Lmn;
					if (fabs(tmp.Zmn)<=zmn_min1.vect[m]) { // looking for min|zmn'|
						if(n!=zmn_min1.vect[m]) {
							zmn_min2.vect[m]=zmn_min1.vect[m];
							zmn_min2_idx.vect[m]=zmn_min1_idx.vect[m]; // downgrade
							zmn_min1.vect[m]=fabs(tmp.Zmn); 
							zmn_min1_idx.vect[m]=n;//tmp.position;
						} else {
							zmn_min2.vect[m]=fabs(tmp.Zmn); zmn_min2_idx.vect[m]=n;
						}
					} else if(fabs(tmp.Zmn)<zmn_min2.vect[m]) {
						zmn_min2.vect[m]=fabs(tmp.Zmn); zmn_min2_idx.vect[m]=n;
					}
					col_index[n].Update(nc, tmp);
				}
				z.vect[n]=y.vect[n]+SumLm.vect[n];
			}
			cout << zmn_min1 << zmn_min2 << zmn_min1_idx << zmn_min2_idx;

/* STEP 3 : Hard Decision and stopping criterion test ... ********************/
			for(n=0;n<N;n++) 
				c.vect[n]=z.vect[n]>0?1:0;
			cout << " t = " << t << " c = "<< c;
			WVector check;
			check=(*this)*c;
			bool lock;
			lock=false;
			for(m=0;m<M;m++) {
				if(check.vect[m]) { 
					lock=true; 
					break; 
				}
			}
			if(lock) continue; else break; // stopping criterion is valid
		} // t

		if (t>=MAX_ITERATION) {
#ifdef VERBOSE_DECODING
			cerr << endl << "erronous but max iteration reached at " << o << "(max="<<MAX_ITERATION<<")"<<endl;
#else 
#ifdef WIN_DOS
			if (printstate) {
				char junkie=getch();
				cerr << endl << "erronous but max iteration reached at " << o << "(max="<<MAX_ITERATION<<")"<<endl;
			}
#endif // win_dos
#endif // verbose_decoding #1
			prev_decode_status.vect[status_idx]=DECODING_FAILED;
		} else   {
			prev_decode_status.vect[status_idx]=DECODING_SUCCESSFUL;

#ifdef VERBOSE_DECODING
			cerr << endl << t << " iteration and no more errors at" << o << "(max="<<MAX_ITERATION<<")"<<endl;
#else
#ifdef WIN_DOS
			if (printstate) {
				printstate=false;
				cerr << endl << t << " iteration and no more errors at" << o << "(max="<<MAX_ITERATION<<")"<<endl;
			}
#endif// win_dos

#endif // verbose_decoding #2
		}

		for(int shuffle=M;shuffle<N;shuffle++){
			output.vect[o++]=c.vect[G_cols.vect[shuffle]];
		}


#ifdef LDPC_LLR_INCLUDED

		llr.insert(z, b);

#endif
	}
	return output;
}

void WLDPCParityCheck::mackay(){
	cerr << "./make-pchk mackay.txt "<<line << " " << col << " ";
	for(register int m=0;m<line;m++){
		for(register int n=1;n<=row_index[m].count;n++){
			cerr << m << ":" << row_index[m].Get(n) << " " ;
		}
	}
	cerr.flush();
}

ostream &operator<<(ostream &os, WLDPCParityCheck &b){
//output basic parameter
	os << b.nb_bit_per_col << endl << b.no4cycle << endl << b.seed ;
// output H
	//os << (WBinarySparseMatrix)b;
	os << b.line << "," << b.col << endl;
	for(register int m=0;m<b.line;m++){
		os << m << "," << b.row_index[m].count << " : " ;
		for(register int n=1;n<=b.row_index[m].count;n++){
			os << b.row_index[m].Get(n) << " " ;
		}
		os << endl;
	}
// output G
	int gtype=b.G_type ;
	os  << b.G << b.G_cols << b.G_rows<< gtype << endl ;
	return os;
}

ifstream &operator>>(ifstream &is, WLDPCParityCheck &b) {
	char junk;
	int line, col, ijunk, count, elt;
	//input basic parameter 
	is >> b.nb_bit_per_col >> b.no4cycle >> b.seed;
	//input H
	is >> line >> junk >> col;
	b.~WLDPCParityCheck();
	b.col=col; b.line=line;
	b.col_index = new TSortedList<CSparseMatrixColInformation>[col];
	b.row_index = new TSortedList<int>[line];

	for(register int m=0;m<line;m++) {
		is >> ijunk >> junk >> count >> junk;
		for(register int n=0;n<count;n++) {
			is >> elt;
			b.Set(m,elt);
		}
	}
	//input G
	int gtype;
	is  >> b.G >> b.G_cols >> b.G_rows;
	is >>gtype;
	b.G_type=(make_method)gtype; 

	return is;
}











