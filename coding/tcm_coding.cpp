#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools/utilitis.h"
#define max(x,y)  (((x)>(y)) ? (x):(y))

double log_add(double a,double b)// calcul log(exp(a)+exp(b))
{
	double c;
	c=max(a,b)+log(1+exp(-fabs(a-b)));
	return c;
}

long_matrice back_ward_treillis(long_matrice FW_treillis)
{
	long_matrice BW(FW_treillis.line,FW_treillis.col,0);
	int i,j,k;
//	BW=init_long_mat(FW_treillis.line,FW_treillis.col,0);
	for (j=0; j<BW.col ; j++)
		for (k=0; k<BW.line;k++)
			for(i=0; i<BW.line; i++)
				if (FW_treillis.mat[i][j]==k)
					BW.mat[k][j]=i;
	return BW;
}

//----------------------------------------------------------------------------
// calcul des symboles de tail quand on ferme un treillis en TCM
// entrée : etat_fin : l'état dans lequel on se trouve à partir duquel on ferem le treillis
//          trans : matrice de transitions
//

long_vector calcul_tail(int etat_fin, long_matrice trans)
{
	int i,N,L,j,ind_min,min;
	//long_vector tail;
	N=(long)floor(log10(trans.col)/log10(2.0)+0.5); // nombre de bits à l'entrée du codeur
	L=(long)floor(log10(trans.line)/log10(2.0)+0.5); // nombre de registres
  long_vector tail(L-N+1,0);
	for (i=0; i<L-N+1 ; i++){ // L-N+1 est le nombre de symboles à ajouter pour terminer le treillis à l'état zéro
		ind_min=0;
		for (min=trans.mat[etat_fin][0],j=1; j<trans.col; j++)
			if (trans.mat[etat_fin][j]<min){
				min=trans.mat[etat_fin][j];
				ind_min=j;
			}
		tail.vect[i]=ind_min;
		etat_fin=min;
	}
	return tail;
}
			


/*--------------------------------------------------------------------------------------
								codeur en TCM m/(m+1)
L'entrée est un séquence de bits de taille m.L dont on prend par paquets de m. La sortie est
une séquence de symboles da taille L dans l'alphabet {0..2^(m+1)-1}
par contre si le treillis doit être fermé, la taille sera L+(N° de reg)-m+1=
L+Longueur de contreinte -m
Entrées:
	x: séquence d'entrée à coder. ATTENTION: LA TAILLE DE X PEUT CHANGER EN APPELANT CETTE
	                                         FONCTION SI FERMETURE==1
	trans : matrice de transissions
	output : matrice de sorties
	fermeture : si 0 le treillis reste ouvert sinon il est fermé à l'état zéro.
Sortie:
	y : Séquence codée
---------------------------------------------------------------------------------------*/
long_vector tcm_encoder(long_vector &x, long_matrice trans, long_matrice output, int fermeture)
{
	int m,max_state,max_in,i,mL,j,sym_in,state,next_state, LC,pos;
	int *p_y, *p_nx,*p_x;
//	long_vector y,tail,nouv_x;
	max_in=trans.col; // nombre de symboles en entrée du codeur
	m=(long)floor(log10(trans.col)/log10(2.0)+0.5); // nombre de bits à l'entrée du codeur
	// max_in=2^m
	LC=(long)floor(log10(trans.line)/log10(2.0)+0.5)+1; // Longueur de contreinte

  long_vector y(x.taille/m+(fermeture!=0)*(LC-m),0); //pour ajouter eventuellement les bits de tail

	max_state=trans.line;
	mL=x.taille;
	state=0;
	p_y=y.vect;
	for(i=0;i<mL; i+=m){
		sym_in=0;
		for (j=0; j<m; j++) sym_in=(sym_in<<1) | x.vect[i+j];
		next_state=trans.mat[state][sym_in];
		*p_y++=output.mat[state][sym_in];
		state=next_state;
	}
	if (fermeture==0) return y;
	long_vector nouv_x(x.taille+(LC-m)*m,0);
	p_nx=nouv_x.vect; p_x=x.vect;
	for (i=0; i<x.taille; i++) // copy de x dans nouveau x afin d'y ajouter les bits de tail
		*p_nx++=*p_x++;
	long_vector tail(LC-m,0);
	tail=calcul_tail(state,trans);

	for (i=0;i<LC-m; i++){
		sym_in=tail.vect[i];
		next_state=trans.mat[state][sym_in];
		*p_y++=output.mat[state][sym_in]; //calcul des sortie du codeur
		state=next_state;
		pos=1<<(m-1);// un 1 en position MSB sur les m bits en entrée
		for (j=0; j<m; j++){
			*p_nx++=((pos & sym_in)!=0);// ajout des bits du symbole de tail un par un
			pos=pos>>1; // avancer vers LSB
		}
	}
	//liberer(x);
//	liberer (tail);
	x=nouv_x; // on retourne le nouveau vecteur x
	return y; // retourner les symboles codés
}

/*---------------------------------------------------------------------------
*  C'est un modulateur 8PSK qui prend un vecteur de données entre 0 et 7
*  et qui retourne un vecteur complexe contenant la constelation d'un 8PSK
*  Entrée: Un vecteur entier des symboles
*  Sortie: Un vecteur complexe 8PSK
*/

cmplx_vect mapper_8PSK(long_vector x)
{
	int i;
	cmplx_vect y(x.taille);
	double re[8]={1.0,0.7071,0.0,-0.7071,-1.0,-0.7071,0.0,0.7071};
	double im[8]={0.0,0.7071,1.0,0.7071,0.0,-0.7071,-1.0,-0.7071};
	for (i=0; i<x.taille; i++){
		y.vect[i].re=re[x.vect[i]];
		y.vect[i].im=im[x.vect[i]];
	}
	return y;
}

/*-----------------------------------------------------------------------------
*  C'est un mapper qui prend un vecteur entier des symboles à envoyer et
*  qui retourne un vecteur complexe des symboles de canal
*  Entrée: Un vecteur entier des symboles à envoyer entre 0 et M-1
*          un vecteur complexe contenant des symboles de canal
*  Sortie: Un vecteur complexe de symboles de canal
*/
cmplx_vect mapper(long_vector x, cmplx_vect scat)
{
	int i;
	cmplx_vect y(x.taille);
	for (i=0; i<x.taille; i++)
		y.vect[i]=scat.vect[x.vect[i]];
	return y;
}

/*-----------------------------------------------------------------------------
*  Génère un vecteur complexe utilisable pour la fonction "mapper" dans
*  le cas d'une modulation 8PSK
*/

cmplx_vect scat_8PSK()
{
	int i;
	cmplx_vect scat(8);
	double re[8]={1.0,0.7071,0.0,-0.7071,-1.0,-0.7071,0.0,0.7071};
	double im[8]={0.0,0.7071,1.0,0.7071,0.0,-0.7071,-1.0,-0.7071};
	for (i=0; i<8; i++){
		scat.vect[i].re=re[i];
		scat.vect[i].im=im[i];
	}
	return scat;
}

cmplx_vect scat_Shift8PSK()
{
	int i;
	cmplx_vect scat(8);
	double re[8]={0.9232,0.3872,-0.3872,-0.9232,-0.9232,-0.3872,0.3872,0.9232};
	double im[8]={0.3872,0.9232,0.9232,0.3872,-0.3872,-0.9232,-0.9232,-0.3872};
	for (i=0; i<8; i++){
		scat.vect[i].re=re[i];
		scat.vect[i].im=im[i];
	}
	return scat;
}
/*-----------------------------------------------------------------------------
*   Entrelaceur. Entrelacement sur des symboles (chaque m bits est un symbole)
*	Entrée: in: vecteur entier des bits à entrelacer
*	        intrl: vecteur entier d'entrelacement
*           m: nombre de bits dans un symbole
*   Sortie: vecteur entier des bits entrelacés
*/
long_vector symbol_interleaver(long_vector in, long_vector intrl,int m)
{
	int i,j;
	long_vector out(in.taille);
	for (i=0; i<intrl.taille; i++)
		for (j=0; j<m ; j++)
			out.vect[intrl.vect[i]*m+j]=in.vect[i*m+j];
	return out;
}

// matrice interleaver, le résultat se fait diréctement sur le vecteur d'entrée

void sym_interleaver(matrice in, long_vector intrl)
{
	int i,j,col;
	matrice out;
	out=copy(in);
	for (i=0; i<intrl.taille ; i++){
		col=intrl.vect[i];
		for (j=0; j<in.line ; j++)
			in.mat[j][col]=out.mat[j][i];
	}
	//liberer(out);
}

		

/*----------------------------------------------------------------------------
*   Desentrelaceur. Desentrelacement sur des symboles (chaque m bits est un symbole)
*	Entrée: in: vecteur entier des bits à desentrelacer
*	        intrl: vecteur entier d'entrelacement
*           m: nombre de bits dans un symbole
*   Sortie: vecteur entier des bits desentrelacés
*/
long_vector symbol_deinterleaver(long_vector in, long_vector intrl, int m)
{
	int i,j;
	long_vector out(in.taille);
	for (i=0; i<intrl.taille; i++)
		for (j=0; j<m; j++)
			out.vect[m*i]=in.vect[intrl.vect[i]*m+j];
	return out;
}

void sym_deinterleaver(matrice in, long_vector intrl)
{
	int i,j,col;
	matrice out;
	out=copy(in);
	for (i=0; i<intrl.taille; i++){
		col=intrl.vect[i];
		for (j=0; j<in.line;j++)
			in.mat[j][i]=out.mat[j][col];
	}
	//liberer(out);
}

//-----------------------------------------------------------------------------

cmplx_vect symbol_deinterleaver(cmplx_vect in,long_vector intrl,int m)
{
	int i,j;
	cmplx_vect out(in.taille);
	for (i=0; i<intrl.taille; i++)
		for (j=0; j<m; j++)
			out.vect[m*i]=in.vect[intrl.vect[i]*m+j];
	return out;
}

/*----------------------------------------------------------------------------
  a ______o
          \
		    ---- out
  b ______o/

  Selecteur qui prend un symbole du vecteur a et un symbole du vecteur b
  et qui met le résultat dans le vecteur out

*/

void swap_selector(cmplx_vect a, cmplx_vect b, cmplx_vect out)
{
	int i;
	for (i=0; i<a.taille; i+=2){
		out.vect[i]=a.vect[i];
		out.vect[i+1]=b.vect[i+1];
	}
}
/*----------------------------------------------------------------------------
            /o----- a
          /
  r -----o
          \
		    \o----- b

  L'opération inverse de swap-selector
*/
void swap_combiner(cmplx_vect r, cmplx_vect a, cmplx_vect b)
{
	if (r.taille%2 !=0) Err_Message("size of transmitted vector should be even! (ref: swap_combiner)");
	int i;
	cmplx zero(0,0);
	for (i=0; i<a.taille; i+=2){
		a.vect[i]=r.vect[i];
		b.vect[i]=a.vect[i+1]=zero;
		b.vect[i+1]=r.vect[i+1];
	}
}

/*-----------------------------------------------------------------------------
* Entrée:
*     x: vecteur complexe du signal recu
*     trans: matreice des transitions dans le treillis
*     output: matrice des sorties du codeur
*     scat: vecteur complexe représentant la constellation
*     close_switch=1  =>  Vitrebi n'est pas fermé
*     close_switch=0  =>  Viterbi est fermé en etat 0
* Sortie:
*     y: vecteur entier des données décodées
*/
void TCM_decoder_SIHO(cmplx_vect x, long_vector y, long_matrice trans, long_matrice output,
					  cmplx_vect scat, int close_switch)
{
	double	inf,metric,min;
	long		size_total,M,N,t,i,state,symbol,next_state,fin_state,prev_state;
	cmplx	coded_signal;
//	matrice	metric_cum;
//	long_matrice last_state;

	size_total=x.taille+1; //taille de la matrice de viterbi
	M=trans.line; // nombre des etats dans le treillis
	N=trans.col; // nombre de symboles d'information
	inf=1e6;
	matrice metric_cum(M,size_total,inf);
	long_matrice last_state(M,size_total,0);
	metric_cum.mat[0][0]=0.0;
	for (t=0; t<size_total ; t++){
		for (state=0; state<M; state++){
			for (symbol=0; symbol<N; symbol++){
				next_state=trans.mat[state][symbol];
				coded_signal=scat.vect[output.mat[state][symbol]]; //donne le signal attendu correspondant à la branche
//				metric=modul2(cmplx_sub(x.vect[t],coded_signal));
        metric=(x.vect[t]-coded_signal).mod2();
				metric+=metric_cum.mat[state][t];
				if (metric<metric_cum.mat[next_state][t+1]){
					metric_cum.mat[next_state][t+1]=metric;
					last_state.mat[next_state][t+1]=state;
				}
			}
		}
	}

// TRACE BACK
//Trouver le survivant pour le dernier point de treillis
	fin_state=0;
	if (close_switch==1){ //si viterbi n'est pas fermé on selection le survivant
		min=metric_cum.mat[0][size_total-1];
		fin_state=0;
		for (i=1;i<M;i++)
			if (metric_cum.mat[i][size_total-1]<min){
				fin_state=i;
				min=metric_cum.mat[i][size_total-1];
			}
	}
	for (t=size_total-1; t>0 ; t--){
		prev_state=last_state.mat[fin_state][t];
		i=0;
		while(trans.mat[prev_state][i]!=fin_state) i++;
		// maintenant i est égal à la colonne où se trouve l'etat final
		// d'où l'entrée du codeur pour la transition de prev_etat vers fin_etat
		y.vect[t-1]=output.mat[prev_state][i];
		fin_state=prev_state;
	}
}


/*-----------------------------------------------------------------------------
* Entrée:
*     x: vecteur complexe du signal recu
*     app: a priori matrice de taille t*2^k qui donne la probabilité a priori
*          sur chaque symbole d'information pour chaque t
*     trans: matrice des transitions dans le treillis de taille Ms*2^k
*     output: matrice des sorties du codeur de taille Ms*2^k
*     BW_trans: matrice des transitions backward dans le treillis de taille Ms*2^k
*     scat: vecteur complexe représentant la constellation
* Sortie:
*     LLR: matrice des LLR
*/
void TCM_log_MAP_decoder(cmplx_vect x, matrice app, matrice LLR, long_matrice trans,
	 long_matrice BW_trans,long_matrice output,cmplx_vect scat, double sigma2, int fermeture)
{
	matrice alpha,beta;
	long	Ms,N,size_total,etat_dep,etat_arr,t,j,etat,entree,i;
	double inf,som,gama;

  double inv_2_sigma2 = 1.0/(2.0*sigma2);
	size_total=x.taille+1; // taille du treillis
	Ms=trans.line; // nombre d'états dans le treillis
	N=trans.col; // nombre de symboles d'information (2^k)
	inf=1e40;
	//Calcu de la matrice de beta
	if (fermeture==1) {
		beta=matrice(Ms,size_total,-inf);
		beta.mat[0][size_total-1]=0; // on suppose que le treillis est fermé
	}else beta=matrice(Ms,size_total,log(1.0/Ms));
	for (t=size_total-2; t>=0; t--){
		for (etat_dep=0; etat_dep<Ms; etat_dep++){
			for(som=-inf,j=0;j<N;j++){ // j est l'entrée du codeur (symbole d'info)
				etat_arr=trans.mat[etat_dep][j];// from state etat_dep to state trans[etat_de p][etat_arr]
				if(x.vect[t].mod2()==0) gama=app.mat[j][t];
				else gama=app.mat[j][t]-
					(x.vect[t]-scat.vect[output.mat[etat_dep][j]]).mod2()*inv_2_sigma2;
				som=log_add(som,beta.mat[etat_arr][t+1]+gama);
			}
			beta.mat[etat_dep][t]=som;
		}
	}
/*
FILE *fp;
fp=fopen("toto","w");
fprintf(fp,"\nbeta matrix=\n");
for (i=0; i<beta.col;i++){
for (j=0; j<beta.line; j++)
fprintf(fp,"%lf\t",beta.mat[j][i]);
fprintf(fp,"\n");
}
fprintf(fp,"\n-----------------------------\n\n");*/
	//Calcu de la matrice de alpha
	alpha=matrice(Ms,size_total,-inf);
	alpha.mat[0][0]=0;
	for (t=1; t<size_total; t++){
		for (etat_arr=0; etat_arr<Ms; etat_arr++){
			for(som=-inf,j=0;j<N;j++){ // j est l'entrée du codeur (symbole d'info)
				etat_dep=BW_trans.mat[etat_arr][j];// from state BW_trans[etat_arr][j] to state etat_arr
				if(x.vect[t-1].mod2()==0) gama=app.mat[j][t-1];
				else gama=app.mat[j][t-1]-
					(x.vect[t-1]-scat.vect[output.mat[etat_dep][j]]).mod2()*inv_2_sigma2;
				som=log_add(som,alpha.mat[etat_dep][t-1]+gama);				
			}
			alpha.mat[etat_arr][t]=som;
		}
	}

	//Calcul de LLR
	for (t=0;t<x.taille;t++){
		for (entree=0; entree<N; entree++){ // collone entree et ligne i de la matrice trans est l'etat suivant pour l'entrée du codeur égale à entree
			for (som=-inf,etat=0 ; etat<Ms ; etat++){ // pour tous les etats
				etat_arr=trans.mat[etat][entree];
				if(x.vect[t].mod2()==0) gama=app.mat[entree][t];
				else gama=app.mat[entree][t]-
					(x.vect[t]-scat.vect[output.mat[etat][entree]]).mod2()*inv_2_sigma2;
				som=log_add(som,alpha.mat[etat][t]+gama+beta.mat[etat_arr][t+1]);
			}
			LLR.mat[entree][t]=som;
		}
		for (i=N-1; i>=0; i--)
			LLR.mat[i][t]-=LLR.mat[0][t]; //normaliser par rapport à p[ct=0]
	}

	//liberer(beta);
	//liberer(alpha);
}


// TCM Turbo decoding
// référence page 275 Vucetic

void Turbo_TCM_I(cmplx_vect r, matrice app, long_vector bit_dec,long_vector intrl,
				 long_matrice trans, long_matrice BW_trans, long_matrice output,
				 cmplx_vect scat, double sigma2)
{
FILE *fp;
	int Ms,N,m,pos,pos_MSB,i,j,ind_max;
	double max;
//	long_vector sym_dec;
	int *p_bit_dec;
	long_vector sym_dec(r.taille,0);
	Ms=trans.line; // nombre d'états dans le treillis
	N=trans.col; // nombre de symboles d'information (2^k)
	m=(long)floor(log10(N)/log10(2.0)+0.5); // nombre de bits à l'entrée du codeur
	matrice LLR(N,r.taille,0);
	cmplx_vect rec1(r.taille), rec2(r.taille);
/*	rec1=init_cmplx_vect(r.taille);
	rec2=init_cmplx_vect(r.taille);*/
	swap_combiner(r,rec1,rec2);
	TCM_log_MAP_decoder(rec1,app,LLR,trans,BW_trans,output,scat,sigma2,1);// fermé à l'état zéro

	//sub_matrice(LLR,app,app);
  app=LLR-app;
	sym_interleaver(app,intrl);
	TCM_log_MAP_decoder(rec2,app,LLR,trans,BW_trans,output,scat,sigma2,0);//n'est pas férmé

	sym_deinterleaver(app,intrl);
	sym_deinterleaver(LLR,intrl);
	//sub_matrice(LLR,app,app);
  app=LLR-app;

	// symboles décodés
	for (i=0; i<LLR.col; i++){
		max=LLR.mat[0][i];
		ind_max=0;
		for (j=1; j<LLR.line; j++)
			if (LLR.mat[j][i]>max){
				max=LLR.mat[j][i];
				ind_max=j;
			}
		sym_dec.vect[i]=ind_max;
	}
	pos_MSB=1<<(m-1);// un 1 en position MSB sur les m bits en entrée
	p_bit_dec=bit_dec.vect;

	for (i=0;i<sym_dec.taille;i++){
		pos=pos_MSB;
		for (j=0; j<m; j++){
			*p_bit_dec++=((pos & sym_dec.vect[i])!=0);// ajout des bits du symbole de tail un par un
			pos>>=1;
		}
	}

	//liberer(LLR);
	//liberer(rec1);
	//liberer(rec2);
	//liberer(sym_dec);
}
