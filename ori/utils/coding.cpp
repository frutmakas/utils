#include <stdio.h>
#include "utilitis.h"

/*---------------------------------------------------------*/
/* Codeur convolutif systematique
*  les entrees:
*  	matrice du code (RS)sous forme d'un vecteur g (a deux ligne, et
*		une colonne)
*		Le nombre de bits en sortie du codeur: n (n doit etre egal a 2)
*		Longueur de contraint du codeur : k
*		vecteur d'entrees : x
*		param : si 1 demande une terminaison de treillis
*  sortie:
*		un vecteur de la taille x.taille*2 si param=1
*                         et (x.taille+k-1)*2 si param<>1
*     les indices 2*i sont les sorties systématiques
*     les indices 2(i+1) sont les sorties codées
*  copyright:
*         Vahid Meghdadi
*         ENSIL/Universite de Limoges
*         janvier 2000
*         last modification: Avril 2001
*/

long_vector codeur_rsc_2(long_vector x, long_vector g, long k, long param)
{
	long_vector a,code,x_tail;
	long i, j, accu, L_total;
	long *p_x_tail,*p_code;
	x_tail=init_long_vect(k-1,0);
   p_x_tail=x_tail.vect;
   
/* initialiser les etats du registre a decalage egal a zero  */   
   a = init_long_vect(k,0);
 
/* if param = 1 on ajoute des bits a la fin des donnee pour terminer le treillis
*  et mettre l'etat final a zero. Dans ce cas les donnees ajoutees sont sorties
*  dans un vecteur : x_tail */

   if (param == 1) 	L_total=x.taille + k -1;
       else		L_total = x.taille;
	code=init_long_vect(L_total*2,0);
   p_code=code.vect;


/* Pour chaque entree, calculer le bit codé */
   for (j=0;j<L_total ; j++){
	accu = 0;

/* a_k=zigma(g(i)*u_(k-i)) + u(k)  */

	for (i=1 ; i<k ; i++)
		accu += a.vect[i] * ((g.vect[0]>>(k-i-1))&1);
	if (j<x.taille) {
		accu += x.vect[j];
		a.vect[0] = accu & 1;
		*p_code++ = x.vect[j]; /* sortie systématique */
	}
	else {
		*p_x_tail = accu & 1;
		a.vect[0]=0;
		*p_code++ = *p_x_tail++; /* sortie systématique de "tail" */
	}

/* bit code est calcule a partir de deuxieme element de matrice comme un CNS */
	accu = 0;
	for (i=0; i<k ; i++)
		accu += a.vect[i] * ((g.vect[1]>>(k-i-1))&1);
	*p_code++ = accu & 1;

/* Decalage d'etats dans le registre a decalage  */
	for (i=k-1 ; i != 0; i--)
		a.vect[i] = a.vect[i-1];
   }
	liberer (a);
	liberer(x_tail);
	return code;
}

/*------------------------------------------------------------*/
/* Calcul des matrices de backward transition et de backward sortie
*   les entrée:
*		g: vecteur du codeur recursif et systematique à n elements.
*			la premiere ligne est la partie recursive sous forme décimale
*			les ligne suivantes sont les polynomes de calcul de la sortie
*				en décimale
*		k: longueur de contraint (nombre de registre +1
*	les sorties:
*		last_state : en entrée est une variable de type long_vector non initialisé
*			en sortie c'est une matrice (n_state x 2). la premiere colonne
*			est l'etat précedent quand l'entrée est égale à 0.
*			la deuxieme colonne est l'etat precedent quand l'entrée est égale
*			à 1. L'etat actuel est le numéro de ligne.
*		last_out : en entrée est une variable de type long_vector non initialisé
*			en sortie c'est une matrice (n_state x 2). la premiere colonne
*			est la sortie du codeur quand l'entrée était égale à 0 pour ariver à
*			l'etat actuel. La deuxieme colonne est quand l'entrée était 1.
*  copyright:
*         Vahid Meghdadi
*         ENSIL/Universite de Limoges
*         Avril 2001
*/
void backward_treillis_rsc
	(long_vector g, long k, long_matrice &last_state, long_matrice &last_out)
{
	long  m_state, state, i, bit, tmp, m, j,res, a_k;
	m = k-1;
	m_state = 1 << m; /*  2^m  */
	last_state=init_long_mat(m_state,2,0);
	last_out=init_long_mat(m_state,2,0);
	for (state=0 ; state<m_state ; state++)
          for (bit=0 ; bit<2 ; bit++){
               tmp = g.vect[0] & (state<<1);
               for (i=0, res=0 ; i<=m; i++, tmp >>=1)
                    res ^= tmp & 1;
               res = (res ^ bit) | (state <<1);
               a_k=res;
               last_state.mat[state][bit] = res & ~(1 << m);
/* calcul de last_out */
			   last_out.mat[state][bit]=bit;
			   for (j=1; j<g.taille ; j++){
                  tmp = g.vect[j] & a_k;
               	  for (i=0, res=0 ; i<=m; i++, tmp >>=1)
                    	res ^= tmp & 1;
               	  last_out.mat[state][bit] = (last_out.mat[state][bit] <<1) | res;
			   }
          }
}
/*------------------------------------------------------------*/
/* Calcul des matrices de transition et de sortie (forward)
*   les entrée:
*		g: vecteur du codeur recursif et systematique à n elements.
*			la premiere ligne est la partie recursive sous forme décimale
*			les ligne suivantes sont les polynomes de calcul de la sortie
*				en décimale
*		k: longueur de contraint (nombre de registre +1
*	les sorties:
*		next_state : en entrée est une variable de type long_vector non initialisé
*			en sortie c'est une matrice (n_state x 2). la premiere colonne
*			est l'etat suivant quand l'entrée est égale à 0.
*			la deuxieme colonne est l'etat suivant quand l'entrée est égale
*			à 1. L'etat actuel est le numéro de ligne.
*		next_out : en entrée est une variable de type long_vector non initialisé
*			en sortie c'est une matrice (n_state x 2). la premiere colonne
*			est la sortie du codeur correspondant à la branche
*			(N° de ligne) -> (etat suivant) quand l'entrée est égale à 0.
*			La deuxieme colonne est la meme chose quand l'entrée est à 1.
*  copyright:
*         Vahid Meghdadi
*         ENSIL/Universite de Limoges
*         Avril 2001
*/
void forward_treillis_rsc
		(long_vector g,long k, long_matrice &next_state, long_matrice &next_out)
{
	long m,m_state,i,j,bit,tmp,res,a_k,state;
	m=k-1;
	m_state = 1 << m;
	next_state=init_long_mat(m_state,2,0);
	next_out=init_long_mat(m_state,2,0);
	for (state=0 ; state<m_state ; state++)
		for (bit=0 ; bit<2 ; bit++){
			tmp = g.vect[0] & state;
			for (i=0, res=0 ; i<m; i++, tmp >>=1)
                    res ^= tmp & 1;
			res = ((res ^ bit) << m) | state;
			next_state.mat[state][bit]=res>>1;
			a_k=res;
			next_out.mat[state][bit]=bit;
			for (j=1; j<g.taille ; j++){
				tmp = g.vect[j] & a_k;
				for (i=0, res=0 ; i<=m; i++, tmp >>=1)
														res ^= tmp & 1;
				next_out.mat[state][bit] = (next_out.mat[state][bit] <<1) | res;
			}
		}
}

/*-----------------------------------------------------------*/

/* Sous programme de calcul de SOVA
*  Les entrees:
*  - Vecteur du signal recu de la taille n*info_len, les echantillons
*    i*n sont les parties systematiques et les autres sont la partie codee.
*  - info_len: la taille des informations non codee (tail compris)
*  - n : 1/n le rendement du code
*  - k : longueur de contraint [g]n*k
*  - L_a : informations a priori de la taille info_len
*  - last_state : matrice contenant les etats precedant l'etat
*    actuel [n_state*2]
*  - last_out : matrice contenant les bits codee precedents sachant
*    l'etat actuel et le bit d'information
*  - ind_cod : indice de codeur: 1 pour un decodeur bien termine a l'etat 0
*                                2 pour un code ouvert.
*  Copyright Vahid Meghdadi
*  ENSIL / Universite de Limiges
*  Janvier 2000
*  last modification: Avril 2001
*/

vector sova_paquet(vector sig_rec, long info_len, long n, long k, vector L_a,
			 long_matrice last_state, long_matrice last_out, long ind_dec)
{
	long	t, state, n_state, sym0, sym1, state0, state1, i,bit, temp, j, delta;
	double     infini, mc0, mc1, LRV;
	double *p_L_a;
	vector L_uk;
	matrice m_c,delta_m;
	long_matrice bit_prec;
	long_vector bit_est, survivant;
	infini = 1e5;
	n_state = 1 << (k-1); /* 2^(k-1)  */
/* matrice des metriques cumulées initialisée à moins infini */
	m_c = init_matrice(n_state,info_len+1,-infini);


	delta_m = init_matrice(n_state,info_len+1,0);
	bit_prec = init_long_mat(n_state,info_len+1,0);
	survivant = init_long_vect(info_len+1,0);
	bit_est = init_long_vect(info_len,0);

	delta = k*5; /*la profondeur de calcul de revision */
	m_c.mat[0][0]=0;
	p_L_a=L_a.vect;
	for (t=0 ; t<info_len ; t++){
		for (state=0 ; state<n_state ; state++){
			sym0 = last_out.mat[state][0];
			sym1 = last_out.mat[state][1];
			state0 = last_state.mat[state][0];
			state1 = last_state.mat[state][1];
			for (i=0, mc0=0, mc1=0 ; i<n ; sym0>>=1,sym1>>=1,i++){
				mc0 += ((sym0 & 1) != 0) ? sig_rec.vect[t*n+n-1-i]
								: -sig_rec.vect[t*n+n-1-i];
				mc1 += ((sym1 & 1) != 0) ? sig_rec.vect[t*n+n-1-i]
								: -sig_rec.vect[t*n+n-1-i];
			}
			mc0 += m_c.mat[state0][t] - *p_L_a/2.0;
			mc1 += m_c.mat[state1][t] + *p_L_a/2.0;
			if (mc0 > mc1) {
				m_c.mat[state][t+1] = mc0;
				delta_m.mat[state][t+1] = mc0 - mc1;
				bit_prec.mat[state][t+1] = 0;
			} else {
				m_c.mat[state][t+1] = mc1;
				delta_m.mat[state][t+1] = mc1 - mc0;
				bit_prec.mat[state][t+1] = 1;
			}
		}
		p_L_a++;
	}

/* Trace back
*  pour le decodeur No 1 (ind_dec=1) on commence de l'etat zero car
*  le treillis est corectement termine
*  pour le decodeur No 2 (ind_dec=2) on commence par l'etat qui donne
*  une metrique maximale  */

	if (ind_dec==1) survivant.vect[info_len] = 0;
	else {	vector _a;
				_a=col_mat2vect(m_c,info_len);
				survivant.vect[info_len] = indice_max(_a);
				liberer(_a);
	}
	for (t=info_len-1 ; t>=0 ; t--){
		bit_est.vect[t]=bit_prec.mat[survivant.vect[t+1]][t+1];
		survivant.vect[t]=last_state.mat[survivant.vect[t+1]][bit_est.vect[t]];	}

/* Chercher les chemins concurants qui donnent une decision sur le bit
*  detecté differente  */
	L_uk=init_vect(info_len,0);
	for (t=0 ; t<info_len ; t++){
		LRV = infini;
		for (i=0 ; i<delta ; i++)
			if (t+i<info_len){
				bit = 1 - bit_est.vect[t+i];
				temp = last_state.mat[survivant.vect[t+i+1]][bit];
				for (j=i-1 ; j>=0 ; j--){
					bit = bit_prec.mat[temp][(t+j+1)];
					temp = last_state.mat[temp][bit];
				}
				if (bit != bit_est.vect[t])
					if (delta_m.mat[survivant.vect[t+i+1]][t+i+1]<LRV)
						LRV = delta_m.mat[survivant.vect[t+i+1]][t+i+1];
			}
		L_uk.vect[t]=(2*bit_est.vect[t] - 1) * LRV;
	}
	liberer(m_c);
	liberer(delta_m);
	liberer(bit_prec);
	liberer(survivant);
	liberer(bit_est);
	return L_uk;
}
/*----------------------------------------------------------*/
vector matrix_interleaver(vector a,long lin, long col, long param)
{
	long i,j,temp;
	vector out;
	out=init_vect(a.taille,0);
	if (param==2) {  /* deinterleaver */
		temp=lin; lin=col; col=temp;}
	for (i=0;i<lin ; i++)
		for (j=0 ; j<col ; j++)
			out.vect[j*lin+i]=a.vect[i*col+j];
	return out;
}
/*----------------------------------------------------------*/
long_vector matrix_interleaver(long_vector a,long lin, long col, long param)
{
	long i,j,temp;
	long_vector out;
	out=init_long_vect(a.taille,0);
	if (param==2) {  /* deinterleaver */
		temp=lin; lin=col; col=temp;}
	for (i=0;i<lin ; i++)
		for (j=0 ; j<col ; j++)
			out.vect[j*lin+i]=a.vect[i*col+j];
	return out;
}


