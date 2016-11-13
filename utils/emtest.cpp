#include <iostream.h>
#include <fstream.h>

#include "tools/utilitis.h"
#include "modulation/psk.h"
#include "rand/randgen.h"
#include "coding/block/ldpc.h"
#include "tools/binmat.h"
#include "fft/fft.h"
#include <time.h>
#include "fading/fadingtdma.h"
#include "tools/tools.h"
#include "stbc/conversion.h"
#include "modulation/ofdm.h"
#include "tools/pilot.h" 

#include  "tools/superclass.h"
#include "globaldef.h"
#include "tools/interleaver.h"
#ifdef WIN_DOS
#include <conio.h>
#endif

#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 

#ifndef jlkdhfkjhf
#pragma message("compile em")



#define test_em
const int EM_ITER_MAX= 1;

int main(void) {
cout.precision(15);

	TRandBitStatePtr in_seed;
	const int STBC_P = 2, OFDM_K=32, STBC_N=2, STBC_M=1, modulation_set_size=8;
	int input_bit_packet_size = (int) STBC_P*OFDM_K*(int)log2(modulation_set_size);
	int LDPC_M=input_bit_packet_size, LDPC_N=LDPC_M*2;
	int k_factor=1;
	int nb_input_packet = k_factor*input_bit_packet_size;

	WLDPCParityCheck ldpc(LDPC_M, LDPC_N);
	ldpc.make_dense_mixed();
	WVector bit_in(input_bit_packet_size); // taille = 3PK
	bRandBit1(&in_seed, bit_in); 
	WVector ldpc_tx = ldpc.encode(bit_in); // taille = 3NPK
	ZVector psk_tx = z8PSK_mod(ldpc_tx); // taille = NPK = 2PK
	
	UVector<ZVector> ant(STBC_N);

	int L=1;
	int i, j, n;

	int PK = STBC_P*OFDM_K;

	// N=2
	ant.vect[0]=ZVector(PK);
	ant.vect[1]=ZVector(PK);
	int t1=0;
	for(i=0;i<PK;i++) { // symbol interleaved .. even symbol in ant 0, odd symbol in ant1
			ant.vect[0].vect[i]=psk_tx.vect[t1++]; // taille = PK
			ant.vect[1].vect[i]=psk_tx.vect[t1++];
	}

//	int stbc_tx_line = stbc_tx.line;
#define debug1
#ifdef debug1
	ofstream os("ant_debug");
	for(int xxx=0, yyy=0;xxx<2*PK;yyy++) {
		os << psk_tx.vect[xxx++] << ";" << ant.vect[0].vect[yyy]<<endl;
		os << psk_tx.vect[xxx++] << ";" << ant.vect[1].vect[yyy]<<endl;
	}
	os.close();
#endif

#ifndef test_em
	psk_tx.~ZVector(); // on a plus beson de symboles psk
#endif
	
	//ofdm ...	// ATTENTION !! on n'a considere ni les pilotes ni les espacement entre les bloc OFDM
	
	UVector<ZVector> ofdm_ant(STBC_N);

	int data_padded[STBC_N];
	for(n=0;n<STBC_N;n++){
		ofdm_ant.vect[n]=ZVector(PK);
		ofdm_ant.vect[n]=ofdm_encode(ant.vect[n], OFDM_K, data_padded[n]);
/*		for(i=0;i<PK;i+=OFDM_K) { // size = PK
			ofdm_ant.vect[n].insert(fft(ant.vect[n].copy(i,i+OFDM_K)),i);
		}*/

	}
	ofstream os1("ofdm_encode.txt");
	for(int tt=0;tt<OFDM_K;tt++) {
		os1 << ant.vect[0].vect[tt]<<";"<<ofdm_ant.vect[0].vect[tt]<<";"<<ant.vect[1].vect[tt]<<";"<<ofdm_ant.vect[1].vect[tt]<<endl;
	}
	os1.close();
	ant.~UVector(); // symboles stbc sur chaque antenne avant ofdm n'ont plus d'utilit .. 

	// marker1 08/10/2002
	// calcule de attenuation pour chaque trajet et frequence
	typedef UVector<ZMatrix> UZMatrixInVector;
	typedef UVector<dRandUniStatePtr> UdRandUniStatePtrInVector;
	UMatrix<dRandUniStatePtr> tapseed(STBC_N, STBC_M);
	UMatrix<ZMatrix> path(STBC_N, STBC_M);
	const double MAX_DOPPLER=333.0;
	int m;

	for(m=0;m<STBC_M;m++) { //M = antenne reception
		for(int n=0; n< STBC_N;n++) { // N = antenne emmission
			dRandUniInit(&tapseed.mat[n][m], 31212145+2*n+17*n);
			// attenuation reste cte pdt 1 mot ofdm ... 
			path.mat[n][m] = fadingtaps(MAX_DOPPLER, 1, L, STBC_P, 500, &tapseed.mat[n][m]);
		}
	}

	// marker2 22/10/2002 ... attenuation calculee ...
#define debug2
#ifdef debug2
	ofstream os2("path_debug.txt");
	for( m=0;m<STBC_M;m++) { //M = antenne reception
		for( n=0; n< STBC_N;n++) { // N = antenne emmission
			os2 << m << " , " << n << " " << path.mat[n][m].mag() <<endl;
		}
	}
	os2.close();
#endif

	zRandGausStatePtr noiseptr;
	zRandGausInit(&noiseptr,0,convertsnrtosigma2(0.16666,7));

	UVector<ZVector> ofdm_rx(STBC_M);

	DCplx noise, noise_dbg;

	ofdm_rx.vect[0]=ZVector(PK);
	for(int t=0;t<STBC_P;t++) {
		register int t0 = t*OFDM_K;
		noise = zRandGaus(&noiseptr);
		noise_dbg = t==0 ? noise : noise_dbg;
		for(register int k=0;k<OFDM_K;k++) {
			ofdm_rx.vect[0].vect[t0+k] = ofdm_ant.vect[0].vect[t0+k]*path.mat[0][0].mat[t][0].mag()
									   + ofdm_ant.vect[1].vect[t0+k]*path.mat[1][0].mat[t][0].mag()+ noise;
		}
	}
	
/*	os2.open("ofdm_rx.txt");
	os2.precision(12);
	for( tt=0;tt<OFDM_K;tt++) {
		os2 << ofdm_rx.vect[0].vect[tt] <<";"<<ofdm_ant.vect[0].vect[tt] <<";"<<path.mat[0][0].mat[0][0].mag()
			<<";"<<ofdm_ant.vect[1].vect[tt]<<";"<<path.mat[1][0].mat[0][0].mag() << ";" << noise_dbg << endl;
	}
	os2.close();
*/
	ofdm_ant.~UVector(); // symboles transmis n'ont plus d'utilit 
	// marker3 09/10/2002 .. yn calcule ... 
	// arrivant ici .. on a attenue tt signal venant de chaque antenne... (ofdm_rx)
	UVector<ZVector> y(STBC_M);
	y.vect[0]=ofdm_decode(ofdm_rx.vect[0], OFDM_K, data_padded[0]);
	//y.vect[0]=ZVector(PK);
/*	for(i=0;i<PK;i+=OFDM_K) { // taille vecteur  = PK
			y.vect[0].insert(fft(ofdm_rx.vect[0].copy(i,i+OFDM_K)),i);
	}
	*/

	ofdm_rx.~UVector(); // on a plus besoin des symboles ofdm attenue ... 

	// marker4 09/10/2002 : allocation memoire pour superclasses semble correcte
	// c'est avec les symboles y qu'il faut travailler .. 
	
	// que la fete commence avec EM .. yeepie .... 


	ZMatrix W(2*OFDM_K,2);
	for(int k1=0;k1<OFDM_K;k1++) { // calcul de W
		//W.mat[k1][0] = W.mat[OFDM_K+k1][1] = DCplx::exp2RI(1.0, -UTILITIS_M_2_PI*k1/OFDM_K);

		W.mat[k1][0] = W.mat[OFDM_K+k1][1] = DCplx::exp2RI(1.0, 0.0);
	}
	ZMatrix W_hermit = W.hermit(); // ... et son hermitien
	cout << "*************** W *******************" << endl;
	cout << W	<< endl;
	cout << "*************** * *******************" << endl;
	
	ZMatrix X_kappa(OFDM_K, 2*OFDM_K ,DCplx(0.0, 0.0));

	//*************************************
	//faudrait initialiser X_Kappa ici ...
	// init X_kappa
	WVector xest(2,0);

	ZVector xcest(2),xmin(2);


	for(int k0=0;k0<OFDM_K;k0++){
		double vmin=1e20;
		for(int p1=0;p1<8;p1++) {
			for(int p2=0;p2<8;p2++) {
				xcest.vect[0]=psk8.sym[xest.vect[0]];
				xcest.vect[1]=psk8.sym[xest.vect[1]];
				DCplx y1 = (xcest.vect[0]*path.mat[0][0].mat[0][0].mag()+xcest.vect[1]*path.mat[1][0].mat[0][0].mag());
				if((y1-y.vect[0].vect[k0]).mag()<vmin) {
					vmin=(y1-y.vect[0].vect[k0]).mag();
					xmin=xcest;
				}
			}
			vector_increment(xest,8);
		}
		X_kappa.mat[k0][k0]=xmin.vect[0];
		X_kappa.mat[k0][k0+OFDM_K]=xmin.vect[1];
	}

	os2.open("xkappa.txt");

	os2 << X_kappa;
	ZMatrix vvv(OFDM_K,OFDM_K*2);
	for(j=0,i=0;i<OFDM_K;i++) {
		vvv.mat[i][i]=psk_tx.vect[j++];
		vvv.mat[i][i+OFDM_K]=psk_tx.vect[j++];
	}
	os2 << vvv;
	os2 << (vvv-X_kappa).mag();
	cout << (vvv-X_kappa).mag();
	os2.close();

//#define extreme_test

#ifdef extreme_test

#pragma message("Extreme initialise test is used")

	X_kappa =vvv;

#endif

	
	//end init X_kappa
	//*************************************

	// on a stbc_tx_line / odfm_k paquets d'OFDM ..  il faut les traiter paquet par paquet .. 

	//UVector<ZVector> sum_positive(OFDM_K), sum_negative(OFDM_K);
	DVector em_llr(2*input_bit_packet_size), ldpc_llr(2*input_bit_packet_size,0.0);//size 6PK

#ifdef test_em
	ZVector em_est(2*PK);
#endif


	ZMatrix X_kappa_next(X_kappa);
	//X_kappa_next=X_kappa;
	double log_px=0.0;//1.0/8.0;
	for(int iter_em=0;iter_em<EM_ITER_MAX;iter_em++) {
		
		unsigned long em_llr_ptr=0, ldpc_llr_ptr=0;
		for(int p=0;p<STBC_P;p++) {

			for(int k=0;k<OFDM_K;k++) {
				cout << "K = " << k << endl;
				DVector sum_positive(2*(int)log2(modulation_set_size));// taille 6
				DVector sum_negative(2*(int)log2(modulation_set_size));// taille 6
				sum_positive.reset();sum_negative.reset();


				ZMatrix W_f_p(2,2*L);
				for(int w=0;w<2;w++) {
					//W_f_p.mat[w][w]=DCplx::exp2RI(1.0, -UTILITIS_M_2_PI*k/OFDM_K);

					W_f_p.mat[w][w]=DCplx::exp2RI(1.0, 0.0);
				}

				ZMatrix x_k(STBC_N,STBC_M, DCplx()); 
				WVector set_x_k(2,0);

				//search argmin q_i_k
				ZMatrix set_x_k_min(STBC_N, STBC_M, DCplx());

				double arg_min = 1e20;
				for(int p1=0;p1<modulation_set_size;p1++) {
					for(int p2=0;p2<modulation_set_size;p2++) {
						x_k.mat[0][0] = psk8.sym[set_x_k.vect[0]].conj();
						x_k.mat[1][0] = psk8.sym[set_x_k.vect[1]].conj();

						ZMatrix y_mat(OFDM_K,1);
						for(int ym=0;ym<OFDM_K;ym++) {
							y_mat.mat[ym][0]=y.vect[0].vect[p*OFDM_K+ym];
						}

						// calcul de h ... ;
						ZMatrix Sigma_hi_cov = ZMatrix::eye(2); // on a qu'un tap et sa puissance moyene = 1 ==> 1/1 = 1 
						ZMatrix hi_left = W_hermit*X_kappa.hermit()*X_kappa*W+Sigma_hi_cov;
						ZMatrix piv(STBC_N,STBC_N);

						

						//cout <<hi_left;
						gaussj(hi_left, piv); 

						//cout <<hi_left;

						//cout << y_mat;
						ZMatrix h = hi_left*W_hermit*X_kappa.hermit()*y_mat;

						// calcul de sigma_h_*

						ZMatrix Sigma_h=ZMatrix::eye(STBC_N);
						ZMatrix X_kappa_hermit = X_kappa.hermit();
						ZMatrix sigma_h_chapeau =Sigma_h-hi_left*W_hermit*X_kappa_hermit*X_kappa*W*Sigma_hi_cov;

						//cout << sigma_h_chapeau;

						// calcul de q[k]
//#define VERBOSE_EM						

#ifdef VERBOSE_EM

						ofstream ox("hidebug.txt");

						ox.precision(17);

						ox<<W << W_hermit << X_kappa << X_kappa.hermit()<<Sigma_hi_cov << sigma_h_chapeau << y_mat;

						ox<< h;

						ox.close();

#endif



						//if (p==1) cout <<k << "->" << h;
						double qi_k_left = (y_mat.mat[k][0]-(x_k.hermit()*W_f_p*h).mat[0][0]).mod2();
						double qi_k_right = (x_k.hermit()*sigma_h_chapeau*x_k).mat[0][0].re;

						double qi_k = qi_k_left + qi_k_right;
						
						if(qi_k < arg_min) {
							arg_min = qi_k;
							set_x_k_min = x_k;

							//cout <<k<<","<< p1 << "," << p2 <<"->"<< qi_k << "," << arg_min << x_k << set_x_k_min << endl;
						}

						if(set_x_k.vect[0]&4) sum_positive.vect[0]+=exp(-qi_k+log_px); else sum_negative.vect[0]+=exp(-qi_k+log_px);
						if(set_x_k.vect[0]&2) sum_positive.vect[1]+=exp(-qi_k+log_px); else sum_negative.vect[1]+=exp(-qi_k+log_px);
						if(set_x_k.vect[0]&1) sum_positive.vect[2]+=exp(-qi_k+log_px); else sum_negative.vect[2]+=exp(-qi_k+log_px);
						if(set_x_k.vect[1]&4) sum_positive.vect[3]+=exp(-qi_k+log_px); else sum_negative.vect[3]+=exp(-qi_k+log_px);
						if(set_x_k.vect[1]&2) sum_positive.vect[4]+=exp(-qi_k+log_px); else sum_negative.vect[4]+=exp(-qi_k+log_px);
						if(set_x_k.vect[1]&1) sum_positive.vect[5]+=exp(-qi_k+log_px); else sum_negative.vect[5]+=exp(-qi_k+log_px);

						vector_increment(set_x_k, modulation_set_size);
						//cout << set_x_k;
					} // balayage tt les x1  possible
				} /// balayage tt les x2 possible
				//fin de balayage de x

				// x_k min in set_x_k_min ... 

				X_kappa_next.mat[k][k] = set_x_k_min.mat[0][0].conj();
				X_kappa_next.mat[k][OFDM_K+k] = set_x_k_min.mat[1][0].conj();

				em_llr.vect[em_llr_ptr++]=log(sum_positive.vect[0]/sum_negative.vect[1])-ldpc_llr.vect[ldpc_llr_ptr++];
				em_llr.vect[em_llr_ptr++]=log(sum_positive.vect[1]/sum_negative.vect[1])-ldpc_llr.vect[ldpc_llr_ptr++];
				em_llr.vect[em_llr_ptr++]=log(sum_positive.vect[2]/sum_negative.vect[2])-ldpc_llr.vect[ldpc_llr_ptr++];
				em_llr.vect[em_llr_ptr++]=log(sum_positive.vect[3]/sum_negative.vect[3])-ldpc_llr.vect[ldpc_llr_ptr++];
				em_llr.vect[em_llr_ptr++]=log(sum_positive.vect[4]/sum_negative.vect[4])-ldpc_llr.vect[ldpc_llr_ptr++];
				em_llr.vect[em_llr_ptr++]=log(sum_positive.vect[5]/sum_negative.vect[5])-ldpc_llr.vect[ldpc_llr_ptr++];
			} // fin de balayage de k

			X_kappa = X_kappa_next;
			cout << "*********** X_KAPPA **************" << endl;

			for(int w=0;w<OFDM_K;w++) 

				cout << X_kappa.mat[w][w] << " " << X_kappa.mat[w][w+OFDM_K] << endl;

			cout << "*********** ******* **************" << endl;


			cout << "*********** emllr**************" << endl;

			os.open("emllr_debug.txt");

			for(unsigned wj=0;wj<em_llr_ptr;wj++) {

				os << em_llr.vect[wj] << endl;

			}

			os.close();

			cout << em_llr;

			cout << "*********** ******* **************" << endl;



			cout << "EM LLR PTR = " << em_llr_ptr << endl;

#ifndef test_em

#define test_em

#endif
#ifdef test_em
			for(int e=0;e<OFDM_K;e++) {
				em_est.vect[2*p*OFDM_K+2*e]=X_kappa.mat[e][e];
				em_est.vect[2*p*OFDM_K+2*e+1]=X_kappa.mat[e][OFDM_K+e];
			}
#endif
		} // fin de balayage p
	} // fin de l'iteration EM

	// fin d'iteration de EM

	ofstream op("pskdebug.txt");

	for(int ppp=0;ppp<psk_tx.taille;ppp++) {

		op << em_est.vect[ppp] << ";" <<psk_tx.vect[ppp] << endl;

	}

	op.close();
	//

	// conversion proba em -> proba ldpc

	// fin conversion proba em -> proba ldpc

//	cout << em_llr << ldpc_llr;


	WVector output = ldpc.decode(em_llr, ldpc_llr);

	int cnt=0;

	ofstream oss("diff.txt");
	for(j=0;j<output.taille;j++) {

		if(output.vect[j]!=bit_in.vect[j]) cnt++;

		oss << output.vect[j] <<";" <<bit_in.vect[j]<<endl;

	}

	oss.close();


	cout << cnt << " errors out of "<<output.taille<<" ("<<(double)cnt/output.taille<<")" <<endl;
	// conversion proba ldpc --> proba em 
	// fin conversion proba ldpc --> proba em 
	cout << "coucouc " << endl;	
	return 0;
}
#else 

#include <stdlib.h>
#include "modulation/ofdm.h"

#ifdef WINDOS
#include <conio.h>
#endif
#pragma message("hello world")

int main(void) {
	WVector vi(512*3*100);

	TRandBitStatePtr seed;

	bRandBit1(&seed, vi);

	WLDPCParityCheck ldpc(100,200);

	ldpc.make_dense_mixed();

	WVector ldpci = ldpc.encode(vi);

	ZVector psk = z8PSK_mod(ldpci);

	int data_pad;

	ZVector otx = ofdm_encode(psk, 64,data_pad)*DCplx(0.498,0.015)+DCplx(0.05,-0.07);

	
	ZVector orx = ofdm_decode(otx,64,data_pad);

	DVector vo = z8PSK_demod(orx).RZtoNRZ_nip()*0.3187;

	DVector llr(vo.taille);

	WVector xo = ldpc.decode(vo, llr);


	int cnt=0;
	for(int i=0;i<xo.taille;i++) 

		if(xo.vect[i]!=vi.vect[i])  ++cnt;

	cout << cnt << " errors .." << endl;
return 1;


}
#endif 
