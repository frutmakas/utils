#include "tools/all.h"
#include "stbc/conversion.h"
#include "rand/randgen.h"
#include "modulation/psk.h"
//#include "fft/fft.h"
#include "modulation/ofdm.h"
#include "fading/fadingtdma.h"

#include <iostream.h>
#include <fstream.h>

int main(void) {
	const int OFDM_K=8, STBC_M=1, STBC_N=2, STBC_P=2, L=3, Q=8, EM_ITER=6; 
	int i, j, k, n, kfactor=8, modulation_size=2;
	WVector bit_in(modulation_size*OFDM_K*STBC_P*Q*kfactor);
	TRandBitStatePtr seed;
	bRandBit1(&seed, bit_in);
	//ZVector psk_tx=z8PSKGray_mod(bit_in);
	ZVector psk_tx = z4PSKGray_mod(bit_in);

	UVector<ZVector> ant(STBC_N);
	for(n=0;n<STBC_N;n++) { // for each stream 
		ant.vect[n] = ZVector(k_factor*QPK);
		for(int t=0;t<k_factor;t++ ) { // for each basic block
			int t0 = t*QPK;
			for(int i=n, j=0; j < QPK;i+=STC_N,j++) {
				ant.vect[n].vect[t0+j]= psk_tx.vect[t0+i];
			}
		}
	}

	//generation des pilots
	ZMatrix X_pilot(OFDM_K, STBC_N*OFDM_K);
	UVector<ZVector> spilot(STBC_N);
	for(i=0;i<STBC_N;i++) {
		spilot.vect[i]=ZVector(OFDM_K);
		int s = (2*i+3)%4;
		for(int k=0;k<OFDM_K;k++) {
			switch(k%3) {
				case 0 : spilot.vect[i].vect[k]=psk4gray.sym[(s+2)%4]; break;
				case 1 : spilot.vect[i].vect[k]=psk4gray.sym[(s+1)%4];break;
				case 2 : spilot.vect[i].vect[k]=psk4gray.sym[(s+3)%4]; break;
			}
		}
	}
		
	for(n=0;n<STBC_N;n++) {
		for(int k=0;k<OFDM_K;k++) {
			X_pilot.mat[k][n*OFDM_K+k]=spilot.vect[n].vect[k];
		}
	}

	ofstream logger;
	logger.precision(4);
	logger.open("xlog.txt");
	logger << ant.vect[0] << ant.vect[1];
	logger.close();


/****************************************
	// full data to transmit (for references..)

	UVector<ZVector> real_ant(STBC_N);
	for(i=0;i<N;i++) {
		real_ant.vect[i]=insert_pilot(spilot.vect[i], ant.vect[i], OFDM_K*STBC_P*Q);
	}
********************************************/
	// nous avons ant.vect[0].taille symboles a transmettre
	// h reste constante pdt OFDM_K*STBC_P donc on a ant.vect[0].taille/OFDM_K/STBC_P 
	//   coef de h a generer pour chaque trajet ant_j -> ant_i
	double baudrate=1e6;
	double max_doppler;
	double symbol_time = 1.0e6/baudrate;
	max_doppler=(long)floor(3.3356409519815204957557671447492e-4*baudrate); // v = 100km/h, c=299792458m/s

	int tapsize=ant.vect[0].taille/OFDM_K;
	UMatrix <ZMatrix> taps(STBC_N,STBC_M);
	UMatrix <dRandUniStatePtr> tapseed(STBC_N,STBC_M);
	logger.open("hijlog.txt");
	ZMatrix temp(3,3);
	temp.mat[0][0]=sqrt(0.8);
	temp.mat[1][1]=sqrt(0.15);
	temp.mat[2][2]=sqrt(0.05);

	for(j=0;j<STBC_N;j++) {
		for(i=0;i<STBC_M;i++) {
			dRandUniInit( &(tapseed.mat[j][i]), 311254488+3*i-2*j, 0.0, 1.0);
			taps.mat[j][i] = fadingtaps(max_doppler, symbol_time, L, tapsize, 512, &(tapseed.mat[j][i]));
			logger << j << " , " << i << " => " <<taps.mat[j][i]*temp << (taps.mat[j][i]*temp).mag();
		}
	}
	logger.close();

	int NB_Q=ant.vect[0].taille/(OFDM_K*STBC_P*Q);

	ZMatrix W(STBC_N*OFDM_K, STBC_N*L);
	for(k=0;k<OFDM_K;k++) {
		for(int l=0;l<L;l++){
			DCplx c = DCplx::exp2RI(1, -UTILITIS_M_2_PI*k*l/OFDM_K);
			for(int n=0;n<STBC_N;n++) {
				W.mat[n*OFDM_K+k][n*L+l]=c;
			}
		}
	}

	zRandGausStatePtr zptr;
	double snr_sweep=20;
	double variance_sweep=convertsnrtosigma2(1.0/3.0/OFDM_K, snr_sweep);
	zRandGausInit(&zptr, 0, variance_sweep);
	UVector<ZMatrix> zpilot(STBC_M);
	for(i=0;i<STBC_M;i++) {
		zpilot.vect[i]=ZMatrix(OFDM_K,1);
		for(int z=0;z<OFDM_K;z++) zpilot.vect[i].mat[z][0]=zRandGaus(&zptr); 
	}

	double tappower[]={0.8, 0.15, 0.05};

	ZMatrix sighicov(STBC_N*L,STBC_N*L), sighi(STBC_N*L,STBC_N*L);
	for(i=0;i<STBC_N*L;i++) {
		sighicov.mat[i][i]=1.0/sqrt(tappower[i%3]);
		sighi.mat[i][i]=tappower[i%3];
	}
	
	ZMatrix hic_left = STBC_P*OFDM_K*ZMatrix::eye(STBC_N*L)+variance_sweep*sighicov;
	ZMatrix pivot(STBC_N*L,STBC_N*L);
	gaussj(hic_left, pivot);
	ZVector x_em(kfactor*Q*STBC_P*OFDM_K);
	
	for(int nbq=0;nbq<NB_Q;nbq++) { // on a NB_Q pilots ... 
		UVector<ZMatrix> h(STBC_M);
		for(i=0;i<STBC_M;i++) {
			h.vect[i]=ZMatrix(STBC_N*L,1);
			for(j=0;j<STBC_N;j++) {
				for(int l=0;l<L;l++) {
					h.vect[i].mat[j*L+l][0]=taps.mat[j][i].mat[nbq*Q][l]*sqrt(tappower[l]);
				}
			}
		}
		UVector <ZMatrix> ypilot(STBC_M);

		for(i=0;i<STBC_M;i++) {
			ypilot.vect[i]=X_pilot*W*h.vect[i]+zpilot.vect[i]; // taille = K x 1
		}

		UVector<ZMatrix> h_est(STBC_M);
		ZMatrix wxxw = W.hermit()*X_pilot.hermit()*X_pilot*W;
		ZMatrix pilot_pivot(STBC_N*OFDM_K,STBC_N*OFDM_K);
		gaussj(wxxw,pilot_pivot);
		logger.open("hijlog.txt", ios::app);
		logger<<wxxw;
		for(i=0;i<STBC_M;i++) {
			h_est.vect[i]= wxxw*W.hermit()*X_pilot.hermit()*ypilot.vect[i];
			logger << i<< " => " << h_est.vect[i] << h_est.vect[i].mag();
		}
		logger.close();
		// ici on a deja une premiere estimation de h a partir des pilotes
		// on va traiter les STBC_P*Q block d'ofdm
		for(int q=0;q<Q;q++) {
			// extraction de h a partir de profile de canal genere au debut duprogramme
			for(i=0;i<STBC_M;i++) {
				h.vect[i]=ZMatrix(STBC_N*L,1);
				for(j=0;j<STBC_N;j++) {
					for(int l=0;l<L;l++) {
						h.vect[i].mat[j*L+l][0]=taps.mat[j][i].mat[nbq*Q+q][l];
					}
				}
			}

			ZMatrix X_tx(STBC_P*OFDM_K, STBC_N*OFDM_K);
			logger.open("timelog.txt");
			for(n=0;n<STBC_N;n++) {
				for(int k=0;k<OFDM_K;k++) {
					for(int p=0;p<STBC_P;p++) {
						X_tx.mat[p*OFDM_K+k][n*OFDM_K+k]=ant.vect[n].vect[nbq*OFDM_K*STBC_P*Q+q*STBC_P*OFDM_K+k*STBC_P+p]; //#averif
						logger << nbq*OFDM_K*STBC_P*Q+q*STBC_P*OFDM_K+k*STBC_P+p << " ";
					}
				}
			}
			logger << endl;
			logger.close();

			logger.open("xlog.txt", ios::app);
			logger << nbq << q << X_tx;
			logger.close();

			UVector<ZMatrix> zi(STBC_M);

			for(j=0;j<STBC_M;j++){// generation des bruits pour les data
				zi.vect[j]=ZMatrix(STBC_P*OFDM_K,1);
				for(int k=0;k<zi.vect[j].line;k++) {
					zi.vect[j].mat[k][0]=zRandGaus(&(zptr));
				}
			}
			UVector<ZMatrix> y(STBC_M);
			for(i=0;i<STBC_M;i++) {
				y.vect[i]=X_tx*W*h.vect[i]+zi.vect[i]; // taille = PK x 1
			}
			// on a les symboles recus ... il faut maintenant determiner X(0)
			ZMatrix X_kappa(STBC_P*OFDM_K, STBC_N*OFDM_K), X_kappa_next(STBC_P*OFDM_K, STBC_N*OFDM_K);

			for(k=0;k<OFDM_K;k++) { // estimation de X(0)
				ZMatrix x_min(STBC_N,1), x(STBC_N,1);
				ZMatrix  Wf(STBC_N,STBC_N*L);
				for(int l=0;l<L;l++) {
					DCplx wkl = DCplx::exp2RI(1.0,-UTILITIS_M_2_PI*k*l/OFDM_K);
					for(int n=0;n<STBC_N;n++) {
						Wf.mat[n][n*L+l]=wkl;
					}
				}

				// on a tt les elements necessaires pour calculer Xkappa
				double somme_qi_min=1e20;
				for(int p1=0;p1<8;p1++) {
					for(int p2=0;p2<8;p2++) {
						double somme_qi=0.0;

						DCplx c1=psk8gray.sym[p1].conj(), c2=psk8gray.sym[p2].conj();
						for(int m=0;m<STBC_M;m++) {
							x.mat[0][0]=c1;
							x.mat[1][0]=c2;
							for(int p=0;p<STBC_P;p++) {
								if(p!=0) { // stbc g2 constraint
									x.mat[0][0]=-c2.conj();
									x.mat[1][0]= c1.conj();
								}
								double qi = (y.vect[m].mat[p*OFDM_K+k][0]-(x.hermit()*Wf*h_est.vect[m]).mat[0][0]).mod2();
								somme_qi+=qi;
							}
						}
						if(somme_qi<somme_qi_min){
							somme_qi_min=somme_qi;
							x_min.mat[0][0]=psk8gray.sym[p1];
							x_min.mat[1][0]=psk8gray.sym[p2];
						}
					}
				}
				X_kappa_next.mat[k][k]=x_min.mat[0][0];
				X_kappa_next.mat[k][OFDM_K+k]=x_min.mat[1][0];
				X_kappa_next.mat[OFDM_K+k][k]=-x_min.mat[1][0].conj();
				X_kappa_next.mat[OFDM_K+k][OFDM_K+k]=x_min.mat[0][0].conj();
			}
			X_kappa=X_kappa_next;
			// on a deja generer X0 ... on comence l'iteration EM
			for(int iter=0;iter<EM_ITER;iter++) {
				ZMatrix x(STBC_N,1), x_min(STBC_N,1);
				for(int k=0;k<OFDM_K;k++) {
					double qi, qi_min=1e20;
					// construction de Wf'
					ZMatrix  Wf(STBC_N,STBC_N*L);
					for(int l=0;l<L;l++) {
						DCplx wkl = DCplx::exp2RI(1.0,-UTILITIS_M_2_PI*k*l/OFDM_K);
						for(int n=0;n<STBC_N;n++) {
							Wf.mat[n][n*L+l]=wkl;
						}
					}
					for(int p1=0;p1<8;p1++) {
						for(int p2=0;p2<8;p2++) {

							DCplx c1=psk8gray.sym[p1].conj(), c2=psk8gray.sym[p2].conj();
							double somme_qi_mk=0.0;
							for(int m=0;m<STBC_M;m++) {
								x.mat[0][0]=c1;
								x.mat[1][0]=c2;
								h_est.vect[m]=hic_left*W.hermit()*X_kappa.hermit()*y.vect[m];
								ZMatrix sighic=sighi-OFDM_K*STBC_P*hic_left*sighi;
								sighic=W*sighic*W.hermit();
								for(int p=0;p<STBC_P;p++) {
									switch(p) {
									case 1: 
										x.mat[0][0]=-c2.conj();
										x.mat[1][0]=c1.conj();
										break;
									default: break;
									}
									ZMatrix tmp1 = x.hermit()*Wf*h_est.vect[m];
									DCplx tmp2 = y.vect[m].mat[p*OFDM_K+k][0]-tmp1.mat[0][0];
									double qi_left=tmp2.mod2();
									double qi_right=0.0;//(x.hermit()*sighic*x).mat[0][0].re;
									qi=qi_left+qi_right;
									somme_qi_mk+=qi;
								}
							}
							if(somme_qi_mk<qi_min) {
								qi_min=somme_qi_mk;
								x_min.mat[0][0]=psk8gray.sym[p1];
								x_min.mat[1][0]=psk8gray.sym[p2];
							}
						}
					}
					cout << iter << " , " << k << " .... " << endl;
					// on a le min ici .. 
					X_kappa_next.mat[k][k]=x_min.mat[0][0];
					X_kappa_next.mat[k][k+OFDM_K]=x_min.mat[1][0];
					X_kappa_next.mat[k+OFDM_K][k]=-x_min.mat[1][0].conj();
					X_kappa_next.mat[k+OFDM_K][k+OFDM_K]=x_min.mat[0][0].conj();
				}
				X_kappa=X_kappa_next;
			}
			//cout << X_tx-X_kappa;
			// on a fini les iterations ... il faut extraire les informations de X_kappa
			for(k=0;k<OFDM_K;k++) {
				for(n=0;n<STBC_N;n++) {
					x_em.vect[nbq*OFDM_K*STBC_P*Q+q*STBC_P*OFDM_K+k*STBC_P+n] =X_kappa.mat[k][n*OFDM_K+k];
					
				}
			}
		}
	}
	WVector bin_out=z8PSKGray_demod(x_em);
	cout << bin_out-bit_in;



	return 0;
}
