#include "all.h"
#include "coding/block/ldpc.h"
#include <iostream.h>
#include <fstream.h>
#include "tools/ext_matrix.h"

#pragma message(__TIMESTAMP__ " : Compiling " __FILE__ )
#define ok
#ifndef ok
//#define _MATLAB_DEBUG_
#ifdef _MATLAB_DEBUG_
void writematlab(double SNR, ZVector psk_tx, ZMatrix sigma_h, ZMatrix sigma_h_cross, UMatrix<ZMatrix> path, ZMatrix h, UVector<ZMatrix> noise, ZMatrix W, ZMatrix X_pilot) {
	ofstream m("prepare.m");
	m << "%[snr, psk_tx, ant1, ant2, sigma_h, sigma_h_cross, path_1_1, path_2_1, h_curr, noise, W, X_pilot] = prepare() " << endl;
	m << "function [snr, psk_tx, ant1, ant2, sigma_h, sigma_h_cross, path_1_1, path_2_1, h_curr, noise, W, X_pilot] = prepare() " << endl<< endl;
	int i, j;
	m << "snr = " << SNR<< endl ;
	m << "psk_tx = [";
	for(i=0;i<psk_tx.taille;i++) 
		m << psk_tx.vect[i] << " ";	
	m << "];" << endl << endl;

	m << "[m,n]=size(psk_tx);\n[ant1, ant2]=split(psk_tx, n);\n" << endl;

	m << "X_pilot = [";
	for(j=0;j<X_pilot.line;j++) { 
		m << "["; 
		for(i=0;i<X_pilot.col;i++) {
			m << X_pilot.mat[j][i] << " "; 
		} 
		m << "]" << endl; 
	} 
	m << "];" << endl << endl;

	m << "sigma_h = ["; 
	for(j=0;j<sigma_h.line;j++) { 
		m << "["; 
		for(i=0;i<sigma_h.col;i++) {
			m << sigma_h.mat[j][i] << " ";
		} m << "]" << endl ; 
	} 
	m << "];" << endl << endl;
	
	m << "sigma_h_cross = ["; 
	for(j=0;j<sigma_h_cross.line;j++) {
		m << "[";
		for(i=0;i<sigma_h_cross.col;i++) {
			m << sigma_h_cross.mat[j][i] << " ";
		} 
		m << "]" << endl; 
	} 
	m<<endl << "];" << endl << endl;

	m << "path_1_1 = ["; 
	for(j=0;j<path.mat[0][0].line;j++) {
		m << "["; 
		for(i=0;i<path.mat[0][0].col;i++) { 
			m << path.mat[0][0].mat[j][i] << " ";
		}
		m << "]" << endl; 
	}
	m<<endl << "];" << endl << endl;

	m << "path_2_1 = ["; 
	for(j=0;j<path.mat[1][0].line;j++) { 
		m << "["; 
		for(i=0;i<path.mat[1][0].col;i++) { 
			m << path.mat[1][0].mat[j][i] << " "; 
		}
		m << "]" << endl; 
	}
	m<<endl << "];" << endl << endl;

	m << "h_curr = ["; 
	for(j=0;j<h.line;j++) { 
		m << "["; 
		for(i=0;i<h.col;i++) { 
			m << h.mat[j][i] << " "; 
		}
		m << "]" << endl; 
	}
	m<<endl << "];" << endl << endl;

	m << "noise  = ["; 
	for(j=0;j<noise.vect[0].line;j++) { 
		m << "["; 
		for(i=0;i<noise.vect[0].col;i++) { 
			m << noise.vect[0].mat[j][i] << " "; 
		}
		m << "]" << endl; 
	}
	m<<endl << "];" << endl << endl;

	m << "W = ["; 
	for(j=0;j<W.line;j++) { 
		m << "["; 
		for(i=0;i<W.col;i++) { 
			m << W.mat[j][i] << " "; 
		}
		m << "]" << endl; 
	} 
	m<<endl << "];" << endl << endl;

	m.close();

} 
#endif

#define _CONVERGE_DEBUG_

#ifdef _CONVERGE_DEBUG_
#pragma message ("CONVERGE DEBUG IS ON")
#endif

int main(int argc, char** argv) {
	int nk=0;
	int EM_ITER_MAX = 3;
	int Lf = 2;
	double SNR_DB = 5; 
	int omega = 4; // on utilise PSK_GRAY_4
	int log2_omega = (int)log2(omega);
	int STC_M = 1, STC_N=2, STC_P=2, OFDM_K=8, Q=10;

	int LDPC_M = 20;
	int LDPC_N=LDPC_M * STC_N /(STC_N-1); // rendement 1/STC_N
	WLDPCParityCheck ldpc(LDPC_M, LDPC_N);
	ldpc.make_dense_mixed();

	int k_factor = LDPC_M;
	int basic_bit_in_size = Q*STC_P*OFDM_K*log2_omega;
	int bit_in_size = k_factor*basic_bit_in_size;
	//int k_factor2

	WVector bit_in(bit_in_size);
	TRandBitStatePtr bit_seed;
	bit_seed.iseed=5547557;
	bRandBit1(&bit_seed, bit_in);
	WVector ldpc_coded_tx = ldpc.encode(bit_in); // on a k_factor*Q*STC_N*P*K*log2_omega bits
	ZVector psk_tx = z4PSKGray_mod(ldpc_coded_tx); // on a k_factor*Q*NPK symboles
	ZVector x_em(psk_tx.taille);
#ifdef _CONVERGE_DEBUG_
	ofstream conv("converge.log");
	conv << "Starting simulation" << endl;
	conv.close();
#endif
#ifdef _MATLAB_DEBUG_
	ofstream matlab("matlab.log");
	matlab << "# Compiled : " << __TIMESTAMP__ " , FILE : " <<  __FILE__ << endl;
	matlab << "# omega = " << omega << ", STC_M = " << STC_M << ", STC_N = " << STC_N << ", STC_P = " << STC_P << "OFDM_K = " << OFDM_K << endl;
	matlab << "# binary informations " << endl << bit_in;
	matlab << "# LDPC Informations " << endl << ldpc;
	matlab << "# LDPC coded " << endl << ldpc_coded_tx ;
	matlab << "# PSK4 gray  modulated symbols " << endl << psk_tx;
	matlab.close();
#endif

#undef ant_tx_debug_
#ifdef ant_tx_debug_
	for(int xx=0;xx<psk_tx.taille;xx++) {
		psk_tx.vect[xx]=xx;
	}
#endif

	UVector<ZVector> ant_tx(STC_N);
	int PK = STC_P*OFDM_K;
	int QPK = Q*PK;

	int m,n, i, j;

    // generating pilots for further use .. 
	ZMatrix X_pilot(OFDM_K, STC_N*OFDM_K);
	UVector<ZVector> spilot(STC_N);
	for(i=0;i<STC_N;i++) {
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

	for(n=0;n<STC_N;n++) {
		for(int k=0;k<OFDM_K;k++) {
			X_pilot.mat[k][n*OFDM_K+k]=spilot.vect[n].vect[k];
		}
	}

	// end of pilot generation .. 
#ifdef _MATLAB_DEBUG_
	matlab.open("matlab.log", ios::app);
	matlab << "# symbol pilot sequences for each antennas .. " << endl << spilot;
	matlab << "# these pilots are placed in X formed matrix " << endl << X_pilot;
	matlab.close();
#endif

	// Stream splitting 
	int kQPK = k_factor*QPK;
	for(n=0;n<STC_N;n++) { // for each stream 
		ant_tx.vect[n] = ZVector(k_factor*QPK);
	}
	for(int t=0, tt=0;t<kQPK;t++) {
		for(int n=0;n<STC_N;n++, tt++) 
		ant_tx.vect[n].vect[t]= psk_tx.vect[tt];
	}

#ifdef _MATLAB_DEBUG_
	matlab.open("matlab.log", ios::app);
	matlab << "# Data spliting : each antennas will send these informations given the original psk information (pilots excluded)" << endl;
	matlab << ant_tx;
	matlab.close();
#endif

	// Stream splitting end

	// now creating the W matrix ... 
	ZMatrix W(2*OFDM_K, 2*Lf, DCplx());
	double UTILITIS_M_2_PI_OFDM_K = -UTILITIS_M_2_PI/OFDM_K;
	for(i=0;i<OFDM_K;i++) {
		for(j=0;j<Lf;j++) {
			W.mat[OFDM_K+i][Lf+j]=W.mat[i][j]=DCplx::exp2RI(1, UTILITIS_M_2_PI_OFDM_K*i*j);
		}
	}

#ifdef _MATLAB_DEBUG_
	matlab.open("matlab.log", ios::app);
	matlab << "# The W matrix is shown below " << endl << W << endl;
	matlab.close();
#endif

	// we now need to generate all fading values .. 
	//each antennas will transmit K_factor*LDPC_M*QPK symbols .. hence, we need the same number of fading coef.

	UMatrix<ZMatrix> taps(STC_N, STC_M);
	UMatrix<dRandUniStatePtr> tapseed(STC_N, STC_M);
	DMatrix tappower=DMatrix::eye(Lf);
	DMatrix Sig_hcc(STC_N*Lf, STC_N*Lf), Sigma_H(STC_N*Lf, STC_N*Lf);
	
	switch (Lf) {
		case 3 : 
			tappower.mat[2][2]=sqrt(0.05);
			Sig_hcc.mat[Lf+2][Lf+2]=Sig_hcc.mat[2][2]=tappower.mat[2][2]!=0 ? 1.0/(tappower.mat[2][2]*tappower.mat[2][2]) : 0.0;
			Sigma_H.mat[2][2]=Sigma_H.mat[Lf+2][Lf+2] = tappower.mat[2][2]*tappower.mat[2][2];
		case 2: 
			tappower.mat[1][1]=sqrt(0.1);
			Sig_hcc.mat[1][1]=Sig_hcc.mat[Lf+1][Lf+1]=tappower.mat[1][1]!=0 ? 1.0/(tappower.mat[1][1]*tappower.mat[1][1]) : 0.0;
			Sigma_H.mat[1][1]=Sigma_H.mat[Lf+1][Lf+1] = tappower.mat[1][1]*tappower.mat[1][1];
		case 1: 
			tappower.mat[0][0]=sqrt(0.8);
			Sig_hcc.mat[0][0]=Sig_hcc.mat[Lf][Lf]=tappower.mat[0][0]!=0 ? 1.0/(tappower.mat[0][0]*tappower.mat[0][0]) : 0.0;
			Sigma_H.mat[0][0]=Sigma_H.mat[Lf][Lf] = tappower.mat[0][0]*tappower.mat[0][0];
	}

	// end tappower

#ifdef _MATLAB_DEBUG_
	matlab.open("matlab.log", ios::app);
	matlab << "# for this run, the channel has "<<Lf << " taps and the power matrix is " << endl << tappower << endl;
	matlab << "# this bring us to the values of Sigma_h and Sigma_h_cross (Eq 21-22)" << endl << Sigma_H << Sig_hcc << endl;
	matlab.close();
#endif


	double baudrate=1e6;
	double max_doppler;
	double symbol_time = 1.0e6/baudrate;
	max_doppler=calc_doppler(baudrate, 100);
	
	int k_factor_QPK = k_factor*QPK; // optimizer 
	
	double snr_variance = convertsnrtosigma2(SNR_DB);
	double inv_snr_var = 1.0/snr_variance;

	zRandGausStatePtr zptr;
	zRandGausInit(&zptr, 0, snr_variance);

	int training_length=1;

	for(m=0;m<STC_M;m++) {
		for(n=0;n<STC_N;n++) {
			dRandUniInit( &(tapseed.mat[n][m]), 311254488+3*n-2*m, 0.0, 1.0);
			taps.mat[n][m]= fadingtaps(max_doppler, symbol_time, Lf, k_factor*(training_length+Q*STC_P), 512,&(tapseed.mat[n][m]))*tappower;
		}
	}

#ifdef _MATLAB_DEBUG_
	matlab.open("matlab.log", ios::app);
	matlab << "# The Channel parameter is shown below " << endl ;
	matlab << "# SNR_dB : "<< SNR_DB << ", SNR_Variance ; " << snr_variance << endl;
	matlab << "# taps  .. " << endl << taps << endl;
	matlab.close();
#endif

	int tapstime=0, tx_time=0;


	for(int kf=0; kf<k_factor;kf++) { // actually , we have k_factor piloted block
#ifdef NEW_PILOT_SEQUENCING
		WVector bpilot(training_length*OFDM_K*log2_omega);
		ZVector xpilot=z4PSKGray_mod(bpilot);
		ZMatrix wpilot = ColumnPackedVectorAsMatrix(xpilot, OFDM_K, training_length);

		ZMatrix hpilot(OFDM_K, training_length);
		for(int tl=0;tl<training_length;tl++, tapstime++) {
			for(int k=0;k<OFDM_K;k++) {
				hpilot.mat[





#else
		// we have pilots here ... 
		// estimating channel with aide of pilots ...
		// calculating ypilot

		UVector<ZMatrix> ypilot(STC_M), y(STC_M);
		UVector<ZMatrix> zpilot(STC_M), z(STC_M);
		UVector<ZMatrix> h(STC_M), h_est(STC_M);

		int current_symbol_time = kf*QPK;

		for(m=0;m<STC_M;m++) {
			h.vect[m] = ZMatrix(STC_N*Lf,1);
			for(n=0;n<STC_N;n++) {
				for(int l=0;l<Lf;l++) {
					h.vect[m].mat[n*Lf+l][0] = taps.mat[n][m].mat[tapstime][l];
				}
			}
		} // m for h
#ifdef _MATLAB_DEBUG_
		matlab.open("matlab.log", ios::app);
		matlab << "# real h for current time (pilot mode) [taptime = " << tapstime << "]" << endl << h; 
#endif
		for(m=0;m<STC_M;m++) {
			zpilot.vect[m]=ZMatrix(OFDM_K, 1);
			for(int z=0;z<OFDM_K;z++) zpilot.vect[m].mat[z][0]=zRandGaus(&zptr);
			ypilot.vect[m]=X_pilot*W*h.vect[m]+zpilot.vect[m];
		} // m for ypilot

#ifdef _MATLAB_DEBUG_
		matlab <<  "# Noise " << endl << zpilot << endl;
		matlab << "# received signal at rx antenna " << endl << ypilot;
		matlab.close();
#endif

		// with knowledge of pilots, we can estimate h directly
		ZMatrix hcleft = inv_snr_var*W.hermit()*X_pilot.hermit()*X_pilot*W+Sig_hcc;
		ZMatrix pivot(hcleft.line, hcleft.col);

		//cout << hcleft ; 
		gaussj(hcleft, pivot);
	//cout << hcleft ; 
#ifdef _MATLAB_DEBUG_
		matlab.open("matlab.log", ios::app);
		matlab << "# inv(W'X'XW + sigma_h_cross) " << endl << hcleft << endl;
		matlab.close();
#endif

		hcleft = hcleft*W.hermit()*X_pilot.hermit()*inv_snr_var;

#ifdef _MATLAB_DEBUG_
		matlab.open("matlab.log", ios::app);
		matlab << "# inv(W'X'XW + sigma_h_cross)W'X' " << endl << hcleft << endl;
		matlab.close();
#endif

		for(m=0;m<STC_M;m++) {
			h_est.vect[m]=hcleft*ypilot.vect[m];
			//cout << (h_est.vect[m]-h.vect[m]).mag();
		}
#ifdef _MATLAB_DEBUG_
		matlab.open("matlab.log", ios::app);
		matlab << "# estimated channel values " << endl << h_est << endl;
		matlab.close();

		writematlab(snr_variance, psk_tx, Sigma_H, Sig_hcc, taps, h.vect[0], zpilot, W, X_pilot); 
#endif
#ifdef _CONVERGE_DEBUG_
				conv.open("converge.log", ios::app);
				conv << "Estimated channel diff with aide of pilots : " << (h_est.vect[0]-h.vect[0]).mag() << endl;
				conv.close();
#endif
#endif

// JUSQU'ICI .. les calculs se concordent avec MATLAB ... 
 
		for(int q=0;q<Q;q++) { 
			for(int p=0;p<STC_P;p++) {
				//tx_time --> current_symbol_time = kf*QPK + q * PK + p*OFDM_K;
				ZMatrix X_tx(OFDM_K, STC_N*OFDM_K), X_kappa(OFDM_K, STC_N*OFDM_K), X_kappa_next(OFDM_K, STC_N*OFDM_K), Sigma_hchap;
				/*for(int n=0;n<STC_N;n++) {
					int nK= n*OFDM_K;
					for(int k=0;k<OFDM_K;k++) {
						X_tx.mat[k][nK+k] = ant_tx.vect[n].vect[current_symbol_time+k];
					}
				}*/
				int txtime_copy=tx_time;
				for(int kt=0;kt<OFDM_K;kt++) {
					for(int n=0;n<STC_N;n++) {
						X_tx.mat[kt][n*OFDM_K+kt]=ant_tx.vect[n].vect[tx_time];
					}
					tx_time++;
				}

#ifdef _MATLAB_DEBUG_
		matlab.open("matlab.log", ios::app);
		matlab << "# TX time is from " << txtime_copy << " to " << tx_time << endl;
		matlab << "# time control : " <<  kf*QPK + q * PK + p*OFDM_K-txtime_copy << ",  " <<  kf*QPK + q * PK + p*OFDM_K+OFDM_K-tx_time << endl;
		matlab << "# antenna 1 symbols : " << endl << ant_tx.vect[0].copy(txtime_copy, tx_time) << endl;
		matlab << "# antenna 2 symbols : " << endl << ant_tx.vect[1].copy(txtime_copy, tx_time) << endl;
		matlab << "# X_tx " << endl << X_tx << endl;
		
		matlab.close();
#endif // first pass ok with matlab
				for(m=0;m<STC_M;m++) {
					for(int n=0;n<STC_N;n++) {
						for(int l=0;l<Lf;l++){
							h.vect[m].mat[n*Lf+l][0]=taps.mat[n][m].mat[tapstime][l];
						}
					}
				}
				tapstime++;

				for(m=0;m<STC_M;m++) {
					z.vect[m]=ZMatrix(OFDM_K, 1);
					for(int cz=0;cz<OFDM_K;cz++) z.vect[m].mat[cz][0]=zRandGaus(&zptr);
					y.vect[m] = X_tx*W*h.vect[m]+z.vect[m];
				}

#ifdef _MATLAB_DEBUG_
		matlab.open("matlab.log", ios::app);
		matlab << "# taps time" << tapstime -1 << endl;
		matlab << "# real channel parameter " << h << endl;
		matlab << "# estimated channel parameter " << h_est << endl;
		matlab << "# noise  " <<  z << endl; 
		matlab << "# received symbol " << endl << y << endl;
		matlab.close();
#endif
				// il faut maintenant determine X_kappa (0).. algo ML
				for(int k=0;k<OFDM_K;k++) {
					ZMatrix xk(2,1);
					double qi=0.0, qi_min = 1e20; ZMatrix xkmin(2,1); 
					ZMatrix Wfp(2, 2*Lf);
					//calcul de Wfp
					double UTILITIS_M_2_PI_K = -UTILITIS_M_2_PI*k/OFDM_K;
					for(int l=0;l<Lf;l++) {
						Wfp.mat[0][l] = Wfp.mat[1][Lf+l] = DCplx::exp2RI(1,UTILITIS_M_2_PI_K*l);
					}
					//end Wfp
#ifdef _MATLAB_DEBUG_
		matlab.open("matlab.log", ios::app);
		matlab << "# Wf'  " << endl << Wfp << endl;
		matlab << "# current k time : " << k << endl;
		matlab << "# qi behavior : " << endl;
		matlab.close();
#endif
					for(int p1=0;p1<omega;p1++) {
						for(int p2=0;p2<omega;p2++) {
							xk.mat[0][0]= psk4gray.sym[p1].conj();
							xk.mat[1][0]= psk4gray.sym[p2].conj();
							qi = (y.vect[0].mat[k][0] - (xk.hermit()*Wfp*h_est.vect[0]).mat[0][0]).mod2();
#ifdef _MATLAB_DEBUG_
		matlab.open("matlab.log", ios::app);
		matlab << "# (P1, P2) = (" << p1 << "," << p2 << ")" << " (qi / qi_min) = (" << qi << ", " << qi_min << ")"  << endl;
#ifdef _ML_DEBUG_
		matlab << "# current symbol " << endl << xk << endl;
		matlab << "# y[k] = " << y.vect[0].mat[k][0] << endl;
		matlab << "# h_est = " << h_est.vect[0] << endl;
		matlab << "# xk'*wfp*h_est = " << (xk.hermit()*Wfp*h_est.vect[0]) << endl;
		matlab.close();
#endif
#endif
							if(qi<qi_min) {
								qi_min = qi;
								xkmin = xk;
							}
						}
					}
					// on a x_min pour k en conjugue .. 
					X_kappa.mat[k][k] = xkmin.mat[0][0].conj();
					X_kappa.mat[k][OFDM_K+k] = xkmin.mat[1][0].conj();
#ifdef _MATLAB_DEBUG_
		matlab.open("matlab.log", ios::app);
		matlab << "# min X couple symbol found to be : " << xkmin << endl;
		matlab.close();
#endif

					//cout << "k=" << k<<endl;
				}
				// fin de X_kappa(0);
#ifdef _CONVERGE_DEBUG_
				conv.open("converge.log", ios::app);
				conv << "X(0)  diff : " << (X_kappa - X_tx) << endl;
				conv << "channel dif a X0 : " << (h_est.vect[0]-h.vect[0]) << endl;
				conv.close();
#endif
#ifdef _MATLAB_DEBUG_
		matlab.open("matlab.log", ios::app);
		matlab << "# X_kappa(0)" << X_kappa << endl;
		matlab.close();
#endif

				// l'algo em peut commencer ici .. 

				hcleft = inv_snr_var*W.hermit()*X_kappa.hermit()*X_kappa*W+Sig_hcc;
				gaussj(hcleft, pivot);
				hcleft = hcleft*W.hermit()*X_kappa.hermit();

				for(int iter=0;iter<EM_ITER_MAX;iter++) {

					// reestimation de h_est
					//h_est.vect[0] = hcleft*inv_snr_var*y.vect[0];
					//cout << "Diff : " << (h_est.vect[0]-h.vect[0]).mag();
					Sigma_hchap = Sigma_H-hcleft*X_kappa*W*Sigma_H*inv_snr_var;
					// end of reestimation 

					for(int k=0;k<OFDM_K;k++) {
						ZMatrix xk(2,1); 
						double qi=0.0, qi_min = 1e20, qi_right=0.0; 
						ZMatrix xkmin(2,1); 
						ZMatrix Wfp(2, 2*Lf);
						//calcul de Wfp
						double UTILITIS_M_2_PI_K = -UTILITIS_M_2_PI*k/OFDM_K;
						for(int l=0;l<Lf;l++) {
							Wfp.mat[0][l] = Wfp.mat[1][Lf+l] = DCplx::exp2RI(1,UTILITIS_M_2_PI_K*l);
						}
						
						ZMatrix subSigma_hchap(2,2), ztemp;
						ztemp = W*Sigma_hchap*W.hermit();
						subSigma_hchap.mat[0][0]= ztemp.mat[k][k];
						subSigma_hchap.mat[0][1]= ztemp.mat[k][OFDM_K+k];
						subSigma_hchap.mat[1][0]= ztemp.mat[OFDM_K+k][k];
						subSigma_hchap.mat[1][1]= ztemp.mat[OFDM_K+k][OFDM_K+k];

						//end Wfp
						for(int p1=0;p1<omega;p1++) {
							for(int p2=0;p2<omega;p2++) {
								xk.mat[0][0]= psk4gray.sym[p1].conj();
								xk.mat[1][0]= psk4gray.sym[p2].conj();
								// a verifier le sens de matrice y
								qi = (y.vect[0].mat[k][0] - (xk.hermit()*Wfp*h_est.vect[0]).mat[0][0]).mod2();
								qi_right = (xk.hermit()*subSigma_hchap*xk).mat[0][0].re;
								qi+=qi_right;
								
								if(qi<qi_min) {
									qi_min = qi;
									xkmin = xk;
								}
							}
						}
						// on a x_min pour k en conjugue .. 
						X_kappa_next.mat[k][k] = xkmin.mat[0][0].conj();
						X_kappa_next.mat[k][OFDM_K+k] = xkmin.mat[1][0].conj();
					}
					X_kappa = X_kappa_next;
					hcleft = inv_snr_var*W.hermit()*X_kappa.hermit()*X_kappa*W+Sig_hcc;
					gaussj(hcleft, pivot);
					hcleft = hcleft*W.hermit()*X_kappa.hermit();
					h_est.vect[0] = hcleft*inv_snr_var*y.vect[0];
#ifdef _CONVERGE_DEBUG_
				conv.open("converge.log", ios::app);
				conv << "X Diff a l'iter " << iter << (X_kappa-X_tx) ;
				conv << "channel diff a lter : " << iter  << (h_est.vect[0]-h.vect[0]);
				conv.close();
#endif
				}

				// em termine .. 
				current_symbol_time = kf*QPK + q * PK + p*OFDM_K;
				for(int kk=0;kk<OFDM_K;kk++) {
					x_em.vect[nk++]= X_kappa.mat[kk][kk];
					x_em.vect[nk++] = X_kappa.mat[kk][OFDM_K+kk];
				}// end em copy
	
				cout << "current_symbol_time ="<< nk << " P, Q = " << p << "," << q << endl;
			} // p
		} // q
	} // kf
	int counter=0;
	for(i =0; i<x_em.taille;i++)  if(x_em.vect[i]!=psk_tx.vect[i])counter++;
	cout << "ERROR : " << counter << " out of " << x_em.taille << endl;
	cout << "EM ITER = " << EM_ITER_MAX << " , SNR_DB (snr_variance) = " << SNR_DB << "("<<  snr_variance << ")" << endl;
	cout << "Lf = " << Lf << endl;

	WVector txcnt(4,0), rxcnt(4,0);
	for(i = 0; i<x_em.taille;i++) {
		for(int o=0;o<4;o++) {
			if(x_em.vect[i]==psk4gray.sym[o]) rxcnt.vect[o]++;
			if(psk_tx.vect[i]==psk4gray.sym[o]) txcnt.vect[o]++;
		}
	}
	
	for(i=0;i<4;i++) cout << i << " : " << txcnt.vect[i] << " --> " << rxcnt.vect[i] << endl;


	return 1;
}
#else

int main(void) {

	DMatrix A(30,20), U, V, W;
	register int i, j;
	dRandGausStatePtr ptr;
	dRandGausInit(&ptr, 1, 37);
	for(i=0;i<A.line;i++) {
		for (j=0;j<A.col;j++) {
			A.mat[i][j]=dRandGaus(&ptr);
		}
	}

	if(svd(A, U, W, V)!=0) 
		cout << "SVD failed ... " << endl;
	else {
		ofstream svdlog;
		svdlog.precision(17);
		svdlog.open("svd.log");
		svdlog<< A << U << W << V << endl;
		svdlog.close();
		cout << "Done " << endl;
	}

	return 0;


}
#endif