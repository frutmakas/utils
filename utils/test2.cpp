#include "tools/utilitis.h"
#include <fstream.h>
#include "coding/block/ldpc.h"
#include "rand/randgen.h"
#include "modulation/psk.h"
#include "globaldef.h"
#include "fading/fadingtdma.h"
//int MAX_ITERATION = 30;
void main(){
	int M=200, N=400;
	WVector mi(1);
	//mi.vect[0] = 200;
	//mi.vect[1] = 400;
	mi.vect[0] = 1024;
	//mi.vect[3] = 2048;
	
	int x=0;
#ifndef NO_FLAT_FADING
	cout << "fading is simulated" << endl;
	double baudrate=1e6;
	double symbol_time = 1.0e6/baudrate;
	int nb_rays=76;
	//max_doppler=(long)floor(1e-2*baudrate);
	int filter_length=1;

	long max_doppler=(long)floor(3.3356409519815204957557671447492e-4*baudrate); // v = 100km/h, c=299792458m/s
#endif
	for(x=0, M = mi.vect[x], N=M*2;x<mi.taille;N=2*(M=mi.vect[++x])) {
		int sz =(int)floor(3000000.0/((N-M)*3));
		sz *=(N-M)*3;
		cout << "calculating " << sz << " bits for ldpc size = " << M << "x" << N << endl;
		WVector bin_in(sz);
		TRandBitStatePtr seed;
		bRandBit1(&seed, bin_in);
		WLDPCParityCheck ldpc;
		ldpc = WLDPCParityCheck(M,N);
		ldpc.make_dense_mixed();
		WVector bin_enc = ldpc.encode(bin_in);
		cout << "encoding PSK"	<< endl;
		ZVector psk_tx = z8PSK_mod(bin_enc);
#ifndef NO_FLAT_FADING
#pragma message("Compiling with fading influences ")
		dRandUniStatePtr tapseed;

		dRandUniInit(&tapseed, 314975412, 0.0, 1.0);
 		ZMatrix zmtaps = fadingtaps(max_doppler, symbol_time, filter_length, psk_tx.taille,nb_rays, &tapseed);
		DVector dtaps(psk_tx.taille);
		for(register int z=0;z<psk_tx.taille;z++){
			dtaps.vect[z]=zmtaps.mat[z][0].mag();
		}
#endif
		for(int snr=10;snr<20;snr++) {
			double noise_variance = convertsnrtosigma2((double)(N-M)/N, snr);
			zRandGausStatePtr noise_ptr;
			zRandGausInit(&noise_ptr, 0.0,noise_variance); 
			ZVector noise(psk_tx.taille);
			zbRandGaus(&noise_ptr, noise);
#ifndef NO_FLAT_FADING
			cout << " applying channel and noise effect to signal " << endl;
			ZVector psk_rx = dtaps*psk_tx+noise;
#else
#pragma message("Compiling with NO fading influences ")

			cout << " applying noise to signal " << endl;
			ZVector psk_rx = psk_tx+noise;
#endif
			DVector bin_dec = z8PSK_demod(psk_rx).RZtoNRZ_nip();
			cout << " decoding ... " << endl;
			DVector llr;
			WVector bin_out = ldpc.decode(bin_dec, llr);
			cout << "decoding done " << endl;
			int cnt=0;
			for(int i=0;i<bin_out.taille;i++){
				if(bin_in.vect[i]!=bin_out.vect[i]) cnt++;
			}

			ofstream os("new2_ldpc_iteration_result.txt", ios::app);
#ifndef NO_FLAT_FADING
			os << "Fading" << ";";
			cout << "Fading" << ";";
#else
			os << "no Fading"<<";";
			cout << "no Fading"<<";";
#endif
			os << MAX_ITERATION << ";" << M << ";" << N <<  ";" << snr << ";" << noise_variance << ";" << cnt << ";" << bin_out.taille << ";" << (double)cnt / (double) bin_out.taille << endl;
			os.close();
			cout << MAX_ITERATION << ";" << M << ";" << N <<  ";" << snr << ";" << noise_variance << ";" << cnt << ";" << bin_out.taille << ";" << (double)cnt / (double) bin_out.taille << endl;
			if (cnt==0) break;
		}
	}
}
