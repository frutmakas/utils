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

#include "tools/interleaver.h"
#ifdef WIN_DOS
#include <conio.h>

#endif


#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 

#ifndef SPOOKY
const int OFDM_SIZE = 16;

#pragma message("Message :> Compiling code for complete simulation ")
int menu(){
	cout << endl << endl;
	cout << "1) Generate and Save LDPC Generator Matrix" << endl ;
	cout << "2) Generate, Save LDPC Generator Matrix and Launch Simulation"<< endl ;
	cout << "3) Load LDPC Generator Matrix and Launch Simulation"<< endl ;
	cout << "4) Quit"<< endl ;
	cout << endl << "Your Choice : ";
	cout.flush();
	int choice;
	do {
#ifdef WIN_DOS
		choice = getch(); //fgetc(stdin);
#else
		choice = getchar();
#endif
	} while (choice >=0 && choice <=4);
	cout << endl;
	return choice;
}

void Generate(){
	cout << "Status :> Entering LDPC and Chirps Generator" << endl << "Status :> Need LDPC Matrix Dimension"<< endl ;
	int M=0,N=0;
	bool dimension_invalid=true;
	do {
		cout << "Query :> Please enter LDPC Matrix Length and Column (separate them with a space)"<< endl ;
		cout << "Answer :> ";
		cout.flush();
		cin >> M >> N;
		dimension_invalid = M>N || N<=0 ||M<=0;
		if(dimension_invalid) {
			cout << "Error :> Matrix size is invalid !" << endl;
		}
	}while (dimension_invalid);

	char filename[80];
	do {
		cout << "Query :> Please enter a file name to save the data (existing file will be erased)" << endl;
		cout << "Answer :> ";
		cout.flush();
		// fgets(filename, 79,stdin);
		try {
			cin >> filename;
		} catch (...){};

	} while (!strcmp(filename, ""));


#ifndef NO_LDPC

#pragma message("Message :> Using LDPC in this binary")
	cout << "Status :> Generating LDPC Matrix with Size "<<M<<"x"<<N<<" and 3 bit checks with Dense Evencol method. Please wait ..." << endl;

	make_method method=Evencol;
	WLDPCParityCheck ldpc(M,N,3,method);
	method= Dense;
	ldpc.make_dense_mixed(method,1);

#else

#pragma message("Message :> LDPC is not used in this binary")

#endif



	int k_factor;

	cout << "Query :> Enter k factor : "; cout.flush();

	cin >> k_factor;

	cout << "Status :> Generating "<< (N-M)*3*2*OFDM_SIZE*k_factor << " random bits to be used for simulation" << endl;
	WVector bin_in((N-M)*3*2*OFDM_SIZE*k_factor);
	TRandBitStatePtr seed;
	bRandBit1(&seed,bin_in);

#ifndef  NO_LDPC

	cout << "Status :> Encoding generation bit using the previous LDPC matrix "<< endl;
	WVector bin_enc=ldpc.encode(bin_in);

	cout << "Status :> Modulating encoded data using PSK8 "<< endl;

	ZVector psk_tx=z8PSK_mod(bin_enc);

#else 

	cout << "Status :> Modulating data using PSK8 "<< endl;

	ZVector psk_tx=z8PSK_mod(bin_in);

#endif

	cout  << "Status :> Saving all necessary data in " << filename << endl;
	ofstream ofs(filename);
#ifndef  NO_LDPC

	cout << "Status :> Saving LDPC in file" <<  endl;
	ofs << ldpc;

#endif
	cout << "Status :> Saving input bits in file"<< endl;
	ofs << bin_in;
	cout << "Status :> Saving modulated signal in file"<<endl;
	ofs << psk_tx;

	cout << "Status :> Closing file ..." << endl;
	ofs.close();
	cout << "Status :> Preparing to leave generator"<< endl;
	ofs.open("lastfile.dat");
	ofs << filename ; //<< endl;
	ofs.close();
	cout << "Status :> Exiting LDPC Generator"<< endl;
	                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
}

void LoadLDPC(){
	cout << "Status :> Looking for saved raw data file for simulation" << endl;
	char filename[80];
	bool failure=true;
	do {
		do {
			cout << "Query :> Enter raw data filename " << endl;
			cout << "Answer :> ";
			cout.flush();
			//fgets(filename,80,stdin);
			cin>>filename;
		}while (!strcmp(filename,""));
		FILE *f =fopen(filename, "rt");
		failure = f==NULL;
		if (!failure) fclose(f);
	} while (failure);
	
	ofstream ofs("lastfile.dat");
	ofs << filename; //<< endl;
	ofs.close();
	cout << "Status :> Exiting LoadLDPC" << endl;
}

void Simulate(){
	char infile[80];
	cout << "Status :> Entering Simulation " << endl;
	//FILE *fi = fopen("lastfile.dat", "rt");
	ifstream fi("lastfile.dat");
	//fgets(infile, 80, fi);
	//fclose(fi);
	fi >> infile;
	fi.close();
	cout << "Status :> Loading raw data file " << infile << endl;
	ifstream ifs(infile);

#ifndef  NO_LDPC

	WLDPCParityCheck ldpc;
	cout << "Status :> Loading ldpc" << endl;

	ifs >> ldpc;

#endif

	WVector bin_in;
	cout << "Status :> Loading binary information" << endl;
	ifs >> bin_in;
	ZVector psk_tx;
	cout << "Status :> Loading PSK8 modulated signal " << endl;
	ifs >> psk_tx;
	ifs.close();

	long max_doppler;
	double symbol_time;
	long filter_length=1;
	long nb_rays;

	double  baudrate;
	cout << "Query :> Please enter PSK8 signal baud rate" << endl;
	cout << "Answer :> "; cout.flush();
	cin  >> baudrate;

	symbol_time = 1.0e6/baudrate;
	nb_rays=76;
	//max_doppler=(long)floor(1e-2*baudrate);

	max_doppler=(long)floor(3.3356409519815204957557671447492e-4*baudrate); // v = 100km/h, c=299792458m/s
	int sweep_yes=0;
	do {
		cout << "Status :> Determining noise simulation method"<< endl;
		cout << " 1) One shot simulation " << endl;
		cout << " 2) Sweep noise Eb/N0 (dB) level "<<endl;
		cout << " Answer :> "; cout.flush();
		cin >> sweep_yes;
	} while (sweep_yes<1 || sweep_yes>2);
	double ebn0_min, ebn0_max, ebn0_step, ebn0;
	int ebn0_count=1;
	switch (sweep_yes) {
		case 1:  	
			    cout << "Query :> Please enter noise Eb/N0 (dB) level) " << endl;
				cout << "Answer :> "; cout.flush();
				cin >> ebn0;
				break;
		case 2: 
				cout << "Query :> Please enter Eb/N0 (dB) min value" << endl;
				cout << "Answer :> "; cout.flush();
				cin >> ebn0_min;
				cout << "Query :> Please enter Eb/N0 (dB) max value" << endl;
				cout << "Answer :> "; cout.flush();
				cin >> ebn0_max;
				cout << "Query :> Please enter Eb/N0 (dB) step value" << endl;
				cout << "Answer :> "; cout.flush();
				cin >> ebn0_step;
				if(ebn0_min>ebn0_max) swap(ebn0_min, ebn0_max);
				ebn0_count = (int)floor((ebn0_max-ebn0_min)/ebn0_step);
				break;
	}

	double Eb=3.0;

	cout << "Status :> Preparing to simulate"<<endl;
	char report_filename[120];
	do {
		cout << "Query :> Enter report filename "<<endl;
		cout << "Answer :> "; cout.flush();
		cin >> report_filename;
	} while(!strcmp(report_filename, ""));


	double ofdm_noise_normalization = 1.0/3.0/OFDM_SIZE; 

#ifndef NO_LDPC
	ofdm_noise_normalization *= (double)(ldpc.Col()-ldpc.Line())/(double)ldpc.Col();
#endif

#ifndef USE_INTERLEAVER
	cout << "Status :> Generating flat fading coefficient with max_doppler = " << max_doppler << " and symbol time = " << symbol_time << "µs" << endl;

	dRandUniStatePtr tapseed;

	dRandUniInit(&tapseed, 314975412, 0.0, 1.0);
 	ZMatrix zmtaps = fadingtaps(max_doppler, symbol_time, filter_length, psk_tx.taille,nb_rays, &tapseed);
#else
	int interleaver_size;
	cout << "Status :> Using interleaver" << endl;
	do {
		cout << "Query :> Enter Interleaver size " << endl;
		cout << "Answer :> " ; cout.flush();
		cin >> interleaver_size;
	} while (interleaver_size < 2);
	cout << "Status :> Interleaver size is set to " << interleaver_size << endl;
	cout << "Status :> Generating flat fading coefficient with max_doppler = " << max_doppler << " and symbol time = " << symbol_time << "µs" << endl;
 	ZMatrix zmtaps = fadingtaps(max_doppler, symbol_time, filter_length, psk_tx.taille+interleaver_size-psk_tx.taille%interleaver_size,nb_rays);
	ofdm_noise_normalization *= (double)psk_tx.taille/(double)(psk_tx.taille+interleaver_size-psk_tx.taille%interleaver_size);
#endif

	ZVector ofdm_tx(psk_tx.taille);
//	ZVector fft_tmp(OFDM_SIZE);
	cout << "Status :> Preparing ofdm signals "<< endl;
	for(int f=0;f<psk_tx.taille;f+=OFDM_SIZE) {
			ofdm_tx.insert(ifft(psk_tx.copy(f,f+OFDM_SIZE)),f);
	}
	int psk_tx_taille = psk_tx.taille;
	psk_tx.~ZVector();

	for(int x=0; x<=ebn0_count;x++) {
		double noise_variance;
		switch (sweep_yes){
			case 1: 
				noise_variance = convertsnrtosigma2(ofdm_noise_normalization, ebn0);
				break;
			case 2: 
				ebn0 = ebn0_min + ebn0_step*x;
				noise_variance = convertsnrtosigma2(ofdm_noise_normalization, ebn0);
				break;
		}

		cout << "Status :> Starting simulation with " << psk_tx_taille << " symbols at Eb/N0 = "<< ebn0 << endl;

		zRandGausStatePtr noise_ptr;
		zRandGausInit(&noise_ptr,0,noise_variance);

#ifndef USE_INTERLEAVER
#pragma message("Message :> Not Using Interleaver")
		ZVector noise(ofdm_tx.taille);
		cout << "Status :> Generating additive white gaussian noise with noise variance = "<< noise_variance << endl;
		zbRandGaus(&noise_ptr,noise);
		ZVector ofdm_rx(ofdm_tx.taille);
		cout << "Status :> Applying noise and flat fading to transmitted signals " << endl;
		for(register int z=0;z<ofdm_tx.taille;z++) {

#ifndef NO_FLAT_FADING	
#pragma message("Message :> Compiled with Flat Fading influences")
			ofdm_rx.vect[z]=ofdm_tx.vect[z]*zmtaps.mat[z][0].mag()+noise.vect[z]; //
#else

#pragma message("Message :> Compiled with NO Flat Fading influences")
			ofdm_rx.vect[z]=ofdm_tx.vect[z]+noise.vect[z];

#endif
		}

#else // USE_INTERLEAVER is defined
#pragma message("Message :> Using Interleaver")
		CInterleaver interleaver(interleaver_size);
		cout << "Status :> Interleaving OFDM signal" << endl;
		ZVector interleaved_ofdm_tx = interleaver.Apply(ofdm_tx);
		ZVector interleaved_ofdm_rx(interleaved_ofdm_tx.taille);
		ZVector noise(interleaved_ofdm_tx.taille);
		cout << "Status :> Generating additive white gaussian noise with noise variance = "<< noise_variance << endl;
		zbRandGaus(&noise_ptr,noise); 
		cout << "Status :> Applying noise and flat fading to transmitted signals " << endl;
		for(register z=0;z<interleaved_ofdm_tx.taille;z++) {

#ifndef NO_FLAT_FADING	

#pragma message("Message :> Compiled with Flat Fading influences")
			interleaved_ofdm_rx.vect[z]=interleaved_ofdm_tx.vect[z]*zmtaps.mat[z][0].mag()+noise.vect[z];
#else //NO_FLAT_FADING	
#pragma message("Message :> Compiled with NO Flat Fading influences")
			interleaved_ofdm_rx.vect[z]=interleaved_ofdm_tx.vect[z]+noise.vect[z];
#endif //NO_FLAT_FADING	
		}
		cout << "Status :> Rearranging interleaved OFDM signal" << endl;
		ZVector ofdm_rx = interleaver.Extract(interleaved_ofdm_rx).copy(0, ofdm_tx.taille);
#endif // USE_INTERLEAVER
		cout << "Status :> Demodulating OFDM signal " << endl;
		ZVector psk_rx(ofdm_rx.taille);
		for(int i=0;i<psk_tx_taille;i+=OFDM_SIZE) {
				psk_rx.insert(fft(ofdm_rx.copy(i,i+OFDM_SIZE)),i);
		}

#ifdef OFDM_DEBUG
		ofstream off("ofdm_debug.txt");
		for(int o=0;o<OFDM_SIZE;o++) 
			off << ofdm_tx.vect[o] << ";" << ofdm_rx.vect[o] << ";" << psk_tx.vect[o] << ";" << psk_rx.vect[o] << ";" << zmtaps.mat[o][0] <<";" << zmtaps.mat[o][0].mag() << endl;
		off.close();
#endif //OFDM_DEBUG

#ifndef NO_LDPC
		cout << "Status :> Demodulating PSK8 signal " << endl;
		DVector ldpc_rx=z8PSK_demod(psk_rx).RZtoNRZ_nip();
		cout << "Status :> Decoding LDPC sequence (size = "<< bin_in.taille<<")" << endl;
		DVector llr;
		WVector bin_dec=ldpc.decode(ldpc_rx, llr);

#else 

		WVector bin_dec=z8PSK_demod(psk_rx);

#endif
	
		int cnt=0;
		cout << "Status :> Comparing signal and determining BER" << endl;
		for(int c=0;c<bin_dec.taille;c++){
			if (bin_in.vect[c]!=bin_dec.vect[c]) cnt++;
		}

		cout << endl << cnt << " errors out of " << bin_dec.taille << " ("<<((double)cnt)/bin_dec.taille << ")" << endl;

		// report

		time_t datetime;
		time(&datetime);

		cout << " Status :> Saving simulation data " << endl; 

		ofstream os;
		os.open(report_filename, ios::app);
		os << "*********************************************" << endl ; 
		os << "              Simulation result              " << endl;
		os << "*********************************************" << endl ; 
		os << " Program build date : " << __TIMESTAMP__ << endl;
		os << " Source file : " << infile << endl;
		os << " PSK influence is taken into account in SNR calculation" << endl;

		os << " Date & Time : " << ctime(&datetime) << endl;

#ifndef NO_LDPC

		os << " LDPC Matrix size : " << ldpc.Line()<<"x"<<ldpc.Col()<<endl;
		os << " Max iteration : " <<  MAX_ITERATION << endl;

#endif
		os << " Noise Eb/N0 (dB) : " << ebn0 << endl;
		os << " Noise variance : " << noise_variance << endl;
		os << " Baudrate : " << baudrate << endl;
		os << " Symbol time : " << symbol_time<< "µs" << endl;
		os << " Max Doppler : " << max_doppler<< endl;
		os << " Nb of bit simulated : " << bin_in.taille<<endl;
		os << " OFDM Packet size : " << OFDM_SIZE << endl;

#ifdef USE_INTERLEAVER
		os << " Interleaver size : " << interleaver_size<< endl;
#endif
		os << " Errors : " << cnt << endl;
		os << " BER    : " << ((double)cnt)/bin_dec.taille << endl;
#ifdef NO_FLAT_FADING	
		os << " #!! Simulation sans attenuation de canal !!#" << endl;
#endif

#ifdef USE_INTERLEAVER
		os << " #!!     Simulation avec l'entrelaceur    !!#" << endl;
#endif


		os << "*********************************************" << endl ; 
		os.close();
		if (cnt==0) break;
	} // noise sweep loop 

}

int main(int argc, char ** argv){
	cout << endl << "Program source filename : " << __FILE__ << endl << "Build : " << __TIMESTAMP__ << endl;
#ifdef NO_FLAT_FADING
	cout << endl << "No Flat Fading" << endl;
#endif
#ifdef NO_LDPC 
	cout << endl << "No LDPC " << endl;
#endif

#ifdef SHARED_NORMALIZATION
    cout << endl << "Shared normalization" << endl;
#endif
	bool freakout=false;
	do {
		switch(menu()){
			case '1':  // generate, save 
				Generate();
				break;
			case '2': // generate, save and simulate
				Generate();
				Simulate();
				break;
			case '3': // load and calculate
				LoadLDPC();
				Simulate();
				break;
			case '4': // quit
				freakout=true;
				break;
		}		
	} while (!freakout);

	return 1;
}

#else 

#endif
