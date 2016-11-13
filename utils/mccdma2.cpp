#include <tools/all.h>
#include <rand/randgen.h>
#include <modulation/psk.h>
#include <cdma/walsh.h>
#include <cdma/spread.h>
#include <fading/fadingtdma.h>
#include <time.h>
#include <coding/block/ldpc.h>

//#define ldpc_test
//#define use_ldpc

#ifndef ldpc_test
//#define singleuser
#ifndef singleuser
int main() {
	int nb_user_max = 16;
	int snr_min,snr_max, snr_step;
	int cin_nb_block, cin_nb_frame;
//#define automat
#ifndef automat
	cout << "Enter nb of user :"; cout.flush();
	cin >> nb_user_max;

	cout << "Enter nb of frame :"; cout.flush();
	cin >> cin_nb_frame;
	cout << "Enter frame size:"; cout.flush();
	cin >> cin_nb_block;
	cout << "Enter snrmin:"; cout.flush();
	cin >> snr_min;
	cout << "Enter snrmax:"; cout.flush();
	cin >> snr_max;
	cout << "Enter snrstep:"; cout.flush();
	cin >> snr_step;
#else
	nb_user_max=8;
	cin_nb_frame=10;
	cin_nb_block=3;
	snr_min=0;
	snr_max=40;
	snr_step=3;
#endif
	UVector<DVector> code(nb_user_max);
	int i, j, k;//, l, t;
	WVector wtmp(nb_user_max);
	UVector<WVector> bit_in(nb_user_max);
	UVector<ZVector> sym_tx(nb_user_max);
	TConstellation constel = psk8gray;
	double powercoef=sqrt(0.5);
	for(int c=0;c<constel.size;c++) constel.sym[c]*=powercoef;	
	dRandUniStatePtr cseed;
	
	int initseed = time(NULL);
	dRandUniInit(&cseed, initseed, 0,1);

	ZMatrix Qh = ZMatrix::Q_DFT(nb_user_max).hermit();

	// chaque ligne represente la sequence d'etalement pour une sous porteuse
	UVector<DMatrix> ucode(nb_user_max);
	UVector<ZVector> sym_rx(nb_user_max);
	int ldpc_m = 20 , ldpc_n=40;
#ifdef use_ldpc
	UVector<WLDPCParityCheck> ldpc(nb_user_max);
#endif

	TRandBitStatePtr seed;
	seed.iseed=time(NULL);

	char filename[500], comment[5000];
	cout << "Comment : " << endl;
	cin.getline(comment, 5000);
	sprintf(filename, "mccdma-res_%d-%ld.txt", nb_user_max, initseed);
	ofstream os(filename, ios::app);
	int nb_block = 2*cin_nb_block*ldpc_m;
	os << "Time : " << initseed << endl << "nb user max = " << nb_user_max << endl << "nb block = " << nb_block << endl<<" Frame : " << cin_nb_frame <<endl;
	os << "Snr min : " << snr_min << endl << "snr_max : " << snr_max << endl << "snr_step : "<< snr_step<<endl;
	os << "Comment : " << comment << endl;
	cout << "ME TOUCHE PAS !! JE TRAVAILLEEEEEUUUHHHHHH !!!" << endl;

#ifdef use_ldpc
	os << "using ldpc .... " << endl;
#else
	os << "not using ldpc .... " << endl;
#endif

	WMatrix overall_teb(nb_user_max, snr_max), overall_sym(nb_user_max, snr_max);
	ofstream ldpc_out_file;
	//ldpc_out_file.open("ldpcin_file.dat");

	for(i=0;i<nb_user_max; i++) {
		code.vect[i] = dWalshCode(nb_user_max, i);
		wtmp.vect[i]=i;
#ifdef use_ldpc
		ldpc.vect[i]=WLDPCParityCheck(ldpc_m,ldpc_n);
		ldpc.vect[i].make_dense_mixed();
#endif
	}
	
	WMatrix coderep = ColumnCirculant(wtmp);
	
	UVector<ZMatrix> Pu(nb_user_max);

	// en ligne = frequence
	// en colonne = utilisateur
	// ex : coderep.mat[3][4] = code d'etalement a utiliser pour utilisateur 4 sur la sous porteuse 3

	for(i=0;i<nb_user_max;i++) {
		ucode.vect[i]=DMatrix(nb_user_max, nb_user_max);
		Pu.vect[i]=ZMatrix(nb_user_max, nb_user_max);
		for(j=0;j<nb_user_max;j++) { // ligne ==> la porteuse
			for(k=0;k<nb_user_max;k++) { // colonne == chip de la porteuse
				ucode.vect[i].mat[j][k] = code.vect[coderep.mat[j][i]].vect[k];
				Pu.vect[i].mat[j][coderep.mat[j][i]] = 1.0;
			}
		}
	}
	
	// ouff ... on a la matrice des code de chaque user
	zRandGausStatePtr nptr;
	typedef UMatrix<ZMatrix> UMZMatrix;
	typedef UVector<UMZMatrix> TSTBCChannel;
	
	int nb_rx = 1, r=0;
	
	ZMatrix c0transpose = ucode.vect[0].transpose_nip();
	UVector<ZMatrix> usend1(nb_user_max), usend2(nb_user_max);
	UVector<ZMatrix> sum1(nb_rx), sum2(nb_rx);
	ofstream bitin_file,/* ldpc_out_file,*/ sym_file_out, canal_file_out;
	
	for(int frame=0;frame<cin_nb_frame;frame++) {
		

		cout << "generating and channel coding user symbols " << endl;

		bitin_file.open("bitin_file.dat");
	       	//ldpc_out_file.open("ldpcin_file.dat");
		for(i=0;i<nb_user_max; i++) {
			WVector bit_in_tmp(nb_block*nb_user_max*constel.bit_per_sym);//(int)floor(log2(constel.size))); // deux blocs == un mot astbcm
			bRandBit1(&seed, bit_in.vect[i]);
#ifdef use_ldpc
#pragma message("using ldpc in "__FILE__)
//			cout << "ldpc encoding for user " << i << " with " << bit_in_tmp.taille << " data size" << endl;
			sym_tx.vect[i] = (*constel.modulate_function)(ldpc.vect[i].encode(bit_in_tmp));
			//ldpc_out_file << ldpc.vect[i];
#else
#pragma message("not using ldpc in "__FILE__)
			sym_tx.vect[i] = (*constel.modulate_function)(bit_in_tmp);
#endif
			sym_rx.vect[i] = ZVector(sym_tx.vect[i].taille);
			bitin_file << bit_in_tmp;
		}

		bitin_file.close();
		//ldpc_out_file.close();
		//ldpc.~UVector();
		
		// etalement de tt le monde

		ZMatrix noise;

		for(r=0;r<nb_rx;r++) {
			sum1.vect[r]=ZMatrix(nb_user_max, nb_user_max);
			sum2.vect[r]=ZMatrix(nb_user_max, nb_user_max);
		}


		TSTBCChannel timechannel(nb_user_max), canal(nb_user_max);


		int filter_length = 3;
							
		//ZVector tmp(nb_user_max);
		UVector< UMatrix<DVector> > canal_power(nb_user_max);
		for(i=0;i<nb_user_max;i++) {
			timechannel.vect[i]=UMZMatrix(2, nb_rx);
			canal.vect[i]=UMZMatrix(2, nb_rx);
			canal_power.vect[i] =UMatrix<DVector>(2,nb_rx);
			cout << "Generating user " << i << " channels ... " << endl;
			for(r=0;r<nb_rx;r++) {
				canal_power.vect[i].mat[0][r] = DVector(filter_length);
				dbRandUni(&cseed, canal_power.vect[i].mat[0][r]);
				canal_power.vect[i].mat[0][r].normalizepower();
				timechannel.vect[i].mat[0][r] = fadingtaps(10, 1 ,filter_length, nb_block>>1, 512, &cseed);
				timechannel.vect[i].mat[1][r] = fadingtaps(10, 1 ,filter_length, nb_block>>1, 512, &cseed);
				canal_power.vect[i].mat[1][r] = DVector(filter_length);
				dbRandUni(&cseed, canal_power.vect[i].mat[1][r]);
				canal_power.vect[i].mat[1][r].normalizepower();
			}
			

		}
			
		ZVector tmp(nb_user_max);
		canal_file_out.open("canal.dat");
	       	sym_file_out.open("symtx.dat");
		cout << "saving symbol and channel info in files ... " << endl;

		for(j=0;j<(nb_block>>1);j++) {
			int offset = 2*j*nb_user_max;
			for(i=0;i<nb_user_max;i++) {
				if(j%2==0) {
					for(r=0;r<nb_rx;r++) {
						tmp.reset(DCplx(0,0));
						for(k=0;k<filter_length;k++) tmp.vect[k] = canal_power.vect[i].mat[0][r].vect[k] * timechannel.vect[i].mat[0][r].mat[j>>1][k];
						canal_file_out << Qh*tmp;
						tmp.reset(DCplx(0,0));
						for(k=0;k<filter_length;k++) tmp.vect[k] =  canal_power.vect[i].mat[1][r].vect[k] * timechannel.vect[i].mat[1][r].mat[j>>1][k];
						canal_file_out << Qh*tmp;
					}
				}		
				sym_file_out << diag(sym_tx.vect[i].copy(offset, offset+nb_user_max))*ucode.vect[i];
				sym_file_out << diag(sym_tx.vect[i].copy(offset+nb_user_max, offset+2*nb_user_max))*ucode.vect[i];
			}			
		}

		for(i=0;i<nb_user_max;i++) {
			timechannel.vect[i].~UMatrix();
			sym_tx.vect[i].~ZVector();
		}
		canal_file_out.close();
		sym_file_out.close();
			
		// os << "****************************" << endl << "timed channel " << endl << "************************" << endl;
			// os << timechannel;
		cout << "Simulating SNR ... " << endl;
		ifstream canal_file_in, sym_file_in, bitin_file_in;//, ldpc_in;
		for(int snr=snr_min; snr<=snr_max;snr+=snr_step) {
			zRandGausInit(&nptr, 0, convertsnrtosigma2(snr));
			cout << " je suis au SNR = " << snr << endl << "je dois en faire "<<snr_max<<"!! " << endl;
			cout << " frame = " << frame<< "  out of " << cin_nb_frame << endl;
			cout << "ME TOUCHE PAS !! JE TRAVAILLEEEEEUUUHHHHHH !!!" << endl;
			canal_file_in.open("canal.dat");
		       	sym_file_in.open("symtx.dat");
			for(j=0;j<(nb_block>>1);j++) {
				int offset = 2*j*nb_user_max;
				//cout << "Extraction necesary channel and symbol info from file ... " << endl;
				if(j%2==0) {
					ZVector tmp1, tmp2; //(nb_user_max);
					for(i=0;i<nb_user_max;i++) {
						for(r=0;r<nb_rx;r++) {
							canal_file_in >> tmp1 >> tmp2;
							canal.vect[i].mat[0][r] = diag(tmp1);
							canal.vect[i].mat[1][r] = diag(tmp2);
						}
					}
				}
						// os << "****************************" << endl << "canal freq at j="<<j << endl << "************************" << endl;
						// os << canal;

		//		sum.reset(DCplx(0,0));

				for(i=0;i<nb_user_max;i++) {
					sym_file_in >> usend1.vect[i] >> usend2.vect[i];
				}

						// os << "****************************" << endl << "spreaded bloc at j="<<j << endl << "************************" << endl;
						// os << " usend 1 " << endl << usend1 << " usend 1 " << endl << usend2;

				//cout << "Simulation starts ... " << endl;
				for(r=0;r<nb_rx;r++) {
					sum1.vect[r].reset(DCplx(0,0));
					sum2.vect[r].reset(DCplx(0,0));
					for(i=0;i<nb_user_max;i++) {
						sum1.vect[r]+=canal.vect[i].mat[0][r]*usend1.vect[i]-canal.vect[i].mat[1][r]*usend2.vect[i].conj();
						sum2.vect[r]+=canal.vect[i].mat[0][r]*usend2.vect[i]+canal.vect[i].mat[1][r]*usend1.vect[i].conj();
					}
				}

						// os << "****************************" << endl << "sum bloc at j="<<j << endl << "************************" << endl;
						// os << " sum1 " << endl << sum1 << " sum2 " << endl << sum2;


				noise = ZMatrix(nb_user_max, nb_user_max);
				zmRandGaus(&nptr, noise);

				ZMatrix decoded1 = (sum1.vect[0]+noise)*c0transpose;
				//ZMatrix decoded1 = (sum1.vect[0])*c0transpose;
				zmRandGaus(&nptr, noise);

				ZMatrix decoded2 = (sum2.vect[0]+noise)*c0transpose;
				//ZMatrix decoded2 = (sum2.vect[0])*c0transpose;

						// os << "****************************" << endl << "decoded at j="<<j << endl << "************************" << endl;
						// os << " decoded1 " << endl << decoded1 << " decoded2 " << endl << decoded2;

				for(i=0;i<nb_user_max;i++) {
								// os << "*-*-*-*-*-*-*  -**-*-*-* *-*-*- \n working on user " << i << endl;
					for(k=0;k<nb_user_max;k++) {
						ZVector rk = ZVector(2);
						rk.vect[0] = decoded1.mat[k][coderep.mat[k][i]];
						rk.vect[1] = decoded2.mat[k][coderep.mat[k][i]].conj();
										// os << " rk = " << rk;
						ZMatrix hki = ZMatrix(2,2);
						hki.mat[0][0]=canal.vect[i].mat[0][0].mat[k][k];
						hki.mat[0][1]=-canal.vect[i].mat[1][0].mat[k][k];
						hki.mat[1][0]=-hki.mat[0][1].conj();
						hki.mat[1][1]=hki.mat[0][0].conj();
										// os << "+ matrice equiv stbc classique " << hki;

						double hpow = hki.mat[0][0].mod2()+ hki.mat[0][1].mod2();
						ZVector mlleft = hki.hermit()*rk;

						// detection ml stbc
						ZVector symest(2);
						double min = 1e20;
						int e1, e2;
						for(int p1 =0; p1<constel.size;p1++) {
							for(int p2=0;p2<constel.size; p2++) {
								symest.vect[0]=hpow*constel.sym[p1];
								symest.vect[1]=hpow*constel.sym[p2].conj();
								ZVector diff = mlleft-symest;
								double dist = diff.vect[0].mod2()+diff.vect[1].mod2();
								if(dist<min) {
									e1 = p1;
									e2 = p2;
									min = dist;
								}
							}
						}
						// end detection
						sym_rx.vect[i].vect[k+offset] = constel.sym[e1];
						sym_rx.vect[i].vect[k+offset+nb_user_max] = constel.sym[e2];
					}
				}
			}
			canal_file_in.close();
			sym_file_in.close();

		//os.close();

		//	os.open("txrx.txt");

			cout << " Decoding and calculating TEB " << endl; 
			bitin_file_in.open("bitin_file.dat");
		       	//ldpc_in.open("ldpcin_file.dat");
			WVector bitin_tmp;
			//WLDPCParityCheck ldpc_tmp;
			int x=0, teb;
			for(i=0;i<nb_user_max;i++)  {
				
#ifdef use_ldpc
				//ldpc_in >> ldpc_tmp;
				WVector tmp = ldpc.vect[i].decode((*constel.demodulate_function)(sym_rx.vect[i]));
#else
				WVector tmp = (*constel.demodulate_function)(sym_rx.vect[i]);
#endif
				bitin_file_in>>bitin_tmp;
				teb =tmp.diff_count(bitin_tmp);
				overall_teb.mat[i][snr]+=teb;
				overall_sym.mat[i][snr]+=tmp.taille;
				x+=teb;

				// os << "User " << i << " : " << sym_rx.vect[i]<<sym_tx.vect[i] <<endl;
				os << "frame = " << frame << ", snr = " << snr << ", user " << i << " diff : " << teb << " / " << bitin_tmp.taille<<endl;//(*constel.demodulate_function)(sym_tx.vect[i]) << endl;
				cout << "frame = " << frame << ", snr = " << snr << ", user " << i << " diff : " << teb << " / " << bitin_tmp.taille<<endl;//(*constel.demodulate_function)(sym_tx.vect[i]) << endl;
				// cout << tmp - bit_in.vect[i]<<endl;
			}
			bitin_file_in.close();
			//ldpc_in.close();
			if(x==0) break;
		}
	}
		os << " ********* OVERALL ****************************** " <<endl;
	for(i=0;i<nb_user_max;i++) {
		for(int snr=snr_min; snr<snr_max; snr+=snr_step) {
				os << "snr = " << snr << ", user " << i << " diff : " << overall_teb.mat[i][snr] << " / " << overall_sym.mat[i][snr]<<endl;
				cout << " snr = " << snr << ", user " << i << " diff : " << overall_teb.mat[i][snr] << " / " << overall_sym.mat[i][snr]<<endl;
		}
	}
	os << " ****************************************** " <<endl;
	os.close();
/*
        char toto;
        cin >> toto;
*/
	return 1;
}

#else
int main() {
	int nb_user_max = 16;
	int snr_min,snr_max, snr_step;
	int cin_nb_block, cin_nb_frame;
#define automat
#ifndef automat
	cout << "Enter nb of user :"; cout.flush();
	cin >> nb_user_max;

	cout << "Enter nb of frame :"; cout.flush();
	cin >> cin_nb_frame;
	cout << "Enter frame size:"; cout.flush();
	cin >> cin_nb_block;
	cout << "Enter snrmin:"; cout.flush();
	cin >> snr_min;
	cout << "Enter snrmax:"; cout.flush();
	cin >> snr_max;
	cout << "Enter snrstep:"; cout.flush();
	cin >> snr_step;
#else
	nb_user_max=8;
	cin_nb_frame=20;
	cin_nb_block=10;
	snr_min=0;
	snr_max=40;
	snr_step=3;
#endif
	UVector<DVector> code(nb_user_max);
	int i, j, k;//, l, t;
	WVector wtmp(nb_user_max);
	UVector<WVector> bit_in(nb_user_max);
	UVector<ZVector> sym_tx(nb_user_max);
	TConstellation constel = psk8gray;
	dRandUniStatePtr cseed;
	int initseed = time(NULL);
	dRandUniInit(&cseed, initseed, 0,1);

	ZMatrix Qh = ZMatrix::Q_DFT(nb_user_max).hermit();

	// chaque ligne represente la sequence d'etalement pour une sous porteuse
	UVector<DMatrix> ucode(nb_user_max);
	UVector<ZVector> sym_rx(nb_user_max);
//#define use_ldpc
	int ldpc_m =1024 , ldpc_n=2048;
#ifdef use_ldpc
	UVector<WLDPCParityCheck> ldpc(nb_user_max);
#endif

	TRandBitStatePtr seed;

	char filename[5000];
	sprintf(filename, "mccdma-res_%d-%ld.txt", nb_user_max, initseed);
	ofstream os(filename, ios::app);
	int nb_block = 2*cin_nb_block*ldpc_m;
	os << "Time : " << initseed << endl << "nb user max = " << nb_user_max << endl << "nb block = " << nb_block << endl<<" Frame : " << cin_nb_frame <<endl;
	os << "Snr min : " << snr_min << endl << "snr_max : " << snr_max << endl << "snr_step : "<< snr_step<<endl;
	cout << "ME TOUCHE PAS !! JE TRAVAILLEEEEEUUUHHHHHH !!!" << endl;

#ifdef use_ldpc
	os << "using ldpc .... " << endl;
#else
	os << "not using ldpc .... " << endl;
#endif

	WMatrix overall_teb(nb_user_max, snr_max), overall_sym(nb_user_max, snr_max);
	ofstream ldpc_out_file;
	//ldpc_out_file.open("ldpcin_file.dat");

	for(i=0;i==0&&i<nb_user_max; i++) {
		code.vect[i] = dWalshCode(nb_user_max, i);
		wtmp.vect[i]=i;
#ifdef use_ldpc
		ldpc.vect[i]=WLDPCParityCheck(ldpc_m,ldpc_n);
		ldpc.vect[i].make_dense_mixed();
#endif
	}
	
	WMatrix coderep = ColumnCirculant(wtmp);
	
	UVector<ZMatrix> Pu(nb_user_max);

	// en ligne = frequence
	// en colonne = utilisateur
	// ex : coderep.mat[3][4] = code d'etalement a utiliser pour utilisateur 4 sur la sous porteuse 3

	for(i=0;i<nb_user_max;i++) {
		ucode.vect[i]=DMatrix(nb_user_max, nb_user_max);
		Pu.vect[i]=ZMatrix(nb_user_max, nb_user_max);
		for(j=0;j<nb_user_max;j++) { // ligne ==> la porteuse
			for(k=0;k<nb_user_max;k++) { // colonne == chip de la porteuse
				ucode.vect[i].mat[j][k] = code.vect[coderep.mat[j][i]].vect[k];
				Pu.vect[i].mat[j][coderep.mat[j][i]] = 1.0;
			}
		}
	}
	
	// ouff ... on a la matrice des code de chaque user
	zRandGausStatePtr nptr;
	typedef UMatrix<ZMatrix> UMZMatrix;
	typedef UVector<UMZMatrix> TSTBCChannel;
	
	int nb_rx = 1, r=0;
	
	ZMatrix c0transpose = ucode.vect[0].transpose_nip();
	UVector<ZMatrix> usend1(nb_user_max), usend2(nb_user_max);
	UVector<ZMatrix> sum1(nb_rx), sum2(nb_rx);
	ofstream bitin_file,/* ldpc_out_file,*/ sym_file_out, canal_file_out;
	
	for(int frame=0;frame<cin_nb_frame;frame++) {
		

		cout << "generating and channel coding user symbols " << endl;

		bitin_file.open("bitin_file.dat");
	       	//ldpc_out_file.open("ldpcin_file.dat");
		for(i=0;i==0&&i<nb_user_max; i++) {
			WVector bit_in_tmp(nb_block*nb_user_max*3);//(int)floor(log2(constel.size))); // deux blocs == un mot astbcm
			bRandBit1(&seed, bit_in.vect[i]);
#ifdef use_ldpc
#pragma message("using ldpc in "__FILE__)
//			cout << "ldpc encoding for user " << i << " with " << bit_in_tmp.taille << " data size" << endl;
			cout << i<<endl;
			sym_tx.vect[i] = (*constel.modulate_function)(ldpc.vect[i].encode(bit_in_tmp));
			//ldpc_out_file << ldpc.vect[i];
#else
#pragma message("not using ldpc in "__FILE__)
			sym_tx.vect[i] = (*constel.modulate_function)(bit_in_tmp);
#endif
			sym_rx.vect[i] = ZVector(sym_tx.vect[i].taille);
			bitin_file << bit_in_tmp;
		}

		bitin_file.close();
		//ldpc_out_file.close();
		//ldpc.~UVector();
		
		// etalement de tt le monde

		ZMatrix noise;

		for(r=0;r<nb_rx;r++) {
			sum1.vect[r]=ZMatrix(nb_user_max, nb_user_max);
			sum2.vect[r]=ZMatrix(nb_user_max, nb_user_max);
		}


		TSTBCChannel timechannel(nb_user_max), canal(nb_user_max);


		int filter_length = 3;
							
		//ZVector tmp(nb_user_max);

		for(i=0;i==0&&i<nb_user_max;i++) {
			timechannel.vect[i]=UMZMatrix(2, nb_rx);
			canal.vect[i]=UMZMatrix(2, nb_rx);
			cout << "Generating user " << i << " channels ... " << endl;
			for(r=0;r<nb_rx;r++) {
				timechannel.vect[i].mat[0][r] = fadingtaps(10, 1 ,filter_length, nb_block>>1, 512, &cseed);
				timechannel.vect[i].mat[1][r] = fadingtaps(10, 1 ,filter_length, nb_block>>1, 512, &cseed);
			}
		}
			
		ZVector tmp(nb_user_max);
		canal_file_out.open("canal.dat");
	       	sym_file_out.open("symtx.dat");
		cout << "saving symbol and channel info in files ... " << endl;

		for(j=0;j<(nb_block>>1);j++) {
			int offset = 2*j*nb_user_max;
			for(i=0;i==0&&i<nb_user_max;i++) {
				if(j%2==0) {
					for(r=0;r<nb_rx;r++) {
						tmp.reset(DCplx(0,0));
						for(k=0;k<filter_length;k++) tmp.vect[k] = timechannel.vect[i].mat[0][r].mat[j>>1][k];
						canal_file_out << Qh*tmp;
						tmp.reset(DCplx(0,0));
						for(k=0;k<filter_length;k++) tmp.vect[k] = timechannel.vect[i].mat[1][r].mat[j>>1][k];
						canal_file_out << Qh*tmp;
					}
				}		
				sym_file_out << diag(sym_tx.vect[i].copy(offset, offset+nb_user_max))*ucode.vect[i];
				sym_file_out << diag(sym_tx.vect[i].copy(offset+nb_user_max, offset+2*nb_user_max))*ucode.vect[i];
			}			
		}

		for(i=0;i==0&&i<nb_user_max;i++) {
			timechannel.vect[i].~UMatrix();
			sym_tx.vect[i].~ZVector();
		}
		canal_file_out.close();
		sym_file_out.close();
			
		// os << "****************************" << endl << "timed channel " << endl << "************************" << endl;
			// os << timechannel;
		cout << "Simulating SNR ... " << endl;
		ifstream canal_file_in, sym_file_in, bitin_file_in;//, ldpc_in;
		for(int snr=snr_min; snr<=snr_max;snr+=snr_step) {
			zRandGausInit(&nptr, 0, convertsnrtosigma2(snr));
			cout << " je suis au SNR = " << snr << endl << "je dois en faire "<<snr_max<<"!! " << endl;
			cout << " frame = " << frame<< "  out of " << cin_nb_frame << endl;
			cout << "ME TOUCHE PAS !! JE TRAVAILLEEEEEUUUHHHHHH !!!" << endl;
			canal_file_in.open("canal.dat");
		       	sym_file_in.open("symtx.dat");
			for(j=0;j<(nb_block>>1);j++) {
				int offset = 2*j*nb_user_max;
				//cout << "Extraction necesary channel and symbol info from file ... " << endl;
				if(j%2==0) {
					ZVector tmp1, tmp2; //(nb_user_max);
					for(i=0;i==0&&i<nb_user_max;i++) {
						for(r=0;r<nb_rx;r++) {
							canal_file_in >> tmp1 >> tmp2;
							canal.vect[i].mat[0][r] = diag(tmp1);
							canal.vect[i].mat[1][r] = diag(tmp2);
						}
					}
				}
						// os << "****************************" << endl << "canal freq at j="<<j << endl << "************************" << endl;
						// os << canal;

		//		sum.reset(DCplx(0,0));

				for(i=0;i==0&&i<nb_user_max;i++) {
					sym_file_in >> usend1.vect[i] >> usend2.vect[i];
				}

						// os << "****************************" << endl << "spreaded bloc at j="<<j << endl << "************************" << endl;
						// os << " usend 1 " << endl << usend1 << " usend 1 " << endl << usend2;

				//cout << "Simulation starts ... " << endl;
				for(r=0;r<nb_rx;r++) {
					sum1.vect[r].reset(DCplx(0,0));
					sum2.vect[r].reset(DCplx(0,0));
					for(i=0;i==0&&i<nb_user_max;i++) {
						sum1.vect[r]+=canal.vect[i].mat[0][r]*usend1.vect[i]-canal.vect[i].mat[1][r]*usend2.vect[i].conj();
						sum2.vect[r]+=canal.vect[i].mat[0][r]*usend2.vect[i]+canal.vect[i].mat[1][r]*usend1.vect[i].conj();
					}
				}

						// os << "****************************" << endl << "sum bloc at j="<<j << endl << "************************" << endl;
						// os << " sum1 " << endl << sum1 << " sum2 " << endl << sum2;
				
				noise = ZMatrix(nb_user_max, nb_user_max);
				zmRandGaus(&nptr, noise);

				ZMatrix decoded1 = (sum1.vect[0]+noise)*c0transpose;
				zmRandGaus(&nptr, noise);

				ZMatrix decoded2 = (sum2.vect[0]+noise)*c0transpose;

						// os << "****************************" << endl << "decoded at j="<<j << endl << "************************" << endl;
						// os << " decoded1 " << endl << decoded1 << " decoded2 " << endl << decoded2;

				for(i=0;i==0&&i<nb_user_max;i++) {
								// os << "*-*-*-*-*-*-*  -**-*-*-* *-*-*- \n working on user " << i << endl;
					for(k=0;k<nb_user_max;k++) {
						ZVector rk = ZVector(2);
						rk.vect[0] = decoded1.mat[k][coderep.mat[k][i]];
						rk.vect[1] = decoded2.mat[k][coderep.mat[k][i]].conj();
										// os << " rk = " << rk;
						ZMatrix hki = ZMatrix(2,2);
						hki.mat[0][0]=canal.vect[i].mat[0][0].mat[k][k];
						hki.mat[0][1]=-canal.vect[i].mat[1][0].mat[k][k];
						hki.mat[1][0]=-hki.mat[0][1].conj();
						hki.mat[1][1]=hki.mat[0][0].conj();
										// os << "+ matrice equiv stbc classique " << hki;

						double hpow = hki.mat[0][0].mod2()+ hki.mat[0][1].mod2();
						ZVector mlleft = hki.hermit()*rk;

						// detection ml stbc
						ZVector symest(2);
						double min = 1e20;
						int e1, e2;
						for(int p1 =0; p1<constel.size;p1++) {
							for(int p2=0;p2<constel.size; p2++) {
								symest.vect[0]=hpow*constel.sym[p1];
								symest.vect[1]=hpow*constel.sym[p2].conj();
								ZVector diff = mlleft-symest;
								double dist = diff.vect[0].mod2()+diff.vect[1].mod2();
								if(dist<min) {
									e1 = p1;
									e2 = p2;
									min = dist;
								}
							}
						}
						// end detection
						sym_rx.vect[i].vect[k+offset] = constel.sym[e1];
						sym_rx.vect[i].vect[k+offset+nb_user_max] = constel.sym[e2];
					}
				}
			}
			canal_file_in.close();
			sym_file_in.close();

		//os.close();

		//	os.open("txrx.txt");

			cout << " Decoding and calculating TEB " << endl; 
			bitin_file_in.open("bitin_file.dat");
		       	//ldpc_in.open("ldpcin_file.dat");
			WVector bitin_tmp;
			//WLDPCParityCheck ldpc_tmp;
			int x=0, teb;
			for(i=0;i==0&&i<nb_user_max;i++)  {
				
#ifdef use_ldpc
				//ldpc_in >> ldpc_tmp;
				WVector tmp = ldpc.vect[i].decode((*constel.demodulate_function)(sym_rx.vect[i]));
#else
				WVector tmp = (*constel.demodulate_function)(sym_rx.vect[i]);
#endif
				bitin_file_in>>bitin_tmp;
				teb =tmp.diff_count(bitin_tmp);
				overall_teb.mat[i][snr]+=teb;
				overall_sym.mat[i][snr]+=tmp.taille;
				x+=teb;

				// os << "User " << i << " : " << sym_rx.vect[i]<<sym_tx.vect[i] <<endl;
				os << "frame = " << frame << ", snr = " << snr << ", user " << i << " diff : " << teb << " / " << bitin_tmp.taille<<endl;//(*constel.demodulate_function)(sym_tx.vect[i]) << endl;
				cout << "frame = " << frame << ", snr = " << snr << ", user " << i << " diff : " << teb << " / " << bitin_tmp.taille<<endl;//(*constel.demodulate_function)(sym_tx.vect[i]) << endl;
				// cout << tmp - bit_in.vect[i]<<endl;
			}
			bitin_file_in.close();
			//ldpc_in.close();
			if(x==0) break;
		}
	}
		os << " ********* OVERALL ****************************** " <<endl;
	for(i=0;i==0&&i<nb_user_max;i++) {
		for(int snr=snr_min; snr<snr_max; snr+=snr_step) {
				os << "snr = " << snr << ", user " << i << " diff : " << overall_teb.mat[i][snr] << " / " << overall_sym.mat[i][snr]<<endl;
				cout << " snr = " << snr << ", user " << i << " diff : " << overall_teb.mat[i][snr] << " / " << overall_sym.mat[i][snr]<<endl;
		}
	}
	os << " ****************************************** " <<endl;
	os.close();
/*
        char toto;
        cin >> toto;
*/
	return 1;
}

#endif
#else
int main() {
	WLDPCParityCheck ldpc(1024,2048);
	TRandBitStatePtr seed;
	zRandGausStatePtr nseed;
	make_method meth=Dense;
	
	WVector bit_in(psk8gray.bit_per_sym*1000*1024);
	bRandBit1(&seed, bit_in);
	ldpc.make_dense_mixed(meth);
	cout << "modulating ... " << endl;
	ZVector sym = (*psk8gray.modulate_function)(ldpc.encode(bit_in));
	ZVector noise(sym.taille);
	cout << "snr simulation ... " << endl;
	for(int snr=0;snr<40;snr+=3) {
		cout << "snr = " << snr << endl;
		zRandGausInit(&nseed, 0, convertsnrtosigma2(1.0/6.0,snr));
		zbRandGaus(&nseed, noise);
		cout << "decoding .... " << endl;
		WVector bitrx = ldpc.decode((*psk8gray.demodulate_function)(noise+sym));
		int diffcnt = bitrx.diff_count(bit_in);
		cout << "Snr = " << snr << " : err = " << diffcnt << " / " << bitrx.taille << endl;
		if(diffcnt==0) break;
	}
	return 0;
	
}
#endif

