<<<<<<< MIMOStbc.cpp
#include <all.h>
#include <fstream.h>
#include <time.h>
//#include <exception.h>

#define NOT_DECOMPOSABLE  1
#define DECOMPOSED        0

#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 

int CreateChannelFilter(const ZMatrix& HQ, ZMatrix &PQ, ZMatrix& VQ, ZMatrix &W) { 
	//checked with matlab
	// status OK .. potential problem with singular matrix
	if(HQ.line<=0 || HQ.col<=0 || HQ.mat==NULL) throw CUtilitisException(EMPTY_MATRIX);
	if(HQ.line!=HQ.col) throw CUtilitisException(INVALID_MATRIX_DIMENSION);
	if (HQ.line==2) {
           VQ=HQ;
           W=ZMatrix::eye(2);
           return NOT_DECOMPOSABLE;
        }

	ZMatrix A = SubMatrix(HQ, 0,0, HQ.line-2, HQ.col-2);
	ZMatrix B = SubMatrix(HQ, 0, HQ.col-2, HQ.line-2, 2);
	ZMatrix G = SubMatrix(HQ, HQ.line-2,0, 2, HQ.col-2);
	ZMatrix D = SubMatrix(HQ, HQ.line-2, HQ.col-2,2,2);
	
	//cout << HQ << A << B << G << D;

	W = ZMatrix::eye(HQ.line);
	ZMatrix W2 = -B*D.inv_nip();
	InsertSubMatrixIntoMatrix(W, W2, 0, HQ.col-2);
	ZMatrix W3 = -G*A.inv_nip();
	InsertSubMatrixIntoMatrix(W, W3, HQ.line-2, 0);
//	wtemp.~ZMatrix();

	PQ=A+W2*G;
	VQ=W3*B+D;

	//cout << W2 << W3;

	return DECOMPOSED;
}

int decoupleUsers(const ZMatrix YQ, const ZMatrix &HQ, ZMatrix &YQp, ZMatrix &YQv, 
				  ZMatrix &PQ, ZMatrix &VQ) {
	if (YQ.line<=0 || YQ.col<=0 || YQ.mat==NULL ||
		HQ.line<=0 || HQ.col<=0 || HQ.mat==NULL) throw CUtilitisException(EMPTY_MATRIX);
	if(YQ.line!=HQ.col) throw CUtilitisException(INCOMPATIBLE_MATRIX_DIMENSION);
	if(HQ.line!=HQ.col) CUtilitisException(INVALID_MATRIX_DIMENSION);

	ZMatrix W;
	if (CreateChannelFilter(HQ, PQ, VQ, W)==NOT_DECOMPOSABLE) {
           YQv=YQ;
           return NOT_DECOMPOSABLE;
        }
	ZMatrix newYQ = W*YQ;
	YQp=SubMatrix(newYQ, 0, 0, YQ.line-2, 1);
	YQv=SubMatrix(newYQ, YQ.line-2, 0, 2, 1);
	return DECOMPOSED;
}

ZVector decodeASTBC(const UVector<ZVector> &r, const UMatrix<ZMatrix> &H, const TConstellation &constel) {
	int blocksize = r.vect[0].taille;
	int nb_rx=r.taille;
	int i, t, frame;
	ZMatrix xest(2,1), h(2,2), rt(2,1), xdet(2,1);
	ZVector X(blocksize);
	for(t=0, frame=0;t<blocksize;t+=2, frame++) {
		double ht_min=1e20;
		for (int p1=0; p1<constel.size;p1++) {
			for (int p2=0; p2<constel.size;p2++) {
				xest.mat[0][0]=constel.sym[p1];
				xest.mat[1][0]=constel.sym[p2];
				double mlsum=0;
				for (i=0; i < nb_rx; i++) {
					h.mat[0][0]=H.mat[0][i].mat[frame][0];//(stbc_frame, :, i);
					h.mat[0][1]=H.mat[1][i].mat[frame][0];//(stbc_frame, :, i);
					h.mat[1][0]=-h.mat[0][1].conj();
					h.mat[1][1]=h.mat[0][0].conj();
					double rho = h.mat[0][0].mod2()+h.mat[0][1].mod2();
					rt.mat[0][0]=r.vect[i].vect[t];
					rt.mat[1][0]=r.vect[i].vect[t+1].conj();
					DMatrix mlmat=(h.hermit()*rt-rho*xest).mag();
					mlsum += ColSum(mlmat).mat[0][0];
				}
				if(mlsum<ht_min) {
					ht_min=mlsum;
					xdet=xest;
				}
			}
		}
		X.vect[t]=xdet.mat[0][0];
		X.vect[t+1]=xdet.mat[1][0];
	}
	return X;
}

ZMatrix decodePreparedSingleSTBCBlock(const UVector<ZMatrix> &r, const UVector<ZMatrix> &h, const TConstellation &constel) {
	double ht_min=1e20;
	ZMatrix xest(2,1), xdet;
	int  nb_rx = r.taille, i;
	for (int p1=0; p1<constel.size;p1++) {
		for (int p2=0; p2<constel.size;p2++) {
			xest.mat[0][0]=constel.sym[p1];
			xest.mat[1][0]=constel.sym[p2];
			double mlsum=0;
			for (i=0; i < nb_rx; i++) {
				double rho = h.vect[i].mat[0][0].mod2()+h.vect[i].mat[0][1].mod2();
				DMatrix mlmat=(h.vect[i].hermit()*r.vect[i]-rho*xest).mag();
				mlsum += ColSum(mlmat).mat[0][0];
			}
			if(mlsum<ht_min) {
				ht_min=mlsum;
				xdet=xest;
			}
		}
	}
	return  xdet;
}

UVector<ZVector> propagateASTBCInFadingChannel(const ZMatrix &X, int nb_rx, double SNR, UMatrix<ZMatrix> &H, dRandUniStatePtr *tapseed) {
	H = UMatrix<ZMatrix>(2, nb_rx);
	UVector<ZVector> r(nb_rx);
	double doppler = 10, symbolT_us =1;
	int blocksize = X.line;
	int i, j, frame;
	for(i=0;i<2;i++) {
		for(j=0;j<nb_rx; j++) {
			H.mat[i][j]=fadingtaps(doppler, symbolT_us, 1, blocksize>>1, 512, tapseed);
		}
	}
	
	ZVector noise(blocksize);
	zRandGausStatePtr noise_ptr;
	zRandGausInit(&noise_ptr ,0, convertsnrtosigma2(1, SNR));

	for(i=0;i<nb_rx; i++) {
		r.vect[i] = ZVector(blocksize);
		zbRandGaus(&noise_ptr, noise);
		for(j=0, frame=0; j<blocksize; j+=2, frame++) {
			r.vect[i].vect[j]=H.mat[0][i].mat[frame][0]*X.mat[j][0]+H.mat[1][i].mat[frame][0]*X.mat[j][1]+noise.vect[j];
			r.vect[i].vect[j+1]=H.mat[0][i].mat[frame][0]*X.mat[j+1][0]+H.mat[1][i].mat[frame][0]*X.mat[j+1][1]+noise.vect[j+1];
		}
	}
	return r;
}


void singleusermain(void) {
	ZMatrix h(8,8);
/*	for(int i=0;i<8;i++) {
		for(int j=0;j<8;j++) {
			h.mat[i][j]=i*8+j;
		}
	}*/
	zRandGausStatePtr ptr;
	zRandGausInit(&ptr, 14, 30);
	zmRandGaus(&ptr, h);
	ofstream matlab("d:\\MATLAB6p5\\work\\fromC.m");
	matlaboutput(matlab, "h", h, 15);

	ZMatrix p,v,w;
	CreateChannelFilter(h, p, v, w);
	matlaboutput(matlab, "cp", p, 15);
	matlaboutput(matlab, "cv", v, 15);
	matlaboutput(matlab, "cw", w, 15);
	matlab.close();
}

//void multiusermain(void)  {
void main(void)  {
	cout << "Compiled on " <<  __TIMESTAMP__ << endl;

	int  nb_user_min=3, nb_user_max=3, nb_user_step=4, Nu;
	double snr_min=100, snr_max=100, snr_step=10;

	int nb_packet=100;
	int optimal_frame_size=100, i;
	typedef UMatrix<ZMatrix> UMatrixZMatrix;
	TRandBitStatePtr seed;
	dRandUniStatePtr tapseed;
	dRandUniInit(&tapseed, time(NULL), 0, 1);

	cout << "How many packets ?"; cout.flush();
	cin >> nb_packet;

	cout << "optimal frame size?"; cout.flush();
	cin >> optimal_frame_size;

	cout << "How many user min?"; cout.flush();
	cin >> nb_user_min;

	cout << "How many user max?"; cout.flush();
	cin >> nb_user_max;

	cout << "How many user step?"; cout.flush();
	cin >> nb_user_step;

	cout << "snr min ?"; cout.flush();
	cin >> snr_min;

	cout << "snr max ?"; cout.flush();
	cin >> snr_max;

	cout << "snr step ?"; cout.flush();
	cin >> snr_step;


	TConstellation  constellation = psk4gray;
	int bits_per_symbol=(int)log2(constellation.size);
	int bits_generated=0;
	for(int nb_user=nb_user_min;nb_user<=nb_user_max; nb_user+=nb_user_step) {
		for(double SNR=snr_min;SNR<=snr_max;SNR+=snr_step) {
			UVector<WVector> bit_in(nb_user);
			UVector<ZMatrix> Xtx(nb_user);
			UVector <ZVector> Xdet(nb_user);
			UVector<ZVector> r(nb_user);
			UVector<UMatrixZMatrix> H(nb_user);
			WVector TEB(nb_user, 0);
			TEB.reset(0);
			bits_generated=0;
			for(int frame=0;frame<nb_packet<<2;frame+=optimal_frame_size){
				int current_packet_size=((nb_packet<<2)-frame)>optimal_frame_size?optimal_frame_size:((nb_packet<<2)-frame);
				for(Nu=0;Nu<nb_user; Nu++) {
					bit_in.vect[Nu] = WVector(current_packet_size<<2);
					bRandBit1(&seed, bit_in.vect[Nu]);
					Xtx.vect[Nu]=G2_type2((*constellation.modulate_function)(bit_in.vect[Nu]), bit_in.vect[Nu].taille/bits_per_symbol);
					Xdet.vect[Nu]=ZVector(Xtx.vect[Nu].line);
				}
				bits_generated+=bit_in.vect[0].taille;
				for(Nu=0;Nu<nb_user; Nu++) {
					if(Nu==0) {
						r  =propagateASTBCInFadingChannel(Xtx.vect[Nu], nb_user, SNR, H.vect[Nu], &tapseed);
					} else {
						r +=propagateASTBCInFadingChannel(Xtx.vect[Nu], nb_user, SNR, H.vect[Nu], &tapseed);
					}
				}
				// signal propagated .. arrived at rx antennas
				// int decoupleUsers(const ZMatrix YQ, const ZMatrix &HQ, ZMatrix &YQp, ZMatrix &YQv,
				//		  ZMatrix &PQ, ZMatrix &VQ) {
				//ofstream os("debugmimostbc.txt");
				//os << r << Xtx << H << "***************************" <<endl ;
				int htime=0;
				for(i=0;i<Xtx.vect[0].line;i+=2, htime++) {
					//preparing YQ and HQ
					ZMatrix YQ(nb_user*2,1);
					ZMatrix HQ(2*nb_user, 2*nb_user);
					for(Nu=0;Nu<nb_user;Nu++) { // user sweep
						YQ.mat[2*Nu][0]=r.vect[Nu].vect[i];
						YQ.mat[2*Nu+1][0]=r.vect[Nu].vect[i+1].conj();
						for(int nu=0;nu<nb_user; nu++) { // rxantenna sweep
							HQ.mat[2*nu][2*Nu] = H.vect[Nu].mat[0][nu].mat[htime][0];
							HQ.mat[2*nu][2*Nu+1] = H.vect[Nu].mat[1][nu].mat[htime][0];
							HQ.mat[2*nu+1][2*Nu] = -H.vect[Nu].mat[1][nu].mat[htime][0].conj();
							HQ.mat[2*nu+1][2*Nu+1] = H.vect[Nu].mat[0][nu].mat[htime][0].conj();
						}
					}
					// on a HQ et YQ initial .
//					os << YQ << HQ << endl;
//					os.close();
					// on decode les nb_user utilisateurs
					ZMatrix YQp, YQv, PQ, VQ;
					for(Nu=nb_user-1; Nu>=0; Nu--) {
						decoupleUsers(YQ, HQ, YQp, YQv, PQ, VQ);
						double ht_min=1e20;
						ZMatrix xest(2,1), xdet;
						for (int p1=0; p1<constellation.size;p1++) {
							for (int p2=0; p2<constellation.size;p2++) {
								xest.mat[0][0]=constellation.sym[p1];
								xest.mat[1][0]=constellation.sym[p2];
								double rho = VQ.mat[0][0].mod2()+VQ.mat[0][1].mod2();
								DMatrix mlmat=(VQ.hermit()*YQv-rho*xest).mag();
								double mlsum = ColSum(mlmat).mat[0][0];
								if(mlsum<ht_min) {
									ht_min=mlsum;
									xdet=xest;
								}
							}
						}
						//on a xdet = symbol detecté pour user NU
						Xdet.vect[Nu].vect[i]=xdet.mat[0][0];
						Xdet.vect[Nu].vect[i+1]=xdet.mat[1][0];
						// preparation de YQ et HQ pour iteration suivante
						YQ=YQp;
						HQ=PQ;
					}
					// decodage terminé pour temps i
				}
				//decodage du frame terminé
				for(Nu=0;Nu<nb_user;Nu++) {
					TEB.vect[Nu]+=bit_in.vect[Nu].diff_count((*constellation.demodulate_function)(Xdet.vect[Nu]));
				}
			}//frame end
			ofstream osres("mimostbc_fading_result.txt", ios::app);
			int errorsum;
			errorsum=0;
			for(Nu=0;Nu<nb_user;Nu++) {
				errorsum+=TEB.vect[Nu];
				cout << "(system_nb_user, SNR, user #, generated_bits, userN_error, userN_TEB) =" << nb_user <<","<< SNR <<","<< Nu <<","<< bits_generated <<","<< TEB.vect[Nu] << "," << (double)TEB.vect[Nu]/(double)bits_generated << ")" << endl;
				osres << "(system_nb_user, SNR, user #, generated_bits, userN_error, userN_TEB) =" << nb_user <<","<< SNR <<","<< Nu <<","<< bits_generated <<","<< TEB.vect[Nu] << "," << (double)TEB.vect[Nu]/(double)bits_generated << ")" << endl;
			}
			osres.close();
			if(errorsum==0) break;
		} //SNR end
	} //nb_user end
	int dummy;
	//cin >> dummy;
}

void xmain(void) {
	cout << "Compiled on " <<  __TIMESTAMP__ << endl;
	cout << "How many packets ?"; cout.flush();
	int nbpacket;
	cin >> nbpacket;
	cout << "How many rx antenna max?"; cout.flush();
	int nb_rx_max;
	cin >> nb_rx_max;
	cout << "snr min ?"; cout.flush();
	double snr_min;
	cin >> snr_min;
	cout << "snr max ?"; cout.flush();
	double snr_max;
	cin >> snr_max;
	cout << "snr step ?"; cout.flush();
	double snr_step;
	cin >> snr_step;

	ofstream os("resultmimostbc.txt");
	os << " nb packet : " <<nbpacket<<endl;
	os << " nb_rx_max : " <<nb_rx_max <<endl;
	os << " snrmin: " <<snr_min<<endl;
	os << " snrmax: " <<snr_max<<endl;
	os << " snrstep: " <<snr_step<<endl;
	os.precision(15);

	int nb_bit=2*2*nbpacket;
	WVector bit_in(nb_bit);
	TRandBitStatePtr seed;
	bRandBit1(&seed, bit_in);
	//ZVector psk_in = z4PSKGray_mod(bit_in);
	TConstellation constellation = psk4gray;
	ZVector psk_in = (*constellation.modulate_function)(bit_in);
//	WVector bitout= (*constellation.demodulate_function)(psk_in);
	//cout << psk_in << bit_in-bitout << endl; 
	ZMatrix X = G2_type2(psk_in, psk_in.taille);
	psk_in.~ZVector();
	UMatrix<ZMatrix> H;
//	int diff;
	dRandUniStatePtr tapseed;
	dRandUniInit(&tapseed, time(NULL), 0, 1);

	for(int nb_rx=1; nb_rx<=nb_rx_max; nb_rx++) {
		for(double snr=snr_min; snr<=snr_max; snr+=snr_step) {
			//UVector<ZVector> r = propagateASTBCInFadingChannel(X, nb_rx, snr, H);
			WVector bit_out = (*constellation.demodulate_function)(decodeASTBC(propagateASTBCInFadingChannel(X, nb_rx, snr, H, &tapseed), H, constellation));
			int diff=bit_in.diff_count(bit_out);
			cout.precision(5);
			cout << "(nb_rx, snr, diff, teb) = ("<< nb_rx << "," << snr << "," << diff << "," <<  (double)diff/(double)bit_in.taille<<")"<<endl;
			os << "(nb_rx, snr, diff, teb) = ("<< nb_rx << "," << snr << "," << diff << "," <<  (double)diff/(double)bit_in.taille<<")"<<endl;
			if (diff==0) break;
		}
	}
}
=======
#include <all.h>
#include <fstream.h>
#include <time.h>
//#include <exception.h>

#define NOT_DECOMPOSABLE  1
#define DECOMPOSED        0

#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 

int CreateChannelFilter(const ZMatrix& HQ, ZMatrix &PQ, ZMatrix& VQ, ZMatrix &W) { 
	//checked with matlab
	// status OK .. potential problem with singular matrix
	if(HQ.line<=0 || HQ.col<=0 || HQ.mat==NULL) throw CUtilitisException(EMPTY_MATRIX);
	if(HQ.line!=HQ.col) throw CUtilitisException(INVALID_MATRIX_DIMENSION);
	if (HQ.line==2) {
           VQ=HQ;
           W=ZMatrix::eye(2);
           return NOT_DECOMPOSABLE;
        }

	ZMatrix A = SubMatrix(HQ, 0,0, HQ.line-2, HQ.col-2);
	ZMatrix B = SubMatrix(HQ, 0, HQ.col-2, HQ.line-2, 2);
	ZMatrix G = SubMatrix(HQ, HQ.line-2,0, 2, HQ.col-2);
	ZMatrix D = SubMatrix(HQ, HQ.line-2, HQ.col-2,2,2);
	
	//cout << HQ << A << B << G << D;

	W = ZMatrix::eye(HQ.line);
	ZMatrix W2 = -B*D.inv_nip();
	InsertSubMatrixIntoMatrix(W, W2, 0, HQ.col-2);
	ZMatrix W3 = -G*A.inv_nip();
	InsertSubMatrixIntoMatrix(W, W3, HQ.line-2, 0);
//	wtemp.~ZMatrix();

	PQ=A+W2*G;
	VQ=W3*B+D;

	//cout << W2 << W3;

	return DECOMPOSED;
}

int decoupleUsers(const ZMatrix YQ, const ZMatrix &HQ, ZMatrix &YQp, ZMatrix &YQv, 
				  ZMatrix &PQ, ZMatrix &VQ) {
	if (YQ.line<=0 || YQ.col<=0 || YQ.mat==NULL ||
		HQ.line<=0 || HQ.col<=0 || HQ.mat==NULL) throw CUtilitisException(EMPTY_MATRIX);
	if(YQ.line!=HQ.col) throw CUtilitisException(INCOMPATIBLE_MATRIX_DIMENSION);
	if(HQ.line!=HQ.col) CUtilitisException(INVALID_MATRIX_DIMENSION);

	ZMatrix W;
	if (CreateChannelFilter(HQ, PQ, VQ, W)==NOT_DECOMPOSABLE) {
           YQv=YQ;
           return NOT_DECOMPOSABLE;
        }
	ZMatrix newYQ = W*YQ;
	YQp=SubMatrix(newYQ, 0, 0, YQ.line-2, 1);
	YQv=SubMatrix(newYQ, YQ.line-2, 0, 2, 1);
	return DECOMPOSED;
}

ZVector decodeASTBC(const UVector<ZVector> &r, const UMatrix<ZMatrix> &H, const TConstellation &constel) {
	int blocksize = r.vect[0].taille;
	int nb_rx=r.taille;
	int i, t, frame;
	ZMatrix xest(2,1), h(2,2), rt(2,1), xdet(2,1);
	ZVector X(blocksize);
	for(t=0, frame=0;t<blocksize;t+=2, frame++) {
		double ht_min=1e20;
		for (int p1=0; p1<constel.size;p1++) {
			for (int p2=0; p2<constel.size;p2++) {
				xest.mat[0][0]=constel.sym[p1];
				xest.mat[1][0]=constel.sym[p2];
				double mlsum=0;
				for (i=0; i < nb_rx; i++) {
					h.mat[0][0]=H.mat[0][i].mat[frame][0];//(stbc_frame, :, i);
					h.mat[0][1]=H.mat[1][i].mat[frame][0];//(stbc_frame, :, i);
					h.mat[1][0]=-h.mat[0][1].conj();
					h.mat[1][1]=h.mat[0][0].conj();
					double rho = h.mat[0][0].mod2()+h.mat[0][1].mod2();
					rt.mat[0][0]=r.vect[i].vect[t];
					rt.mat[1][0]=r.vect[i].vect[t+1].conj();
					DMatrix mlmat=(h.hermit()*rt-rho*xest).mag();
					mlsum += ColSum(mlmat).mat[0][0];
				}
				if(mlsum<ht_min) {
					ht_min=mlsum;
					xdet=xest;
				}
			}
		}
		X.vect[t]=xdet.mat[0][0];
		X.vect[t+1]=xdet.mat[1][0];
	}
	return X;
}

ZMatrix decodePreparedSingleSTBCBlock(const UVector<ZMatrix> &r, const UVector<ZMatrix> &h, const TConstellation &constel) {
	double ht_min=1e20;
	ZMatrix xest(2,1), xdet;
	int  nb_rx = r.taille, i;
	for (int p1=0; p1<constel.size;p1++) {
		for (int p2=0; p2<constel.size;p2++) {
			xest.mat[0][0]=constel.sym[p1];
			xest.mat[1][0]=constel.sym[p2];
			double mlsum=0;
			for (i=0; i < nb_rx; i++) {
				double rho = h.vect[i].mat[0][0].mod2()+h.vect[i].mat[0][1].mod2();
				DMatrix mlmat=(h.vect[i].hermit()*r.vect[i]-rho*xest).mag();
				mlsum += ColSum(mlmat).mat[0][0];
			}
			if(mlsum<ht_min) {
				ht_min=mlsum;
				xdet=xest;
			}
		}
	}
	return  xdet;
}

UVector<ZVector> propagateASTBCInFadingChannel(const ZMatrix &X, int nb_rx, double SNR, UMatrix<ZMatrix> &H, dRandUniStatePtr *tapseed) {
	H = UMatrix<ZMatrix>(2, nb_rx);
	UVector<ZVector> r(nb_rx);
	double doppler = 10, symbolT_us =1;
	int blocksize = X.line;
	int i, j, frame;
	for(i=0;i<2;i++) {
		for(j=0;j<nb_rx; j++) {
			H.mat[i][j]=fadingtaps(doppler, symbolT_us, 1, blocksize>>1, 512, tapseed);
		}
	}
	
	ZVector noise(blocksize);
	zRandGausStatePtr noise_ptr;
	zRandGausInit(&noise_ptr ,0, convertsnrtosigma2(1, SNR));

	for(i=0;i<nb_rx; i++) {
		r.vect[i] = ZVector(blocksize);
		zbRandGaus(&noise_ptr, noise);
		for(j=0, frame=0; j<blocksize; j+=2, frame++) {
			r.vect[i].vect[j]=H.mat[0][i].mat[frame][0]*X.mat[j][0]+H.mat[1][i].mat[frame][0]*X.mat[j][1]+noise.vect[j];
			r.vect[i].vect[j+1]=H.mat[0][i].mat[frame][0]*X.mat[j+1][0]+H.mat[1][i].mat[frame][0]*X.mat[j+1][1]+noise.vect[j+1];
		}
	}
	return r;
}


void singleusermain(void) {
	ZMatrix h(8,8);
/*	for(int i=0;i<8;i++) {
		for(int j=0;j<8;j++) {
			h.mat[i][j]=i*8+j;
		}
	}*/
	zRandGausStatePtr ptr;
	zRandGausInit(&ptr, 14, 30);
	zmRandGaus(&ptr, h);
	ofstream matlab("d:\\MATLAB6p5\\work\\fromC.m");
	matlaboutput(matlab, "h", h, 15);

	ZMatrix p,v,w;
	CreateChannelFilter(h, p, v, w);
	matlaboutput(matlab, "cp", p, 15);
	matlaboutput(matlab, "cv", v, 15);
	matlaboutput(matlab, "cw", w, 15);
	matlab.close();
}

//void multiusermain(void)  {
void main(void)  {
	cout << "Compiled on " <<  __TIMESTAMP__ << endl;

	int  nb_user_min=3, nb_user_max=3, nb_user_step=4, Nu;
	double snr_min=100, snr_max=100, snr_step=10;

	int nb_packet=100;
	int optimal_frame_size=100, i;
	typedef UMatrix<ZMatrix> UMatrixZMatrix;
	TRandBitStatePtr seed;
	dRandUniStatePtr tapseed;
	dRandUniInit(&tapseed, time(NULL), 0, 1);

	cout << "How many packets ?"; cout.flush();
	cin >> nb_packet;

	cout << "optimal frame size?"; cout.flush();
	cin >> optimal_frame_size;

	cout << "How many user min?"; cout.flush();
	cin >> nb_user_min;

	cout << "How many user max?"; cout.flush();
	cin >> nb_user_max;

	cout << "How many user step?"; cout.flush();
	cin >> nb_user_step;

	cout << "snr min ?"; cout.flush();
	cin >> snr_min;

	cout << "snr max ?"; cout.flush();
	cin >> snr_max;

	cout << "snr step ?"; cout.flush();
	cin >> snr_step;


	TConstellation  constellation = psk4gray;
	int bits_per_symbol=(int)log2(constellation.size);
	int bits_generated=0;
	for(int nb_user=nb_user_min;nb_user<=nb_user_max; nb_user+=nb_user_step) {
		for(double SNR=snr_min;SNR<=snr_max;SNR+=snr_step) {
			UVector<WVector> bit_in(nb_user);
			UVector<ZMatrix> Xtx(nb_user);
			UVector <ZVector> Xdet(nb_user);
			UVector<ZVector> r(nb_user);
			UVector<UMatrixZMatrix> H(nb_user);
			WVector TEB(nb_user, 0);
			TEB.reset(0);
			bits_generated=0;
			for(int frame=0;frame<nb_packet<<2;frame+=optimal_frame_size){
				int current_packet_size=((nb_packet<<2)-frame)>optimal_frame_size?optimal_frame_size:((nb_packet<<2)-frame);
				for(Nu=0;Nu<nb_user; Nu++) {
					bit_in.vect[Nu] = WVector(current_packet_size<<2);
					bRandBit1(&seed, bit_in.vect[Nu]);
					Xtx.vect[Nu]=G2_type2((*constellation.modulate_function)(bit_in.vect[Nu]), bit_in.vect[Nu].taille/bits_per_symbol);
					Xdet.vect[Nu]=ZVector(Xtx.vect[Nu].line);
				}
				bits_generated+=bit_in.vect[0].taille;
				for(Nu=0;Nu<nb_user; Nu++) {
					if(Nu==0) {
						r  =propagateASTBCInFadingChannel(Xtx.vect[Nu], nb_user, SNR, H.vect[Nu], &tapseed);
					} else {
						r +=propagateASTBCInFadingChannel(Xtx.vect[Nu], nb_user, SNR, H.vect[Nu], &tapseed);
					}
				}
				// signal propagated .. arrived at rx antennas
				// int decoupleUsers(const ZMatrix YQ, const ZMatrix &HQ, ZMatrix &YQp, ZMatrix &YQv,
				//		  ZMatrix &PQ, ZMatrix &VQ) {
				//ofstream os("debugmimostbc.txt");
				//os << r << Xtx << H << "***************************" <<endl ;
				int htime=0;
				for(i=0;i<Xtx.vect[0].line;i+=2, htime++) {
					//preparing YQ and HQ
					ZMatrix YQ(nb_user*2,1);
					ZMatrix HQ(2*nb_user, 2*nb_user);
					for(Nu=0;Nu<nb_user;Nu++) { // user sweep
						YQ.mat[2*Nu][0]=r.vect[Nu].vect[i];
						YQ.mat[2*Nu+1][0]=r.vect[Nu].vect[i+1].conj();
						for(int nu=0;nu<nb_user; nu++) { // rxantenna sweep
							HQ.mat[2*nu][2*Nu] = H.vect[Nu].mat[0][nu].mat[htime][0];
							HQ.mat[2*nu][2*Nu+1] = H.vect[Nu].mat[1][nu].mat[htime][0];
							HQ.mat[2*nu+1][2*Nu] = -H.vect[Nu].mat[1][nu].mat[htime][0].conj();
							HQ.mat[2*nu+1][2*Nu+1] = H.vect[Nu].mat[0][nu].mat[htime][0].conj();
						}
					}
					// on a HQ et YQ initial .
//					os << YQ << HQ << endl;
//					os.close();
					// on decode les nb_user utilisateurs
					ZMatrix YQp, YQv, PQ, VQ;
					for(Nu=nb_user-1; Nu>=0; Nu--) {
						decoupleUsers(YQ, HQ, YQp, YQv, PQ, VQ);
						double ht_min=1e20;
						ZMatrix xest(2,1), xdet;
						for (int p1=0; p1<constellation.size;p1++) {
							for (int p2=0; p2<constellation.size;p2++) {
								xest.mat[0][0]=constellation.sym[p1];
								xest.mat[1][0]=constellation.sym[p2];
								double rho = VQ.mat[0][0].mod2()+VQ.mat[0][1].mod2();
								DMatrix mlmat=(VQ.hermit()*YQv-rho*xest).mag();
								double mlsum = ColSum(mlmat).mat[0][0];
								if(mlsum<ht_min) {
									ht_min=mlsum;
									xdet=xest;
								}
							}
						}
						//on a xdet = symbol detecté pour user NU
						Xdet.vect[Nu].vect[i]=xdet.mat[0][0];
						Xdet.vect[Nu].vect[i+1]=xdet.mat[1][0];
						// preparation de YQ et HQ pour iteration suivante
						YQ=YQp;
						HQ=PQ;
					}
					// decodage terminé pour temps i
				}
				//decodage du frame terminé
				for(Nu=0;Nu<nb_user;Nu++) {
					TEB.vect[Nu]+=bit_in.vect[Nu].diff_count((*constellation.demodulate_function)(Xdet.vect[Nu]));
				}
			}//frame end
			ofstream osres("mimostbc_fading_result.txt", ios::app);
			int errorsum;
			errorsum=0;
			for(Nu=0;Nu<nb_user;Nu++) {
				errorsum+=TEB.vect[Nu];
				cout << "(system_nb_user, SNR, user #, generated_bits, userN_error, userN_TEB) =" << nb_user <<","<< SNR <<","<< Nu <<","<< bits_generated <<","<< TEB.vect[Nu] << "," << (double)TEB.vect[Nu]/(double)bits_generated << ")" << endl;
				osres << "(system_nb_user, SNR, user #, generated_bits, userN_error, userN_TEB) =" << nb_user <<","<< SNR <<","<< Nu <<","<< bits_generated <<","<< TEB.vect[Nu] << "," << (double)TEB.vect[Nu]/(double)bits_generated << ")" << endl;
			}
			osres.close();
			if(errorsum==0) break;
		} //SNR end
	} //nb_user end
	int dummy;
	//cin >> dummy;
}

void xmain(void) {
	cout << "Compiled on " <<  __TIMESTAMP__ << endl;
	cout << "How many packets ?"; cout.flush();
	int nbpacket;
	cin >> nbpacket;
	cout << "How many rx antenna max?"; cout.flush();
	int nb_rx_max;
	cin >> nb_rx_max;
	cout << "snr min ?"; cout.flush();
	double snr_min;
	cin >> snr_min;
	cout << "snr max ?"; cout.flush();
	double snr_max;
	cin >> snr_max;
	cout << "snr step ?"; cout.flush();
	double snr_step;
	cin >> snr_step;

	ofstream os("resultmimostbc.txt");
	os << " nb packet : " <<nbpacket<<endl;
	os << " nb_rx_max : " <<nb_rx_max <<endl;
	os << " snrmin: " <<snr_min<<endl;
	os << " snrmax: " <<snr_max<<endl;
	os << " snrstep: " <<snr_step<<endl;
	os.precision(15);

	int nb_bit=2*2*nbpacket;
	WVector bit_in(nb_bit);
	TRandBitStatePtr seed;
	bRandBit1(&seed, bit_in);
	//ZVector psk_in = z4PSKGray_mod(bit_in);
	TConstellation constellation = psk4gray;
	ZVector psk_in = (*constellation.modulate_function)(bit_in);
//	WVector bitout= (*constellation.demodulate_function)(psk_in);
	//cout << psk_in << bit_in-bitout << endl; 
	ZMatrix X = G2_type2(psk_in, psk_in.taille);
	psk_in.~ZVector();
	UMatrix<ZMatrix> H;
//	int diff;
	dRandUniStatePtr tapseed;
	dRandUniInit(&tapseed, time(NULL), 0, 1);

	for(int nb_rx=1; nb_rx<=nb_rx_max; nb_rx++) {
		for(double snr=snr_min; snr<=snr_max; snr+=snr_step) {
			//UVector<ZVector> r = propagateASTBCInFadingChannel(X, nb_rx, snr, H);
			WVector bit_out = (*constellation.demodulate_function)(decodeASTBC(propagateASTBCInFadingChannel(X, nb_rx, snr, H, &tapseed), H, constellation));
			int diff=bit_in.diff_count(bit_out);
			cout.precision(5);
			cout << "(nb_rx, snr, diff, teb) = ("<< nb_rx << "," << snr << "," << diff << "," <<  (double)diff/(double)bit_in.taille<<")"<<endl;
			os << "(nb_rx, snr, diff, teb) = ("<< nb_rx << "," << snr << "," << diff << "," <<  (double)diff/(double)bit_in.taille<<")"<<endl;
			if (diff==0) break;
		}
	}
}
>>>>>>> 1.1.2.7
