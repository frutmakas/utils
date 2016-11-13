/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/07/19 20:46:58 $
 * $Revision: 1.1.2.11 $
 * $Id: sfbcofdmldpc.cpp,v 1.1.2.11 2004/07/19 20:46:58 syed Exp $
 *******************************************************************/
 
#include <iostream.h>
#include <stdlib.h>
#include <modulation/psk.h>
#include <coding/block/ldpc.h>
#include <tools/all.h>
#include <rand/randgen.h>
#include <cdma/walsh.h>
#include <fading/fadingtdma.h>
#include <coding/block/ldpc.h>
#include <time.h>

// WARNING: UTILITIS IO "REPROGRAMMED"

//#define nonoisefortest

#define DONT_USE_LDPC	
#define DONT_USE_INTERLEAVER
#define DONT_SAVE_RAW_DATA

int main() {
	const int nb_user=1, nb_user_max=2, nb_porteuse=128, nb_frame=16; 
	const int snr_min=0, snr_max=25, snr_step=3;
	const int interleaver_size=16000;
	const int flmax=5, flmin=4, flstep=3;
	double max_doppler=75;
	double symbol_time=64;
	int LDPC_N = 400, LDPC_M=200;
	int nb_rays = 1024;



	int timestamp=time(NULL)/*1751147*/, filter_length=2;
	char bfname[500];
	sprintf(bfname, "%d-sfbc-simul_%d_%d.log", timestamp, nb_porteuse, flmin);
	ofstream logfile(bfname, ios::app);
#ifndef DONT_USE_LDPC

#define ldpcfileout
    make_method emethod=Evencol, dmethod=Dense;

#ifdef ldpcfileout

	WLDPCParityCheck xldpc(LDPC_M, LDPC_N, 3, emethod);
    xldpc.make_dense_mixed(dmethod,1);
    ofstream ldpc_file("ldpc_code.dat");
    ldpc_file<<xldpc;
    ldpc_file.close();
#else
	WLDPCParityCheck xldpc;
    ifstream ldpc_file("ldpc_code.dat");
    ldpc_file>>xldpc;
    ldpc_file.close();
#endif	

#endif // dont use LDPC

dRandUniStatePtr tapseed;
dRandUniInit(&tapseed, timestamp, 0,1);
TRandBitStatePtr bseed;
bseed.iseed=timestamp;
#ifndef DONT_USE_INTERLEAVER
UVector<CInterleaver> interleaver(nb_user);
for(int ui=0;ui<nb_user;ui++) {
    interleaver.vect[ui]=CInterleaver(interleaver_size);
}
#endif
for(filter_length=flmin;filter_length<=flmax;filter_length+=flstep) {
    int frame_size=5;
    cout << "filter length  = " << filter_length << endl;
	int nb_samples = frame_size;
	logfile<< "Timestamp = " << timestamp<< endl;
	logfile<< "nb user : " << nb_user << endl;
	logfile<< "code length : " << nb_user_max << endl;
	logfile<< "no porteuse : " << nb_porteuse << endl;
	logfile<< "filter length : " << filter_length << endl;
	logfile<< "frame_size : " << frame_size << endl;
	logfile<< "max_doppler : " << max_doppler << endl; //=50;
	logfile<< "symbol_time : " << symbol_time << endl; //=1;
	logfile<< "nb_rays : " << nb_rays << endl; // = 1024;
	logfile<< "nb frame : " << nb_frame << endl; // = 1024;
	logfile << "snr min : " << snr_min << endl;
	logfile << "snr max : " << snr_max << endl;
	logfile << "snr step : " << snr_step << endl;
#ifndef DONT_USE_INTERLEAVER
	logfile<< "interleaver size : " << interleaver_size << endl;
#endif
#ifndef DONT_USE_LDPC
	logfile<< "LDPC size : " << LDPC_M <<"x"<<LDPC_N << endl;
#endif

	WMatrix wer(nb_user, snr_max);

	if(nb_porteuse%2) {
		cerr << "Error : nb _porteuse doit etre multiple de 2" << endl;
		return -1;
	}
	WMatrix odiff(nb_user, snr_max), bdiff(nb_user, snr_max);
	WVector nb_bit_total(snr_max);

	double coef = 1.0;

//	ZMatrix Qh = ZMatrix::Q_IDFT(nb_porteuse, filter_length, coef);
	ZMatrix Qh = ZMatrix::Q_DFT(nb_porteuse, filter_length, coef);
	char varname[5000];
	int u,k,l;
	DMatrix walsh(nb_user_max, nb_user_max);
	DVector dtmp;
	UVector<DMatrix> code(nb_user_max);
	for(u=0;u<nb_user_max;u++) {
		code.vect[u] = DMatrix(nb_porteuse, nb_user_max);
		dtmp = dWalshCode(nb_user_max, u);
		for(int l=0;l<nb_user_max;l++) {
			walsh.mat[u][l]=dtmp.vect[l];
			for(int k=0;k<nb_porteuse;k++) {
				code.vect[u].mat[k][l]=dtmp.vect[l];
			}
		}
	}
	logfile<< "ouch .. !" << endl;

	WMatrix nb_bit_simulated(nb_user, snr_max);
	TConstellation constel = psk8;
	ZVector xxxtemp(constel.size);
	for(int xxx=0;xxx<constel.size;xxx++) {
	   xxxtemp.vect[xxx]=constel.sym[xxx];
    }
    //WVector wwwtemp((*constel.demodulate_function)(xxxtemp));
    DMatrix dddtemp = LinePackedVectorAsMatrix((*constel.demodulate_function)(xxxtemp).RZtoNRZ_nip(), constel.size, constel.bit_per_sym);
    UVector<DVector> symmap = SplitMatrixByLine(dddtemp);
    int interleaver_sub_size=interleaver_size, isize=0;
	// balance
	if((nb_porteuse*frame_size*(LDPC_N-LDPC_M))%interleaver_size==0) {
		interleaver_sub_size=1;
	}
	
    if(interleaver_size%nb_porteuse==0) {
		isize = interleaver_size/nb_porteuse;
		if (interleaver_sub_size>isize) interleaver_sub_size=isize;
	}
    if(interleaver_size%frame_size==0) {
		isize= interleaver_size/frame_size;
		if (interleaver_sub_size>isize) interleaver_sub_size=isize;
	} 
    if(interleaver_size%(LDPC_N-LDPC_M)==0) {
		isize = interleaver_size/(LDPC_N-LDPC_M);
		if (interleaver_sub_size>isize) interleaver_sub_size=isize;
	}
 
	
	// end balance try
	int nb_bit = constel.bit_per_sym*nb_porteuse*frame_size*(LDPC_N-LDPC_M)*interleaver_sub_size;	
	logfile<< "interleaver sub size : " << interleaver_sub_size << endl;
	logfile<< "nb_bit : " << nb_bit<<endl;
	cout << "interleaver sub size : " << interleaver_sub_size << endl;
	cout << "nb_bit : " << nb_bit<<endl;
	
	double sqrt_05=sqrt(0.5);

	typedef UVector<ZMatrix> TTapsPath;
	typedef UVector<DVector> TTapsPower;
	UVector<TTapsPath> taps(nb_user);
	UVector<TTapsPower> tpower(nb_user);
    zRandGausStatePtr nptr;
	ofstream tapsofile[nb_user][2];
    ifstream tapsifile[nb_user][2];
    ofstream dataofile[nb_user];
    ifstream dataifile[nb_user];
    ofstream bitrxofile[snr_max][nb_user];
    ifstream bitrxifile[snr_max][nb_user];

    int sym_size;
    /*
#ifndef  DONT_USE_LDPC
        nb_samples = nb_bit*2*nb_frame/(constel.bit_per_sym*nb_porteuse*2);
#else
        nb_samples = nb_bit*nb_frame/(constel.bit_per_sym*nb_porteuse*2);
#endif
//   */
    nb_samples = nb_bit*2*nb_frame/(constel.bit_per_sym*nb_porteuse);
    logfile<< "total nb of samples " << nb_samples << endl;
    for(u=0;u<nb_user;u++) {
    	tpower.vect[u]=TTapsPower(2);
        taps.vect[u]=TTapsPath(2);
        for(int p=0;p<2;p++) {
				cout << "Creating channel "<<p<< " for user " << u << endl;
				tpower.vect[u].vect[p]=DVector(filter_length,1);
//				dbRandUni(&tapseed, tpower.vect[u].vect[p]);
				tpower.vect[u].vect[p].normalizepower();
				logfile << "user " << u << ", path " << p << " : " << tpower.vect[u].vect[p] << endl;
                ZMatrix tapstemp=fadingtaps(max_doppler,symbol_time,filter_length,nb_samples,nb_rays, &tapseed, tpower.vect[u].vect[p], 20);
				sprintf(varname, "%d-path_U%dP%d-fr%d-l%d.dat", timestamp, u,p, nb_frame,filter_length);
				tapsofile[u][p].open(varname);
				cout << "saving channel in file .. " <<  varname << endl;
				for (int t=0;t<nb_samples;t++) {
				    tapsofile[u][p] << filter_length << " : [ ";
				    for(int l=0;l<filter_length;l++) {
                        tapsofile[u][p] << tapstemp.mat[t][l] << " ";
                    }
                    tapsofile[u][p] << "]" << endl;
                }
                tapsofile[u][p].close();
         }
    }
    cout << "opening channel files .. " << endl;
    for(u=0;u<nb_user;u++) {
        for(int p=0;p<2;p++) {
            sprintf(varname, "%d-path_U%dP%d-fr%d-l%d.dat", timestamp, u,p, nb_frame,filter_length);
            tapsifile[u][p].open(varname);
        }
    }    
    
	for (int frame=0;frame<nb_frame;frame++) {
	    logfile<< "Processing frame " << frame << endl;
        sprintf(varname, "%d-bitin_fr%d-l%d.dat", timestamp, frame, filter_length);
	    ofstream bitofile(varname);
	    for(u=0;u<nb_user;u++) {
            int snr_idx=0;
            double snr=snr_min;
   	        for(snr=snr_min, snr_idx=0; snr<snr_max; snr+=snr_step, snr_idx++) {
		       sprintf(varname, "%d_bitrx_fr%d_snr%d_u%d_l%d.dat", timestamp, frame, snr_idx, u, filter_length);
    		   bitrxofile[snr_idx][u].open(varname);
#ifndef DONT_USE_LDPC
    		   bitrxofile[snr_idx][u] << nb_bit*2/*LDPC_RATE*/ << " : [ "; 
#else
               bitrxofile[snr_idx][u] << nb_bit << " : [ "; 
#endif    		   
    		   bitrxofile[snr_idx][u].flush();
            }
		}
	    
		for (u=0;u<nb_user;u++) {
            WVector bitin=WVector(nb_bit);
			bRandBit1(&bseed, bitin);
#ifndef DONT_USE_LDPC

#ifndef DONT_USE_INTERLEAVER
			ZVector symtx = interleaver.vect[u].Apply((*constel.modulate_function)(xldpc.encode(bitin))*sqrt_05);
#else //DONT_USE_INTERLEAVER
			ZVector symtx = (*constel.modulate_function)(xldpc.encode(bitin))*sqrt_05;
#endif  // DONT_USE_INTERLEAVER

#else // DONT_USE_LDPC

#ifndef DONT_USE_INTERLEAVER
			ZVector symtx = interleaver.vect[u].Apply((*constel.modulate_function)(bitin)*sqrt_05);
#else // DONT_USE_INTERLEAVER
            ZVector symtx = (*constel.modulate_function)(bitin)*sqrt_05;
#endif // DONT_USE_INTERLEAVER

#endif // DONT_USE_LDPC
			frame_size=symtx.taille/nb_porteuse;
			sym_size=symtx.taille;
			sprintf(varname, "%d-symtx_U%d-fr%d-l%d.dat", timestamp, u, frame,filter_length);
			dataofile[u].open(varname);
			cerr << "saving psk_sym in file .. " << endl;
			for (int t=0;t<frame_size;t++) {
    		    dataofile[u] << nb_porteuse << " : [ ";
    		    int xptr=t*nb_porteuse;
			    for(int l=0;l<nb_porteuse;l++) {
                       dataofile[u] << symtx.vect[xptr+l] << " ";
                }
                dataofile[u] << "]" << endl;
            }
            dataofile[u].close();
            bitofile << bitin << endl;
        }
        bitofile.close();
        for(u=0;u<nb_user;u++) {
			sprintf(varname, "%d-symtx_U%d-fr%d-l%d.dat", timestamp, u, frame, filter_length);
			dataifile[u].open(varname);
		}
 
//		int sptr = 0;
		ZVector currenttaps(filter_length);
		ZMatrix rx1(nb_porteuse, nb_user_max), rxu1(nb_porteuse, nb_user_max);//,	rx2=0;
		typedef UVector<ZVector> TSymRx;
		UVector<TSymRx> sym_rx(snr_max);
		WMatrix sptr(snr_max, nb_user,0);
		for(int s=snr_min, si=0;s<snr_max;si++, s+=snr_step) {
			sym_rx.vect[si] = TSymRx(nb_user);
			for(u=0;u<nb_user;u++) {
				sym_rx.vect[si].vect[u] = ZVector(interleaver_size);
			}
		}

        cout << "starting transmission and stressing " << endl;
		for(int f=0;f<frame_size;f++) {
			ZMatrix S1, S2;
			UVector<TTapsPath>L(nb_user);
			ZVector sym_read;
			for (u=0;u<nb_user;u++) {
				L.vect[u]=TTapsPath(2);
				dataifile[u] >> sym_read;
				S1 = diag(sym_read)*code.vect[u];
				S2 = ZMatrix(S1.line, S1.col);
				// permute ligne S2... 
				for (k=0;k<nb_porteuse;k+=2) {
					for(l=0;l<S1.col;l++) {
						S2.mat[k][l]   = -S1.mat[k+1][l].conj();
						S2.mat[k+1][l] =  S1.mat[k][l].conj();
					}
				}
				for (int p=0;p<2;p++) {
					tapsifile[u][p] >> currenttaps;
					L.vect[u].vect[p]=diag(Qh*currenttaps);
				}
				rxu1 = L.vect[u].vect[0]*S1+L.vect[u].vect[1]*S2;
				rx1=(u!=0) ? rx1+rxu1:rxu1;
			} //tranmission end

			ZMatrix noise(rx1.line, rx1.col), nrx1;
			int snr_idx=0;
			double snr=snr_min;
			for(snr=snr_min, snr_idx=0;snr<snr_max;snr+=snr_step, snr_idx++) {
#ifndef nonoisefortest
				zRandGausInit(&nptr, 0, convertsnrtosigma2(snr));
    			zmRandGaus(&nptr, noise);
				nrx1 = rx1+noise;
#else
                nrx1=rx1;
#endif
				ZMatrix y1 = nrx1*walsh.transpose_nip();
				//% on a les messages stbc sur chaque colonne par user
				for (u=0;u<nb_user;u++) {
					ZMatrix zy =ExtractColXFromMatrix(y1, u);
					for (k=1;k<nb_porteuse;k+=2) {
						 zy.mat[k][0]=zy.mat[k][0].conj();
					}
					ZMatrix z(2,1);

					for (k=0;k<nb_porteuse;k+=2) {
						z.mat[0][0] = L.vect[u].vect[0].mat[k+1][k+1].conj()*zy.mat[k][0]+
									  L.vect[u].vect[1].mat[k][k]*zy.mat[k+1][0];
						z.mat[1][0] = L.vect[u].vect[0].mat[k][k]*zy.mat[k+1][0]-
									  L.vect[u].vect[1].mat[k+1][k+1].conj()*zy.mat[k][0];
						DCplx Lh = L.vect[u].vect[0].mat[k+1][k+1].conj()*L.vect[u].vect[0].mat[k][k]+
							       L.vect[u].vect[1].mat[k+1][k+1].conj()*L.vect[u].vect[1].mat[k][k];
                    
						double diffmin1=1e20, diff1;
						double diffmin2=1e20, diff2;
						int e1, e2, p1;

						for(p1=0;p1<constel.size;p1++) {
							diff1 = (z.mat[0][0]-Lh*constel.sym[p1]).mag();
							diff2 = (z.mat[1][0]-Lh*constel.sym[p1].conj()).mag();
							if (diff1 < diffmin1) {
								diffmin1 = diff1;
								e1 = p1;		
							}
							if (diff2 < diffmin2) { 
								diffmin2 = diff2;
								e2 = p1;
							}
						}
						sym_rx.vect[snr_idx].vect[u].vect[sptr.mat[snr_idx][u]++]=constel.sym[e1];
						if (sptr.mat[snr_idx][u]>=interleaver_size) {
#ifndef DONT_USE_INTERLEAVER
                          ZVector xirx = interleaver.vect[u].Extract(sym_rx.vect[snr_idx].vect[u]);
#else
                          ZVector xirx = sym_rx.vect[snr_idx].vect[u];
#endif

#ifndef DONT_USE_LDPC
						  DVector wtmp = (*constel.demodulate_function)(xirx).RZtoNRZ_nip();
#else
						  WVector wtmp = (*constel.demodulate_function)(xirx);
#endif
						  sptr.mat[snr_idx][u]=0;
                          for (int w=0;w<wtmp.taille;w++) bitrxofile[snr_idx][u] << wtmp.vect[w] << " ";
		                }
		                sym_rx.vect[snr_idx].vect[u].vect[sptr.mat[snr_idx][u]++]=constel.sym[e2];
						if (sptr.mat[snr_idx][u]>=interleaver_size) {
#ifndef DONT_USE_INTERLEAVER
                          ZVector xirx = interleaver.vect[u].Extract(sym_rx.vect[snr_idx].vect[u]);
#else
                          ZVector xirx = sym_rx.vect[snr_idx].vect[u];
#endif
#ifndef DONT_USE_LDPC
						  DVector wtmp = (*constel.demodulate_function)(xirx).RZtoNRZ_nip();
#else
						  WVector wtmp = (*constel.demodulate_function)(xirx);
#endif
						  sptr.mat[snr_idx][u]=0;
                          for (int w=0;w<wtmp.taille;w++) bitrxofile[snr_idx][u] << wtmp.vect[w] << " ";
		                }
					}
				}
			} //snr
//			sptr=sptr+nb_porteuse;        
		}//frame_size
		int osnr_idx=0;
		double osnr=snr_min;
		for(osnr=snr_min, osnr_idx=0;osnr<snr_max;osnr+=snr_step, osnr_idx++) {
		  for(u=0;u<nb_user;u++) {
		    bitrxofile[osnr_idx][u] << "]" <<endl;
    		bitrxofile[osnr_idx][u].close();
    		sprintf(varname, "%d_bitrx_fr%d_snr%d_u%d_l%d.dat", timestamp, frame, osnr_idx, u, filter_length);
    		bitrxifile[osnr_idx][u].open(varname);
		  }
		}
		cout << "End processing frame " << frame << endl;
		for(u=0;u<nb_user;u++) {        
          dataifile[u].close(); 
#ifdef DONT_SAVE_RAW_DATA
          sprintf(varname, "%d-symtx_U%d-fr%d-l%d.dat", timestamp, u, frame, filter_length);
          remove(varname);
#endif          
        }
//		logfile<< "sptr =  " << sptr << endl;
		logfile<< "Rframesize=" << frame_size << endl;
		
		int snr_idx=0;
		double snr=snr_min;
    	bdiff.reset(0);
 	    cout << "Decoding and  evaluating performance for differnet SNR" << frame << endl;
 		sprintf(varname, "%d-bitin_fr%d-l%d.dat", timestamp, frame, filter_length);
 	    ifstream bitifile(varname);
	    
		for (u=0;u<nb_user;u++) {
		      cout << "#"; cout.flush();
		      WVector bit_in_read;
		      bitifile >> bit_in_read;
			for(snr=snr_min, snr_idx=0;snr<snr_max; snr+=snr_step, snr_idx++) {
			     cout << "*"; cout.flush();
                //DVector dout = (*constel.demodulate_function)(sym_rx.vect[snr_idx].vect[u]).RZtoNRZ_nip();				
#ifndef DONT_USE_LDPC
                DVector dout; 
                bitrxifile[snr_idx][u] >> dout;
            	DVector llr;
            	WVector bout=xldpc.decode(dout, llr);
				bdiff.mat[u][snr_idx]=bout.diff_count(bit_in_read);
				nb_bit_simulated.mat[u][snr_idx]+=bout.taille;
#else
                WVector bout; 
                bitrxifile[snr_idx][u] >> bout;
				bdiff.mat[u][snr_idx]=bout.diff_count(bit_in_read);
				nb_bit_simulated.mat[u][snr_idx]+=bout.taille;
#endif
			}
			cout<< endl;
		}
		bitifile.close();
#ifdef DONT_SAVE_RAW_DATA
        sprintf(varname, "%d-bitin_fr%d-l%d.dat", timestamp, frame, filter_length);
        remove(varname);
#endif		
		for(snr=snr_min, snr_idx=0;snr<snr_max;snr+=snr_step, snr_idx++) {
		  for(u=0;u<nb_user;u++) {
    		bitrxifile[snr_idx][u].close();
#ifdef DONT_SAVE_RAW_DATA    		
    		sprintf(varname, "%d_bitrx_fr%d_snr%d_u%d_l%d.dat", timestamp, frame, snr_idx, u, filter_length);
    		remove(varname);
#endif    		
		  }
		}
		odiff=odiff+bdiff;
		logfile<< " intermediate odiff at " << frame << endl <<odiff<<endl;
		cout << " intermediate odiff at " << frame << endl <<odiff<<nb_bit_simulated<<endl;
	}
    for(u=0;u<nb_user;u++) {
         for(int p=0;p<2;p++) {
             tapsifile[u][p].close();
 		}
	}

	logfile<< "ODIFF = " << odiff << endl << "NB_bit_sim = " << nb_bit_simulated << endl;	
    cout << "ODIFF = " << odiff << endl << "NB_bit_sim = " << nb_bit_simulated << endl;	
}
logfile.close();
cout << "simulation done .. . " <<endl;
	getchar();
	char toto[78];
	cin >> toto;
	return 1;
}
