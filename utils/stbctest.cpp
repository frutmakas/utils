/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/08/04 12:23:46 $
 * $Revision: 1.1.2.1 $
 * $Id: stbctest.cpp,v 1.1.2.1 2004/08/04 12:23:46 syed Exp $
 ********************************************************************/
 
#include <iostream>
using namespace std;
#include <all.h>
#include <time.h>
//#define nonoise
int main(){
    const int rxmax_c = 3, rxmin=1, rxstep=1, nspb=75000, nbblock=500;
    const double snr_min=0, snr_max=25,snr_step=2;
    TConstellation constel = psk8;
    int nbpb = constel.bit_per_sym*nspb;
    WVector bit(nbpb);
    double sqrt_1_2 = 1.0/sqrt(2.0);
    for(int s=0;s<constel.size;s++) constel.sym[s]*=sqrt_1_2;
    int timestamp=time(NULL);
    TRandBitStatePtr bseed;
    bseed.iseed = time(NULL);
    zRandGausStatePtr zptr;
    dRandUniStatePtr dptr;
    dRandUniInit(&dptr, time(NULL));
	for(int rxmax=rxmin;rxmax<=rxmax_c;rxmax+=rxstep) {
/*		UMatrix<ZMatrix> taps(2, rxmax);
		for(int tx=0;tx<2;tx++) {
			for (int rx=0; rx<rxmax;rx++) {
				cout << " Creating channel .. " << tx << " <-> " << rx << endl;
				//taps.mat[tx][rx] = ZMatrix((nbblock*nspb)>>1, 1, DCplx(tx,rx));//fadingtaps(75, 1, 1, (nbblock*nspb)>>1, 1024, &dptr);
				taps.mat[tx][rx] = fadingtaps(75, 1, 1, (nbblock*nspb)>>1, 1024, &dptr);
				//cout << taps.mat[tx][rx];
			}    
		}   */
	/*    int t=0; 
		int errsum=0;*/
		int cptr = 0, cptr2=0;
		WVector errcnt(100);
		for(int f=0;f<nbblock;f++) {
       		UMatrix<ZMatrix> taps(2, rxmax);
    		for(int tx=0;tx<2;tx++) {
    			for (int rx=0; rx<rxmax;rx++) {
    				cout << " Creating channel .. " << tx << " <-> " << rx << endl;
    				//taps.mat[tx][rx] = ZMatrix((nbblock*nspb)>>1, 1, DCplx(tx,rx));//fadingtaps(75, 1, 1, (nbblock*nspb)>>1, 1024, &dptr);
    				taps.mat[tx][rx] = fadingtaps(75, 1, 1, (nspb)>>1, 1024, &dptr);
    				//cout << taps.mat[tx][rx];
    			}    
    		}   
		    cptr = cptr2=0;
			cout << "+ frame  : " << f << endl;
			bRandBit1(&bseed, bit);
			ZVector symtx = (*constel.modulate_function)(bit);
			//cout << symtx << endl;
			ZMatrix stbctx(nspb, 2);
			cout << "+ + creating stbc symbols .. " << endl;
			int i;
			for(i=0;i<symtx.taille;i+=2) {
				stbctx.mat[i][0]=symtx.vect[i];
				stbctx.mat[i+1][0]=symtx.vect[i+1];
				stbctx.mat[i][1]=-symtx.vect[i+1].conj();
				stbctx.mat[i+1][1]=symtx.vect[i].conj();
			}    
			//cout << stbctx ;
			ZMatrix stbcrx(stbctx.line,rxmax);
			cout << "+ + combining stbc symbols .. " << endl;

			for(i=0;i<stbctx.line;i++) {
				for(int r=0;r<rxmax;r++) {
					stbcrx.mat[i][r]=0;
					for(int tx=0;tx<2;tx++) {
						stbcrx.mat[i][r]+=stbctx.mat[i][tx]*taps.mat[tx][r].mat[cptr][0];
					}
				}
				if (i&1) cptr++;

			}            
			cout << "+ + adding noise .." << endl;
			ZMatrix noise(stbcrx.line, rxmax);
			int cptr3 = cptr2;
			int snr_idx=0;
			for(double snr=snr_min; snr<=snr_max;snr+=snr_step,snr_idx++) {
				cptr2 = cptr3;
				cout << "+ + + " << snr <<"dB .. " << endl;
				zRandGausInit(&zptr, 0, convertsnrtosigma2(snr));
				zmRandGaus(&zptr, noise);
	#ifndef nonoise
				 ZMatrix rxsym = noise+stbcrx;
	#else
				 ZMatrix rxsym = stbcrx;
	#endif
				DVector modh(rxmax);
				ZVector r1(rxmax),r2(rxmax), symrx(stbcrx.line);
				for(int i=0;i<stbcrx.line;i+=2) {
					double min1=1e20, min2=1e20, e1,e2;
					int p1,p2;
					for(int r=0;r<rxmax;r++) {
						modh.vect[r] = taps.mat[0][r].mat[cptr2][0].mod2()+taps.mat[1][r].mat[cptr2][0].mod2();
						r1.vect[r] = taps.mat[0][r].mat[cptr2][0].conj()*rxsym.mat[i][r]+taps.mat[1][r].mat[cptr2][0]*rxsym.mat[i+1][r].conj();
						r2.vect[r] = taps.mat[0][r].mat[cptr2][0]*rxsym.mat[i+1][r].conj()-taps.mat[1][r].mat[cptr2][0].conj()*rxsym.mat[i][r];
					}    
					for(int p=0;p<constel.size;p++) {
						e1=e2=0;
						for(int r=0;r<rxmax;r++) {
							e1 += (r1.vect[r]-modh.vect[r]*constel.sym[p]).mag();
							e2 += (r2.vect[r]-modh.vect[r]*constel.sym[p].conj()).mag();
						}    
						if(e1<min1) {
							p1 = p;
							min1=e1;
						}    
						if(e2<min2) {
							p2 = p;
							min2=e2;
						}    
					}    
					symrx.vect[i] = constel.sym[p1];
					symrx.vect[i+1] = constel.sym[p2];
					cptr2++;
				}    
				//cout << symrx - symtx<< endl;
				WVector bitrx = (*constel.demodulate_function)(symrx);
				//cout << bitrx - bit << endl;
				int err = bitrx.diff_count(bit);
				errcnt.vect[snr_idx]+=err;
				cout << "+ + + + " << err << " errors " << endl;
				if (err == 0) break;
			}    
			cout << "+ + end block "<< f << endl;
		}
		int snr_idx=0;
		char varname[450];
		sprintf(varname, "%d_psk8res_rx%d.log", timestamp, rxmax); 
		ofstream os(varname, ios::app);
		for(double snr=snr_min; snr<=snr_max;snr+=snr_step,snr_idx++) {
			os << rxmax << ";" << snr<< ";" <<errcnt.vect[snr_idx] << ";" << nspb*nbblock << endl;
		}
		os.close();
	}
	system("pause");
    return 1;
}        
        
        
            

            
        
                                    
        
        
        
        
    
    
