/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/08/05 20:08:30 $
 * $Revision: 1.1.2.1 $
 * $Id: ldpctest3.cpp,v 1.1.2.1 2004/08/05 20:08:30 syed Exp $
 ********************************************************************/
#include <iostream>
#include <fstream>
using namespace std;

#include <modulation/psk.h>
#include <coding/block/ldpc.h>
#include <tools/all.h>
#include <fading/fadingtdma.h>
#ifndef WIN_DOS
#include <sys/time.h>
#else
#include <time.h>
#endif
//int main1(TConstellation &constel, int M, int N, int base, int kf, double snrmin=-5, double snrmax=15, double snrstep=0.5);

int main1(TConstellation &constel, int M, int N, int base, int kf, double snrmin=-5, double snrmax=15, double snrstep=0.5){
    int k = base/constel.bit_per_sym;
    WLDPCParityCheck ldpc(M,N,3);
    ldpc.make_dense_mixed();
    TRandBitStatePtr bseed;
    bseed.iseed=time(NULL);
    zRandGausStatePtr zseed;
    dRandUniStatePtr dseed;
    dRandUniInit(&dseed, time(NULL));
    int timestamp = time(NULL);
    int nbbit = constel.bit_per_sym*(N-M)*k;
    ofstream os("res.txt", ios::app);
    os << " ************** "<< constel.id <<" ********************* " << endl;
    os << " timestamp = " << timestamp << endl;
//    os << " doppler " << max_doppler << endl;
//    os << " symtime " << symbol_time_us << endl;
    os << "channel : pure WGN     " << endl;
    os << "LDPC " << M << "x" << N << endl;    
    os << "************************"<<endl;
    os << "xx;snr;diff;cnt;teb"<<endl;
    
    for (int xx=0;xx<kf;xx++) {
        WVector bitin(nbbit);
        bRandBit1(&bseed, bitin);
        cout << "created " << nbbit << " bits .. " << endl;
        ZVector psktx = (*constel.modulate_function)(ldpc.encode(bitin));
        cout << "created modulated symbols .. " << endl;
     //   int max_doppler=75, symbol_time_us=10, nbray=1024;
        cout << "bypassing channel ... " << endl;
        ZVector noise(psktx.taille);
        ZVector pskrxne(psktx.taille), pskrxe(psktx.taille);
        bool out= false;
        for(double snr=snrmin;snr<=snrmax && out==false;snr+=snrstep) {
            zRandGausInit(&zseed, 0, convertsnrtosigma2(double(N-M)/(M*constel.bit_per_sym), snr));
            zbRandGaus(&zseed, noise);
            for(register int x=0; x< psktx.taille;x++) {
                pskrxe.vect[x]=psktx.vect[x]+noise.vect[x];
            }    
            WVector bitrxe = ldpc.decode((*constel.demodulate_function)(pskrxe).RZtoNRZ_nip());
            int diff=bitrxe.diff_count(bitin);
            out = diff==0;
            cout << "@" << xx << " SNR = " << snr << ", error = " << diff << " out of " << bitrxe.taille << "( TEB= " << (double)(diff)/bitrxe.taille << ")" << endl;
            os << xx << ";" << snr << ";" << diff << ";" << bitrxe.taille << ";" << (double)(diff)/bitrxe.taille << endl;
            if (out) cout << "no more errors @" << xx << endl;
        }
    }    
    os.close();
    cout << "DONE .. ." <<endl;
//    getchar();
    return 1;
}    
        
    
    
int main() {
    int M=2000,N=4000;
    cout << "BPSK 2000x4000 ... " << endl;
    int nbbit=24, kf=10;
    main1(bpsk,M,N, nbbit, kf);
    return 1;
}
    
