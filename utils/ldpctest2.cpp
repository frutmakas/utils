/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/08/05 20:08:29 $
 * $Revision: 1.1.2.3 $
 * $Id: ldpctest2.cpp,v 1.1.2.3 2004/08/05 20:08:29 syed Exp $
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

int main1(TConstellation &constel, int M, int N, int base, int kf) {
#define bypasschannel    
//#define noldpc
//    int M=200, N=400;
//    TConstellation constel = bpsk;
    int k = base/constel.bit_per_sym;
    WLDPCParityCheck ldpc(M,N,3);
    ldpc.make_dense_mixed();
    TRandBitStatePtr bseed;
    bseed.iseed=time(NULL);
    zRandGausStatePtr zseed;
    dRandUniStatePtr dseed;
    dRandUniInit(&dseed, time(NULL));

    int nbbit = constel.bit_per_sym*(N-M)*k;
    WVector bitin(nbbit);
    bRandBit1(&bseed, bitin);
    cout << "created " << nbbit << " bits .. " << endl;
    char buffer[5000];
    int timestamp = time(NULL);
    sprintf(buffer, "taps%ld.m", timestamp);
    ofstream mo(buffer);
#ifndef noldpc
    ZVector psktx = (*constel.modulate_function)(ldpc.encode(bitin));
#else
    ZVector psktx = (*constel.modulate_function)(bitin);
#endif    
    cout << "created modulated symbols .. " << endl;
    int max_doppler=75, symbol_time_us=10, nbray=1024;
#ifndef bypasschannel
    ZMatrix taps = fadingtaps(max_doppler, symbol_time_us, 1, psktx.taille, nbray, &dseed);
    cout << "created channel taps .. " << endl;
//ofstream &matlaboutput(ofstream &os, const char* varname, const ZMatrix &m, int precision) {
    matlaboutput(mo, buffer, taps, 15);
    mo << "figure; semilogy(abs("<<buffer<<").^2);" << endl;
    mo.close();
#else
    cout << "bypassing channel ... " << endl;
//    ZMatrix taps(psktx.taille, 1, DCplx(1,0)); //fadingtaps(max_doppler, symbol_time_us, 1, psktx.taille, nbray, &dseed);
#endif
    ZVector noise(psktx.taille);
    ZVector pskrxne(psktx.taille), pskrxe(psktx.taille);
    ofstream os("res.txt", ios::app);
    os << " ************** "<< constel.id <<" ********************* " << endl;
    os << " timestamp = " << timestamp << endl;
    os << " doppler " << max_doppler << endl;
    os << " symtime " << symbol_time_us << endl;
#ifdef bypasschannel
    os << "channel : pure WGN     " << endl;
#endif  
#ifdef noldpc
    os << " no ldpc used .. " << endl;
#else 
    os << "LDPC " << M << "x" << N << endl;    
#endif
    double snrmin=-5, snrmax=15, snrstep=0.5;
    bool out= false;
    for(double snr=snrmin;snr<=snrmax && out==false;snr+=snrstep) {
#ifndef noldpc
        zRandGausInit(&zseed, 0, convertsnrtosigma2(double(N-M)/(M*constel.bit_per_sym), snr));
#else
        zRandGausInit(&zseed, 0, convertsnrtosigma2(1.0/(constel.bit_per_sym), snr));
#endif        
        zbRandGaus(&zseed, noise);
        for(register int x=0; x< psktx.taille;x++) {
#ifndef bypasschannel            
//            pskrxne.vect[x]=psktx.vect[x]*taps.mat[x][0]+noise.vect[x];
            pskrxe.vect[x]=pskrxne.vect[x]*taps.mat[x][0].conj()/taps.mat[x][0].mod2();
#else
            pskrxe.vect[x]=psktx.vect[x]+noise.vect[x];
#endif
        }    
#ifndef noldpc
  //      WVector bitrxne = ldpc.decode((*constel.demodulate_function)(pskrxne).RZtoNRZ_nip());
        WVector bitrxe = ldpc.decode((*constel.demodulate_function)(pskrxe).RZtoNRZ_nip());
#else
//        WVector bitrxne = (*constel.demodulate_function)(pskrxne);
        WVector bitrxe = (*constel.demodulate_function)(pskrxe);
#endif        
        int diff=bitrxe.diff_count(bitin);
        out = diff==0;
        cout << "SNR=" << snr << endl;
        os  << "E SNR = " << snr << ", error = " << diff << " out of " << bitrxe.taille << "( TEB= " << (double)(diff)/bitrxe.taille << ")" << endl;
/*        diff=bitrxe.diff_count(bitin);
        out = out && (diff==0);
        os  << "E  SNR = " << snr << ", error = " << diff << " out of " << bitrxe.taille << "( TEB= " << (double)(diff)/bitrxe.taille << ")" << endl;
  */
        if (out) cout << "no more errors" <<endl;
    }
    os.close();
    cout << "DONE .. ." <<endl;
//    getchar();
    return 1;
}    
        
    
    
int main() {
  /*  cout << "BPSK ... " << endl;
    main1(bpsk);
    cout << "QPSK ... " << endl;
    main1(psk4);
    cout << "QPSK Gray... " << endl;
    main1(psk4gray);*/
    cout << "8PSK 200x400... " << endl;
    int M=200, N=400;
    int nbbit=24000;
    main1(psk8,M,N, nbbit);
    cout << "8PSK Gray 200x400... " << endl;
    main1(psk8gray,M,N, nbbit);
    M=2000;N=4000;
    cout << "BPSK 2000x4000 ... " << endl;
    nbbit=2400;
    main1(bpsk,M,N, nbbit);
    cout << "QPSK 2000x4000 ... " << endl;
    main1(psk4,M,N, nbbit);
    cout << "QPSK Gray 2000x4000 ... " << endl;
    main1(psk4gray,M,N, nbbit);    
    cout << "8PSK 2000x4000... " << endl;
    main1(psk8,M,N, nbbit);
    cout << "8PSK Gray 2000x4000... " << endl;
    main1(psk8gray,M,N, nbbit);
        getchar();
    return 1;
}
    
