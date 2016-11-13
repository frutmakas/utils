#include <iostream.h>
#include <stdlib.h>
#include <modulation/psk.h>
#include <coding/block/ldpc.h>
#include <tools/all.h>
#ifndef WIN_DOS
#include <sys/time.h>
#else
#include <time.h>
#endif

int main(void) {
    make_method method=Evencol;
    int M=1000, N=2000;
    char buffer[5000];
    sprintf(buffer, "ldpc_%d_%d_%d.txt", M, N, 3);
	WLDPCParityCheck ldpc(M,N,3,method);
	method= Dense;
	ldpc.make_dense_mixed(method,1);
    ofstream os(buffer);
    os << ldpc;
    os.close();
    return 0;
}
int maimn(void) {
    WVector bit_in(80000);
    TRandBitStatePtr seed;
    seed.iseed=time(NULL);
    bRandBit1(&seed, bit_in);
    make_method method=Evencol;
	WLDPCParityCheck ldpc(200,400,3,method);
	method= Dense;
	ldpc.make_dense_mixed(method,1);
	WVector btmp = ldpc.encode(bit_in);
    for (int t=0;t<20000;t++) {
        int xx=int(rand()/(RAND_MAX+1.0)*btmp.taille);
        cout << "toggle " << xx << endl;
        btmp.vect[xx]=btmp.vect[xx]==0?1:0;
    }
	DVector ldpc_rx=btmp.RZtoNRZ_nip();
	cout << "Status :> Decoding LDPC sequence (size = "<< bit_in.taille<<")" << endl;
	DVector llr;
	WVector bin_dec=ldpc.decode(ldpc_rx, llr);
	cout << bin_dec.diff_count(bit_in) << " errors out of " << bin_dec.taille << endl;
#ifndef	WIN_DOS
	getchar();
#endif
	return 1;
}

