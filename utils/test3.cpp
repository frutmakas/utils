//#include "tools/all.h"
#include "tools/utilitis.h"
#include "tools/tools.h"
#include "modulation/psk.h"
#include "modulation/ofdm.h"
#include "coding/block/ldpc.h"
#include "rand/randgen.h"

int main(void){
	WVector input(3*256);
	TRandBitStatePtr s;
	bRandBit1(&s, input);
	ZVector psk_tx=z8PSK_mod(input);
	int d;
	ZVector ofdm_tx=ofdm_encode(psk_tx, 16, d);
	ZVector n(ofdm_tx.taille);
	zRandGausStatePtr np;
	double v=convertsnrtosigma2(1.0/3.0/16,5);
	zRandGausInit(&np,0, v);
	zbRandGaus(&np, n);

	ZVector y=ofdm_decode(ofdm_tx+n,16,d);

	
	ZVector xest(y.taille);
	
	for(int k=0;k<256;k++) {
		double min=1e20;
		double diff;
		 DCplx s,x;
		for(int i=0;i<8;i++) {
			x=psk8.sym[i];
			diff = (y.vect[k]-x).mod2();
			if(diff<min) {
				min=diff;
				s=x;
			}
		}
		xest.vect[k]=s;
	}
	cout << psk_tx << xest << psk_tx-xest;
	return 1;
}

		