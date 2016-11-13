#include <stdlib.h>
#include <time.h>
#include <tools/all.h>
#include <cdma/walsh.h>
#include <cdma/spread.h>
#include <modulation/psk.h>
#include <coding/block/ldpc.h>
#include <fading/fadingtdma.h>

int xmain() {
	WLDPCParityCheck ldpc(20,40);
	make_method method=Evencol;
	ldpc.make_dense_mixed(method);
	int i,nb_user_max=8,u;
	UVector<DVector> walsh(nb_user_max);
	for(u=0;u<nb_user_max;u++) {
		walsh.vect[u]=dWalshCode(nb_user_max, u);
	}
	TConstellation constel = psk8gray;

//	cout << walsh;
	TRandBitStatePtr seed;
	seed.iseed=7;
	WVector b(100*20*constel.bit_per_sym);
	bRandBit1(&seed, b);
#define use_ldpc
#ifdef use_ldpc
	DVector llr;
	ZVector tx = (*constel.modulate_function)(ldpc.encode(b));
#else
	cout << "no ldpc ... ! " << endl;
	ZVector tx = (*constel.modulate_function)(b);
#endif
	dRandUniStatePtr tseed;
	dRandUniInit(&tseed, 7, 0,1);
	ZMatrix canal = fadingtaps(10,1,1,tx.taille,512, &tseed);
	cout << "canal ... done " << endl;
	ZVector ctx = SpreadData(tx, walsh.vect[1]);
	for(i=0;i<canal.line;i++){
		for(int j=0;j<walsh.vect[1].taille;j++) {
			ctx.vect[i*walsh.vect[1].taille+j]*=canal.mat[i][0];
		}
	}
	zRandGausStatePtr nseed;
	ZVector ns(ctx.taille);
	cout << "canal ... done " << endl;
	for(int snr=0;snr<40;snr+=3) {
		zRandGausInit(&nseed, 0, convertsnrtosigma2(snr));
		zbRandGaus(&nseed, ns);
		ZVector rx;
		rx=ctx+ns;
		for(i=0;i<canal.line;i++){
			for(int j=0;j<walsh.vect[1].taille;j++) {
				rx.vect[i*walsh.vect[1].taille+j]*=canal.mat[i][0].conj();
			}
		}		
		cout << "canal ... done " << endl;
		ZVector urx = UnspreadData(rx, walsh.vect[1]);
#ifdef use_ldpc
		WVector bt = ldpc.decode((*constel.demodulate_function)(urx));
#else
		WVector bt = (*constel.demodulate_function)(urx);
#endif
		int tmp = bt.diff_count(b);
		cout << "snr = " << snr << " diff = " << tmp <<" / " << bt.taille<<endl;
		if(tmp==0) break;
	}/*
	ZVector c0 = tx.copy(0,40);
	cout << c0;
	ZVector c1 = SpreadData(c0, walsh.vect[1]);
	cout << c1;
	ZVector c2 =  UnspreadData(c1, walsh.vect[1]);
	cout << c2;
	ZVector zdiff =c2-c0;
	cout << zdiff;*/
	return 1;
}

int vmain() {
	int k=1, lm=20, ln=40;
	int lk=ln-lm;
	WVector b(k*lk*3*2);
	TRandBitStatePtr s;
	s.iseed=5454;
	DVector llr;
	bRandBit1(&s, b);
	WLDPCParityCheck ldpc(lm,ln);
	make_method method=Dense;
	ldpc.make_dense_mixed(method);
	WVector tx = ldpc.encode(b);
	srand(time(NULL));
	for(int t=0;t<7;t++) {
		int pos=rand()%tx.taille;
		cout << "pos = " << pos << endl;
		tx.vect[pos]=tx.vect[pos]?0:1;
	}
	
	WVector d=ldpc.decode(tx);
//	DVector llr;
	cout << "diff : " << b.diff_count(d)<<endl;
	return 1;
	TConstellation c;
	zRandGausStatePtr n;
	WVector p;
	ZVector noise4(k*lk*3*2/**2/2*/), noise8(k*lk*2*2/**3/3*/);
	ZVector nnoise4(k*lk*3/**2/2*/), nnoise8(k*lk*2/**3/3*/);
	for(int snr=0;snr<40;snr+=3) {
		int sum=0, di;
		zRandGausInit(&n, 0, convertsnrtosigma2(snr));
		zbRandGaus(&n, noise4);
		zbRandGaus(&n, noise8);
		zbRandGaus(&n, nnoise4);
		zbRandGaus(&n, nnoise8);

//-----------------------

		c = psk4;
		p=ldpc.decode((*c.demodulate_function)(noise4+(*c.modulate_function)(ldpc.encode(b))));
		di = p.diff_count(b);
		sum+=di;
		cout << "l psk4     > snr = " << snr << " : diff = " << di << "/" << p.taille<< "  (" << double(di)/p.taille << ")" << endl;
		
		p=(*c.demodulate_function)(nnoise4+(*c.modulate_function)(b));
		di = p.diff_count(b);
		sum+=di;
		cout << "n psk4     > snr = " << snr << " : diff = " << di << "/" << p.taille<< "  (" << double(di)/p.taille << ")" << endl;

//-----------------------

		c = psk4gray;
		p=ldpc.decode((*c.demodulate_function)(noise4+(*c.modulate_function)(ldpc.encode(b))));
		di = p.diff_count(b);
		sum+=di;
		cout << "l psk4gray > snr = " << snr << " : diff = " << di << "/" << p.taille<< "  (" << double(di)/p.taille << ")" <<  endl;
		p=(*c.demodulate_function)(nnoise4+(*c.modulate_function)(b));
		di = p.diff_count(b);
		sum+=di;
		cout << "n psk4gray > snr = " << snr << " : diff = " << di << "/" << p.taille<< "  (" << double(di)/p.taille << ")" << endl;

//-----------------------

		c = psk8;
		p=ldpc.decode((*c.demodulate_function)(noise8+(*c.modulate_function)(ldpc.encode(b))));
		di = p.diff_count(b);
		sum+=di;
		cout << "l psk8     > snr = " << snr << " : diff = " << di << "/" << p.taille<< "  (" << double(di)/p.taille << ")" <<  endl;
		
		p=(*c.demodulate_function)(nnoise8+(*c.modulate_function)(b));
		di = p.diff_count(b);
		sum+=di;
		cout << "n psk8     > snr = " << snr << " : diff = " << di << "/" << p.taille<< "  (" << double(di)/p.taille << ")" << endl;

//-----------------------

		c = psk8gray;
		p=ldpc.decode((*c.demodulate_function)(noise8+(*c.modulate_function)(ldpc.encode(b))));
		di = p.diff_count(b);
		sum+=di;
		cout << "l psk8gray > snr = " << snr << " : diff = " << di << "/" << p.taille<<  "  (" << double(di)/p.taille << ")" << endl;
		p=(*c.demodulate_function)(nnoise8+(*c.modulate_function)(b));
		di = p.diff_count(b);
		sum+=di;
		cout << "n psk8gray > snr = " << snr << " : diff = " << di << "/" << p.taille<< "  (" << double(di)/p.taille << ")" << endl;

//-----------------------

		if (sum==0) break;
	}

	return 1;

}



int main() {
	WLDPCParityCheck ldpc(10,20);
	ldpc.make_dense_mixed();
	ldpc.mackay();
	int tmp;
	TRandBitStatePtr s;
	s.iseed=77;
	WVector b(10);
	bRandBit1(&s,b);
	ofstream os("utilbsrc.txt");
	for(int i=0;i<b.taille;i++) {
		os << b.vect[i] << endl;
	}
	os.close();

	os.open("utilsbtest.enc");
	WVector e = ldpc.encode(b);
	for(i=0;i<e.taille;i++) {
		if(i!=0&&i%20==0) os << '\n';
		os << e.vect[i];
	}
	os << '\n';

	os.close();

	for(int t=0;t<3;t++) {
		int pos = rand()%e.taille;
		cout << "pos : " << pos << endl;
		e.vect[pos]=e.vect[pos]?0:1;
	}

	//cout << endl;

	os.open("utilsbtest.tx");
	for(i=0;i<e.taille;i++) {
		if(i!=0&&i%20==0) os << '\n';
		os << e.vect[i];
	}
	os << '\n';
	os.close();

//	e.RZtoNRZ();
	WVector tt= ldpc.decode(e);
	cout << " diff = " << tt.diff_count(b) << endl;
	cout << " prev "<< ldpc.prev_decode_status<<endl;
	cout << " e     " << e << endl;
	
	//cin >> tmp;
	//cout << ldpc;
	return 1;
}