#include <all.h>
#include <cdma/walsh.h>
#include <cdma/gold.h>
#include <cdma/spread.h>
#include<time.h>

int main(){
	int nb_user_max=8;

	UVector<DVector> walsh(nb_user_max);
	int i,j;
	for(i=0;i<nb_user_max;i++) {
		walsh.vect[i]=dWalshCode(nb_user_max, i);
	}
	//cout << walsh;
/*	for(i=0;i<nb_user_max;i++) cout << i << " : " << walsh.vect[i].sum() << endl;
	DMatrix W(nb_user_max, walsh.vect[0].taille);
	for(i=0;i<W.line;i++) {
		for(int j=0;j<W.col;j++) {
			W.mat[i][j]=walsh.vect[i].vect[j];
		}
	}
	cout << W;
	cout << W * W.transpose_nip()/nb_user_max;
*/	
	TRandBitStatePtr bseed;
	zRandGausStatePtr nseed;
	dRandUniStatePtr tseed;
	int timestamp = time(NULL);
	dRandUniInit(&tseed, timestamp,0,1);
	TConstellation constel = psk8gray;
	const int bitsize = 50000*(int)floor(log2(constel.size));
	UVector<WVector> userbits(nb_user_max);
	UVector<ZVector> cdmasym(nb_user_max);
	ZVector sym(nb_user_max)  ;
	UVector<ZMatrix> channel(nb_user_max);
	UVector<ZVector> txatt(nb_user_max);
	double nb_1 = 1.0/nb_user_max;
#define const_values
#ifdef const_values
	cout << " coding ... " <<endl;
	
	for(i=0;i<nb_user_max;i++) {
		userbits.vect[i]=WVector(bitsize);
		bRandBit1(&bseed, userbits.vect[i]);
		sym=(constel.modulate_function)(userbits.vect[i])*nb_1;
		channel.vect[i]=fadingtaps(10,1,1,sym.taille, 512, &tseed);
		for(j=0;j<sym.taille;j++) {
			sym.vect[j]*=channel.vect[i].mat[j][0];
		}
		cdmasym.vect[i]=SpreadData(sym, walsh.vect[i]);

	}
	cout << " decoding ... " <<endl;
#endif
	char filename[3000];
	sprintf(filename, "cdma-walsh-%d-%d", nb_user_max, timestamp);
	ofstream os(filename);

	for(int snr=0;snr<42;snr+=3) {
#ifndef const_values
		cout << " coding ... " <<endl;
		for(i=0;i<nb_user_max;i++) {
			userbits.vect[i]=WVector(bitsize);
			bRandBit1(&bseed, userbits.vect[i]);
			sym=(constel.modulate_function)(userbits.vect[i])*nb_1;
			channel.vect[i]=fadingtaps(10,1,1,sym.taille, 512, &tseed);
			for(j=0;j<sym.taille;j++) {
				sym.vect[j]*=channel.vect[i].mat[j][0];
			}
			cdmasym.vect[i]=SpreadData(sym, walsh.vect[i]);

		}
		cout << " decoding ... " <<endl;
#endif
		zRandGausInit(&nseed, 0, convertsnrtosigma2(snr));
		ZVector noise(cdmasym.vect[0].taille);

		zbRandGaus(&nseed, noise);
		//sum = cdmasym.vect[0];
		for(i=0;i<nb_user_max;i++) {
			noise+=cdmasym.vect[i];
		}
		// noised received
		// decoding 
		ZVector uns;
		WVector bitrx;
		for(i=0;i<nb_user_max;i++) {
			uns = UnspreadData(noise, walsh.vect[i]);
			for(j=0;j<uns.taille;j++) {
				uns.vect[j]*=channel.vect[i].mat[j][0].conj();
			}
			bitrx = (*constel.demodulate_function)(uns);
			cout << "snr=" << snr << " , user=" << i << ", diff="<< bitrx.diff_count(userbits.vect[i]) << " out of " << bitrx.taille<<endl;
			os << "snr=" << snr << " , user=" << i << ", diff="<< bitrx.diff_count(userbits.vect[i]) << " out of " << bitrx.taille<<endl;
		}
	}


	os.close();
	return 0;
}