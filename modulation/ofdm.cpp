#include <iostream>
using namespace std;
#include "tools/utilitis.h"
#include "fft/fft.h"
#include "modulation/ofdm.h"
#include "tools/tools.h"
#include <math.h>


ZVector ofdm_encode(const ZVector &data, int fft_size, int &data_padded) {
	// test fft size 
	int cnt=0;
        int i=0,j=sizeof(int)*8;
	int shift=1;
	for(i=0;i<j && cnt<2;i++) {
		if(fft_size&shift)  cnt++;
		shift<<=1;
	}
	if(cnt!=1) throw COFDMException(INVALID_FFT_SIZE);
	// final data size 
	int in_block  = (int)ceil((double)data.taille/fft_size);
	data_padded = in_block*fft_size - data.taille;
	int out_size = in_block*fft_size;
	ZVector output(out_size);
	for(i=0 ;i<data.taille;i+=fft_size) {
		if((i+fft_size)>data.taille) break;
		output.insert(ifft(data.copy(i,i+fft_size)),i);
	}
	if(i<data.taille) {
		ZVector tmp(fft_size, DCplx());
		tmp.insert(data.copy(i, data.taille),0);
		output.insert(ifft(tmp),i);
	}
	return output;
}

ZVector ofdm_decode(const ZVector data, int fft_size, int data_padded) {
	int cnt=0, i=0,j=sizeof(int)*8;
	int shift=1;
	for(i=0;i<j && cnt<2;i++) {
		if(fft_size&shift) cnt++;
		shift<<=1;
	}
	if(cnt!=1) throw COFDMException(INVALID_FFT_SIZE);
	if(data.taille%fft_size!=0) throw COFDMException(INVALID_OFDM_DATA_SIZE_TO_DECODE);
	int in_block = data.taille/fft_size;
	int out_size = (in_block*fft_size);
	ZVector output_tmp(out_size);
	for(i=0;i<data.taille;i+=fft_size) {
		output_tmp.insert(fft(data.copy(i,i+fft_size)),i);
	}
	return output_tmp.copy(0, output_tmp.taille-data_padded);
}



