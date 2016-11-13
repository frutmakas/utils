
#include "tools.h"

double log2(double x) {
	if (x<=0) return 0; 
	return log(x)/log(2);
}

int bin2oct(int *data, int data_left) {
	int x = data_left>=3 ?3 : data_left;
	int oct=0;
	for(int i=0;i<3;i++) {
		if (x>i) oct = (oct << 1) | data[i];
		else oct <<= 1; 
	}
	return oct;
}

void oct2bin(int pskvalue, int *data_out) {
	data_out[0]=(pskvalue & 4)==4 ? 1 : 0;
	data_out[1]=(pskvalue & 2)==2 ? 1:0;
	data_out[2]=(pskvalue & 1)==1 ? 1:0;
}
