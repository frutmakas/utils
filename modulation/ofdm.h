#ifndef _MODULATION_OFDM_H_
#define _MODULATION_OFDM_H_

#include "tools/utilitis.h"
#include "globaldef.h"
#include <iostream>
using namespace std;

class COFDMException {
  public:
	int message;
	COFDMException(int mesg) { 
		message = mesg; 
		cerr << "OFDM module error code : " << message << endl;
	}
};

ZVector ofdm_encode(const ZVector &data, int fft_size, int &data_padded);
ZVector ofdm_decode(const ZVector data, int fft_size, int data_padded);



#endif

