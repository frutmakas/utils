/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/04/05 18:43:07 $
 * $Revision: 1.1.2.8 $
 * $Id: fadingtdma.h,v 1.1.2.8 2004/04/05 18:43:07 syed Exp $
 ********************************************************************/
 
#ifndef _FADING_H_
#define _FADING_H_

#include "globaldef.h"
#include <iostream>
using namespace std;

#ifndef _NSP_INTEL_

#include "tools/utilitis.h"
#include "rand/randgen.h"

#define FD_FLUCT_OUT_OF_RANGE 1280
#define TAPS_AMPLITUDE_SIZE_ERROR 1281
class CFadingTapsException{
	public:

		int ErrorCode;
		CFadingTapsException(const char * filename, int linenumber, int mesg) {
			ErrorCode = mesg;
			cerr << "CFadingTapsException: " << mesg << " at "<< filename << "@" << linenumber << endl;
			printerror();
		} 
		CFadingTapsException(int code):ErrorCode(code){
			cerr << "CFadingTapsException : " << code << endl; //<< " at "<< filename << "@" << linenumber << endl;
			printerror();
		}
protected:
	void printerror() {
		switch (ErrorCode) {
 			case FD_FLUCT_OUT_OF_RANGE : 
				cerr << "CFadingTapsException : The doppler fluctuation must be between 0 and 50" << endl;
				break;
			case TAPS_AMPLITUDE_SIZE_ERROR : 
				cerr << "CFadingTapsException : The channel taps amplitude vector size is not equal to filter length" << endl;
				break;
			default:
				cerr << "CFadingTapsException : Unknown error. Refer to error definition" << endl;
				break;
		}
	}

};

DVector freqd(long nb_rays, double max_doppler, dRandUniStatePtr *tapseed);
DMatrix fgwave(long nb_rays, long filter_length, double symbol_time, dRandUniStatePtr *tapseed);
cmplx_mat tapsgen(double  symbol_T, int filter_length, long nb_samples, long nb_rays, const DVector &fd, const DMatrix &g, dRandUniStatePtr *tapseed);
cmplx_mat fadingmain(double  symbol_time, long filter_length, long nb_samples, long nb_rays, double seed, const DVector &fd, const DMatrix &g, dRandUniStatePtr *tapseed);

cmplx_mat fadingtaps(double max_doppler, double symbol_time, long filter_length, long nb_samples, long nb_rays, dRandUniStatePtr *tapseed, double fd_fluct=20.0);
cmplx_mat fadingtaps(double max_doppler, double symbol_time, long filter_length, long nb_samples, long nb_rays, dRandUniStatePtr *tapseed, const DVector &tapsamplitude, double fd_fluct=20.0);



#else
#error "Header incompatible with Intel Signal Processing Performance Library"

#endif

#endif

