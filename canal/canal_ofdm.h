/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2003/06/02 09:07:01 $
 * $Revision: 1.1.2.5 $
 * $Id: canal_ofdm.h,v 1.1.2.5 2003/06/02 09:07:01 syed Exp $
 ********************************************************************/

#ifndef _CANAL_OFDM_H_
#define _CANAL_OFDM_H_

ZMatrix ofdm_time_frequency_channel_estimation(const ZMatrix &C, const ZMatrix &D, 
											   int training_length, int Lf,  
											   double sigma, double fd /*max_doppler*/, 
											   ZMatrix &nC) ;
DMatrix  MIMO1OFDM(DMatrix &TEB ) ;

#endif
