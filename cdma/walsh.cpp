/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/06/06 09:20:29 $
 * $Revision: 1.1.2.8 $
 * $Id: walsh.cpp,v 1.1.2.8 2006/06/06 09:20:29 syed Exp $
 ********************************************************************/

#include <cdma/walsh.h>
#include <tools/utilitis.h>
#include <tools/tools.h>


DVector dWalshCode(int nb_user_max, int user_index, bool normalize ,int startbit) {
	if (user_index>=nb_user_max) throw CCodeGeneratorException(INDEX_OUT_OF_RANGE);
	if (startbit>1 || startbit < 0) throw CCodeGeneratorException(BAD_START_BIT);
	int reclevel = (int)ceil(log_2(nb_user_max));
	int i, taille=1<<(reclevel);
	WVector code(taille), tmp(code);
	DVector walsh(code);
	code.vect[taille-1]=startbit==1?1:-1;
	//int user_taille =(int)ceil(log_2(user_index));
	for(i=0;i<reclevel;i++) {
		int invbit = (user_index & (1<<i))!=0;
		tmp = code << (1<<i);
		for(int j=taille-(1<<i); j<taille; j++) {
			if(invbit) {
				tmp.vect[j]=code.vect[j]==-1 ? 1 :-1;
			} else {
				tmp.vect[j]=code.vect[j]==-1 ? -1 :1;
			}
		}
		code = tmp;
	}
	double coef=1.0;
	if(normalize) {
		coef = 1.0/sqrt((double)taille);
		return coef*code;
	}
	return code;
/*	for(i=0;i<taille;i++) 
		walsh.vect[i]=coef*(code.vect[i]==0?-1.0:1.0);
	return walsh;	*/
}
