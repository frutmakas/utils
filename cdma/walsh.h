/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/02/08 23:27:37 $
 * $Revision: 1.1.2.4 $
 * $Id: walsh.h,v 1.1.2.4 2004/02/08 23:27:37 syed Exp $
 ********************************************************************/

#ifndef __CDMA_WALSH_H_
#define __CDMA_WALSH_H_
#include <tools/utilitis.h>

#define BAD_START_BIT 7778

class CCodeGeneratorException{
public : 
	int err_mesg;
	CCodeGeneratorException(int mesg) : err_mesg(mesg) {
		cerr << "CCodeGeneratorException : " << mesg << endl;
		printerror();

	}
	CCodeGeneratorException(const char * filename, int linenumber, int mesg) : err_mesg(mesg) {
		cerr << "CCodeGeneratorException : " << mesg << " at "<< filename << "@" << linenumber << endl;
		printerror();
	} 
protected:
	void printerror() {
		switch (err_mesg) {
			case INDEX_OUT_OF_RANGE	: 
				cerr << "CCodeGeneratorException : Index is not within permitted range" << endl;
				break;
			case BAD_START_BIT	: 
				cerr << "CCodeGeneratorException : Start bit must be 0 or 1" << endl;
				break;
			default :
				cerr << "CCodeGeneratorException : Unknown error code. Refer to definition" << endl;
				break;
		}
	}
};

DVector dWalshCode(int nb_user_max, int user_index, bool normalize=true, int startbit=1);

#endif

