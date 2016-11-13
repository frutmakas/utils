#ifndef _SORT_THROW_EXCEPTION_H_
#define _SORT_THROW_EXCEPTION_H_

#define TYPE_QUICKSORT 1
#define TYPE_HEAPSORT 2

#include <string.h>

class CSortException {
	public:
		int type;
		char *error_message;
		CSortException(char *err_mesg, int type){
			error_message = new char[strlen(err_mesg)];
			strcpy(error_message, err_mesg);
			this->type = type;
		}
		~CSortException(){
			delete[] error_message;
		}
};

#endif

