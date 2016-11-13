/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/07/21 10:59:01 $
 * $Revision: 1.1.2.11 $
 * $Id: interleaver.h,v 1.1.2.11 2006/07/21 10:59:01 syed Exp $
 ********************************************************************/
#ifndef _INTERLEAVER_H_
#define _INTERLEAVER_H_

#include "tools/utilitis.h"
#include "rand/randgen.h"
#include <iostream>
#include <fstream>
using namespace std;


#define INVALID_INTERLEAVED_VECTOR_SIZE         14100
#define INVALID_INTERLEAVER_STATE               14101
#define CANNOT_PAD_IN_INPLACE_MODE              14102
#define INUSABLE                                14103
#define BAD_PREDEFINED_MAPPING                  14104
#define CREATE_OPTIMIZED_INTERLEAVER_FAILED     14105

class CInterleaver{
	protected:
		WVector map;
		int padsize;
		dRandUniStatePtr interleaver_seed;
	public:
		CInterleaver(int size=0, int seed=0);
		CInterleaver(const WVector &pmap);
        CInterleaver(char type, int line, int col);
        CInterleaver(const CInterleaver &src);
        ~CInterleaver();
		ZVector Apply(const ZVector &src);
		DVector Apply(const DVector &src);
		WVector Apply(const WVector &src);
        void Apply(char *dest, int *destlen, char *src, int srclen );
		ZVector Extract(const ZVector &src);
		DVector Extract(const DVector &src);
		WVector Extract(const WVector &src);
        void Extract(char *dest, int *destlen, char *src, int srclen );

		void ipApply(const ZVector &src);
		void ipApply(const DVector &src);
		void ipApply(const WVector &src);
		void ipExtract(const ZVector &src);
		void ipExtract(const DVector &src);
		void ipExtract(const WVector &src);

		friend ostream &operator<<(ostream &os, CInterleaver &b);
		friend istream &operator>>(istream &is, CInterleaver &b);
        int Get(int position) { return map.vect[position]; }
        WVector GetMap() { return map; }
        int Size() { return map.taille; }
        bool operator==(const CInterleaver &c)  const; 
        bool operator!=(const CInterleaver &c)  const; 
		CInterleaver &operator=(const CInterleaver &c) ;
        int get_padsize() {return padsize; }
        void set_padsize(int pad) {padsize=pad; }



};

class CInterleaverException{
	public:
		int message;
		CInterleaverException(int mesg) : message(mesg){
            cerr << "CInterleaverException : Error " << mesg << " occured" << endl;
        }

};

ostream &operator<<(ostream &os, CInterleaver &b);
istream &operator>>(istream &is, CInterleaver &b) ;
int interleaver_selftest();
#endif

