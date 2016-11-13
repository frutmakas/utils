/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/07/21 10:57:08 $
 * $Revision: 1.1.2.28 $
 * $Id: ldpc.h,v 1.1.2.28 2006/07/21 10:57:08 syed Exp $
 ********************************************************************/
#ifndef _CODING_BLOCK_LPDC_H_
#define _CODING_BLOCK_LPDC_H_

#include "tools/utilitis.h"
#include "rand/randgen.h"
#include "tools/binsparse.h"
#include "queue/list.h"
#undef HEAPWORKAROUND


//#define MAX_ITERATION 120 
extern int MAX_ITERATION;
typedef enum { Evencol, Evenboth ,Dense, Mixed  } make_method;

/*ostream &operator<<(ostream &os, make_method &m);

istream &operator>>(istream &is, make_method &m);*/


class WLDPCParityCheck : public WBinarySparseMatrix {
protected:
	int nb_bit_per_col, no4cycle;
	BMatrix G;
	make_method G_type;
	dRandUniStatePtr seed; 
public:
	WVector prev_decode_status;
	WLDPCParityCheck(int line, int col, int nb_bit_per_col=3 ,make_method method=Evencol);
	~WLDPCParityCheck();
	WVector G_cols, G_rows;

	BMatrix make_dense_mixed(make_method method=Dense, int verbose=0); 
	void special();
	WLDPCParityCheck::WLDPCParityCheck();
	//void ReadFromFile(const char* filename);
	WVector encode(const WVector &bitsrc);
#ifdef LDPC_LLR_INCLUDED

	WVector decode(const DVector &input, DVector &llr);

#else

	WVector decode(const DVector &input);

#endif

	WVector operator*(const WVector &c);
	int GLine() {
		return G.line;
	}
	int GCol() {
		return G.col;
	}

	void mackay();

	friend ZVector operator*(const ZVector &z, const WLDPCParityCheck &l);
	friend DVector operator*(const DVector &z, const WLDPCParityCheck &l);
	friend WVector operator*(const WVector &z, const WLDPCParityCheck &l);
	friend ostream &operator<<(ostream &os, WLDPCParityCheck &b);
	friend istream &operator>>(istream &is, WLDPCParityCheck &b) ;

};

ZVector operator*(const ZVector &z, const WLDPCParityCheck &l);
DVector operator*(const DVector &z, const WLDPCParityCheck &l);
WVector operator*(const WVector &z, const WLDPCParityCheck &l);


#define ROW_INDEX_OUT_OF_BOUND			32
#define COL_INDEX_OUT_OF_BOUND			33
#define UNKNOWN_DENSE_GENERATOR_METHOD	34
#define INVALID_SOURCE					35
#define ENCODING_PARITY_CHECK_FAILED	36
#define DECODING_SUCCESSFUL				37
#define DECODING_FAILED					38



class CLDPCParityCheckException {
public : int err_mesg;
		 CLDPCParityCheckException(int erm) { err_mesg = erm; }
};

void make_ldpc (dRandUniStatePtr *seed, const WMatrix &parity_check, make_method method, int cb,	int no4cycle);
WMatrix gen_ldpc(int line, int col, int nb_bit_per_col ,make_method method);

ostream &operator<<(ostream &os, WLDPCParityCheck &b);
istream &operator>>(istream &is, WLDPCParityCheck &b) ;

#endif

