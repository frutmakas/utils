/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/08/04 12:31:17 $
 * $Revision: 1.1.2.7 $
 * $Id: binmat.h,v 1.1.2.7 2004/08/04 12:31:17 syed Exp $
 ********************************************************************/
#ifndef _TOOLS_BINMMAT_H_
#define _TOOLS_BINMMAT_H_
/*#include <iostream.h>
#include <fstream.h>*/
#include <iostream>
using namespace std;

class BVector {
public:
	void Force(int position, bool value);
	int real_size, word_count;
	int *vect;
	BVector(int nb_element);
	BVector(const BVector &bv);
	BVector();
	~BVector();
	BVector& operator=(const BVector &bv);
	int operator[](int position) const;
	bool operator==(const BVector &bv) const;
	BVector &operator^=(const BVector &bv);
	bool operator!=(const BVector &bv) const;
	BVector operator+(const BVector &bv) const;
	BVector operator-(const BVector &bv) const;
	BVector operator^(int boolean) const;
	BVector operator^(const BVector &bv) const;
	BVector operator~() const;
	BVector operator!() const;
	BVector operator&(const BVector &bv) const;
	BVector operator|(const BVector &bv) const;
	void Set(int position) const;
	void Reset(int position) const;
	void Toggle(int position) const;
	void Clear() const;
	BVector operator>>(int shift) const;
	BVector operator<<(int shift) const;
};

ostream& operator<<(ostream &os, const BVector& bv);
istream& operator>>(istream &is, BVector &bv);

class BVectorException{
public:
	int mesg;
	BVectorException(int msg):mesg(msg){}
};
class BIMatrix;
class BMatrix{
	public:
		int line, col;
		BVector *mat;
		BMatrix(int line, int col);
		BMatrix();
		BMatrix(const BMatrix &bv);
		~BMatrix();
		bool operator==(const BMatrix &bv) const;
		bool operator!=(const BMatrix &bv) const;
		BMatrix& operator=(const BMatrix &bv);
		BMatrix& operator=(const BIMatrix &bv);
		BMatrix operator+(const BMatrix &bv) const;
		BMatrix operator^(const BMatrix &bv) const;
		BMatrix operator~() const;
		BMatrix operator&(const BMatrix &bv) const;
		BMatrix operator|(const BMatrix &bv) const;
		BMatrix operator*(const BMatrix &bv) const;
		void Set(int line, int col) ;
		void Reset(int line, int col) ;
		void Toggle(int line, int col);
};

class BIMatrix{
	public:
		int line, col;
		BVector *mat;
		BIMatrix(int line, int col);
		BIMatrix();
		BIMatrix(const BIMatrix &bv);
		~BIMatrix();
		bool operator==(const BIMatrix &bv) const;
		bool operator!=(const BIMatrix &bv) const;
		BIMatrix& operator=(const BIMatrix &bv);
		BIMatrix operator+(const BIMatrix &bv) const;
		BIMatrix operator^(const BIMatrix &bv) const;
		BIMatrix operator~() const;
		BIMatrix operator&(const BIMatrix &bv) const;
		BIMatrix operator|(const BIMatrix &bv) const;
		BIMatrix operator*(const BIMatrix &bv) const;
		void Set(int line, int col) ;
		void Reset(int line, int col) ;
		void Toggle(int line, int col);
};


ostream& operator<<(ostream &os, const BMatrix& bv);
istream& operator>>(istream &is, BMatrix &bv);
ostream& operator<<(ostream &os, const BIMatrix& bv);
istream& operator>>(istream &is, BIMatrix &bv);

class BMatrixException{
public:
	int mesg;
	BMatrixException(int msg):mesg(msg){}
};

typedef BMatrixException BIMatrixException;


#endif

