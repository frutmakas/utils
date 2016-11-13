/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: jebon $
 * $Date: 2003/03/31 14:25:56 $
 * $Revision: 1.1.2.12 $
 * $Id: binsparse.h,v 1.1.2.12 2003/03/31 14:25:56 jebon Exp $
 ********************************************************************/
#ifndef _TOOLS_BINSPARSE_h
#define _TOOLS_BINSPARSE_h

#include "queue/list.h"
#include "tools/utilitis.h"
#include "tools/binmat.h"

class CSparseMatrixColInformation{
public : 
	int position; 
	double Zmn, Lmn;
	bool operator<(const CSparseMatrixColInformation  &c) {
		return (position<c.position);
	}
	bool operator<=(const CSparseMatrixColInformation  &c) {
		return (position<=c.position);
	}
	bool operator>(const CSparseMatrixColInformation  &c) {
		return (position>c.position);
	}
	bool operator>=(const CSparseMatrixColInformation  &c) {
		return (position>=c.position);
	}
	bool operator==(const CSparseMatrixColInformation  &c) {
		return (position==c.position);
	}
	bool operator<(int c) {
		return (position<c);
	}
	bool operator<=(int c) {
		return (position<=c);
	}
	bool operator>(int c) {
		return (position>c);
	}
	bool operator>=(int c) {
		return (position>=c);
	}
	bool operator==(int c) {
		return (position==c);
	}

	CSparseMatrixColInformation &operator=(const CSparseMatrixColInformation  &c) {
		this->position=c.position;
		this->Zmn = c.Zmn;
		this->Lmn = c.Lmn;
		return *this;
	}


	CSparseMatrixColInformation(){
		position=0;
		Zmn=Lmn=0.0;
	}
	CSparseMatrixColInformation(const CSparseMatrixColInformation &c){
		position=c.position;
		Lmn=c.Lmn;
		Zmn=c.Zmn;
	}
	CSparseMatrixColInformation(int position, double zmn=0.0, double lmn=0.0) {
		this->position=position;
		Lmn=lmn;
		Zmn=zmn;
	}
};

class WBinarySparseMatrix{
protected : 
	TSortedList<CSparseMatrixColInformation> *col_index;
	TSortedList<int> *row_index;
		int line, col;
		//WVector rowsum, colsum;
	public:
		WBinarySparseMatrix(int line, int col);

		int Line() { return line; }

		int Col() { return col; }

		WBinarySparseMatrix();

		int sparse_decompose(WBinarySparseMatrix **L, WBinarySparseMatrix **U);
		void ReadFromFile(const char* filename);
		~WBinarySparseMatrix();
		void Set(int row, int col);
		void Unset(int row, int col);
		int operator()(int row, int col);
		BMatrix get_matrix();
		BIMatrix get_clmatrix();

		WBinarySparseMatrix &operator=(const WBinarySparseMatrix &r);
		friend ostream &operator<<(ostream &os, WBinarySparseMatrix &b);
};

class WBinarySparseMatrixException{
	public:
		int error_code;
		WBinarySparseMatrixException(int err_code){ error_code=err_code; }
};


ostream &operator<<(ostream &os, WBinarySparseMatrix &b);
ifstream &operator>>(ifstream &is, WBinarySparseMatrix &b) ;
#endif

