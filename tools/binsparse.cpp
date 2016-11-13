/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2006/07/28 11:03:32 $
 * $Revision: 1.1.2.15 $
 * $Id: binsparse.cpp,v 1.1.2.15 2006/07/28 11:03:32 syed Exp $
 ********************************************************************/
#include "tools/binsparse.h"
#include "tools/utilitis.h"
#include "queue/list.h"
#include "globaldef.h"
//#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;


//#pragma comment( user, "Source File : " __FILE__ ". Compiled on " __TIMESTAMP__ ) 

WBinarySparseMatrix::WBinarySparseMatrix(int line, int col){ //:rowsum(line,0), colsum(col, 0) {
	this->col=col; this->line=line;
	col_index = new TSortedList<CSparseMatrixColInformation>[col];
//	colsum = new int[col];
	row_index = new TSortedList<int>[line];
//	rowsum = new int[line];
}

WBinarySparseMatrix::~WBinarySparseMatrix(){
	delete[] row_index;
	delete[] col_index;
//	delete[] rowsum;
//	delete[] colsum;
}

void WBinarySparseMatrix::Set(int row, int col){
	if (col_index[col].isExist(row)!=-1) return;
	col_index[col].Insert(CSparseMatrixColInformation(row));
	//colsum.vect[col]++;
	row_index[row].Insert(col);
	//rowsum.vect[row]++;
}

void WBinarySparseMatrix::Unset(int row, int col){
	if (col_index[col].DeleteByValue(row)!=0){
//		colsum.vect[col]--;
		row_index[row].DeleteByValue(col);
//		rowsum.vect[row]--;
	}
}

int WBinarySparseMatrix::operator()(int row, int col){
	if (row>=this->line|| col >= this->col) throw WBinarySparseMatrixException(INDEX_OUT_OF_RANGE);
	if(row<0 || col <0) throw WBinarySparseMatrixException(INVALID_RANGE);
	return col_index[col].isExist(row)==-1 ? 0 : 1;
}

BMatrix WBinarySparseMatrix::get_matrix(){
	BMatrix b(line, col);
	for(register int i=0;i<col;i++){
		for(register int j=1;j<=col_index[i].count;j++){
			b.Set(col_index[i].Get(j).position,i);
		}
	}
	return b;
}

BIMatrix WBinarySparseMatrix::get_clmatrix(){
	BIMatrix b(line, col);
	for(register int i=0;i<col;i++){
		for(register int j=1;j<=col_index[i].count;j++){
			b.Set(col_index[i].Get(j).position,i);
		}
	}
	return b;
}

WBinarySparseMatrix::WBinarySparseMatrix(){ //:rowsum(0), colsum(0){
	row_index=NULL;col_index=NULL;
}

WBinarySparseMatrix &WBinarySparseMatrix::operator=(const WBinarySparseMatrix& r){
	if (col_index!=NULL) {
		delete[] col_index;
		delete[] row_index;
//		rowsum.~WVector();
//		colsum.~WVector();
	}

	this->col=r.col;
	this->line=r.line;
	col_index = new TSortedList<CSparseMatrixColInformation>[col];
//	colsum = WVector(col,0);
	row_index = new TSortedList<int>[line];
	//colsum = WVector(line,0);

	for(register int i=0;i<col;i++) {
		col_index[i]=r.col_index[i];
	}
	for(register int j=0;j<line;j++) {
		row_index[j]=r.row_index[j];
	}
//	this->colsum = r.colsum;
//	this->rowsum = r.rowsum;
	return *this;
}

// this function is not working !!!!
int WBinarySparseMatrix::sparse_decompose(WBinarySparseMatrix **L, WBinarySparseMatrix **U){
	*L = new WBinarySparseMatrix(this->line, this->line);
	*U = new WBinarySparseMatrix(this->line, this->col);
	int M=this->line, N=this->col;
	int K=M;

	int pr=col*line, e_col=0, e_row=0;
	int *rinv=new int[M];
	int *cinv=new int[N];
	WVector rows(M), cols(N);
	int abandon_when=0, abandon_number=0;
	int cc, cc3;

	// int *rcnt=new int[M]; // strategy = minprod ==> B.rowsum

	int *acnt=new int[M+1]; // abandon_number>0
	WBinarySparseMatrix B;
	B=*this;
	register int i,j;//, c;
	for(i=0;i<M;i++) rows.vect[i]=rinv[i]=i;
	for(j=0;j<N;j++) cols.vect[j]=cinv[j]=j;
	int found, x, cc2, y, r=0, k=0, f, fn;

	int nnf=0;
	for(i=0;i<K;i++){
		found=0;
		for(j=i;j<N;j++){
			x=cols.vect[j];
			cc2=B.col_index[x].count;//colsum.vect[x];
			y=1;
			try {
				for(y=1;y<B.col_index[x].count/*colsum.vect[x]*/;y++) {
//				while((r=B.col_index[x].Get(y))!=-1){
					r=B.col_index[x].Get(y).position;
					if(rinv[r]>=i){
						int cr2 = B.row_index[r].count;//rowsum.vect[r];
						if(!found || cc2==1 || (cc2-1)*(cr2-1) <pr) {
							found=1;
							pr=cc2==1 ?0 : (cc2-1)*(cr2-1);
							e_col=x; //e=e2;
							e_row=r;
							k=j;
						}
					}
//					y++; // e2=nxtincol(e2);
				}
			} catch (...) {
				// on a parcouru toutes les lignes ... 
			}
		}

		if(!found){
			nnf++;
		}

		if(found) {
			if (cinv[e_col]!=k) throw WBinarySparseMatrixException(DECOMPOSE_ERROR);
			cols.vect[k]=cols.vect[i];
			cols.vect[i]=e_col;

			cinv[cols.vect[k]]=k;
			cinv[cols.vect[i]]=i;

			k=rinv[e_row];

			if(k<i) throw WBinarySparseMatrixException(DECOMPOSE_ERROR);
			rows.vect[k]=rows.vect[i];
			rows.vect[i]=r;
			rinv[rows.vect[k]]=k;
			rinv[rows.vect[i]]=i;
		}

		// update L U B
		y=1;
		try {
			while ((f=B.col_index[cols.vect[i]].Get(y).position)>=0) {
				fn = B.col_index[cols.vect[i]].Get(++y).position;
				k=f;
				if(rinv[k]>i) {
					try { // sparse_add_row start
						int n=1, m, a;
						while((m=B.row_index[e_row].Get(n++))>=0) {
							if((a=B.row_index[k].isExist(m))==-1) {
								B.Set(k,m);
							} else {
								B.Unset(k, B.row_index[k].Get(a));
							}
						}
					}catch(...){
					}// sparse_add_row stop
					(*L)->Set(k, i);
				} else if (rinv[k]<i) {
					(*U)->Set(rinv[k], cols.vect[i]);
				} else {
					(*L)->Set(k, i);
					(*U)->Set(i, cols.vect[i]);
				}
				f=fn;
			}
		} catch(...) {

		}
		
		// get rid of all entries in the current column B just to save space ... 
		while(B.col_index[cols.vect[i]].count/* colsum.vect[cols.vect[i]]*/>1) {
			B.Unset(B.col_index[cols.vect[i]].Get(1).position, cols.vect[i]);
		}

		if(abandon_number>0 &&  i==abandon_when) {
			for(k=0;k<M+1;k++){
				acnt[k]=0;
			}
			for(j=0;j<N;j++){
				acnt[B.col_index[j].count/*colsum.vect[j]*/]++;
			}
			cc=abandon_number;
			k=M;
			while(acnt[k]<cc){
				cc-=acnt[k--];
				if(k<0) throw WBinarySparseMatrixException(DECOMPOSE_ERROR);
			}
			cc2=0;
			for(j=0;j<N;j++){
				cc3 = B.col_index[j].count;//colsum.vect[j];
				if(cc3>k || cc3==k && cc> 0){
					if(cc3==k) cc--;
					while(B./*colsum.vect*/col_index[j].count>1){
						B.Unset(B.col_index[j].Get(1).position, j);
					}
					cc2++;
				}
			}
			if(cc2!=abandon_number) throw WBinarySparseMatrixException(DECOMPOSE_ERROR);
		}
		// 
	}


	// get rid of all entries
	for(i=K;i<M;i++){
		while((*L)->row_index/*rowsum.vect*/[rows.vect[i]].count>1){
			(*L)->Unset(rows.vect[i], (*L)->row_index[rows.vect[i]].Get(1));
		}
	}

	delete[] rinv;
	delete[] cinv;
	delete[] acnt;
	return nnf;
}

void WBinarySparseMatrix::ReadFromFile(const char* filename){
	FILE *f=fopen(filename,"rt");
	if(!f) return;
	char  buffer[3000];
	fgets(buffer,3000, f);
	char *token=strtok(buffer, ",");
	int line = atoi(token);
	//fgets(buffer,3000, f);
	token=strtok(NULL, ",");
	int col= atoi(token);
	//construct
	this->col=col; this->line=line;
	col_index = new TSortedList<CSparseMatrixColInformation>[col];
//	colsum = new int[col];
	row_index = new TSortedList<int>[line];
	//
	for(int i =0;i<line;i++){
		fgets(buffer,3000,f);
		token=strtok(buffer,":");
		int ll=atoi(token);
		//int ll=atoi(strtok(token,","));
		//token=strtok(buffer,":")
		while(token!=NULL) {
			token=strtok(NULL, " ");
			if(token) Set(ll, atoi(token));
		}
	}
	fclose(f);
}

ostream &operator<<(ostream &os, WBinarySparseMatrix &b){
	if (os.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	if (b.line==0 || b.col==0) {
		os << "0 : []" << endl;
		return os; 
	}
	os << b.line << "," << b.col << endl;
	for(register int m=0;m<b.line;m++){
		os << m << "," << b.row_index[m].count << " : " ;
		for(register int n=1;n<=b.row_index[m].count;n++){
			os << b.row_index[m].Get(n) << " " ;
		}
		os << endl;
	}
	return os;
}

ifstream &operator>>(ifstream &is, WBinarySparseMatrix &b) {
	if (is.bad()) throw CUtilitisException(__FILE__, __LINE__, BAD_FILE);
	char junk;
	int line, col, elt, ijunk, count;
	is >> line >> junk >> col;
	b.~WBinarySparseMatrix();
	b=WBinarySparseMatrix(line, col);
	for(register int m=0;m<line;m++) {
		is >> ijunk >> junk >> count >> junk;
		for(register int n=0;n<count;n++) {
			is >> elt;
			b.Set(m,elt);
		}
	}
	return is;
}
