// nspmatrice.cpp : Defines the entry point for the console application.
//
#include "../globaldef.h"

#ifdef _NSP_INTEL_

#include "nspdcplx.h"
#include "nspmatrice.h"
#include <stdio.h>

//#define DEBUG


DMatrix::DMatrix(int x, int y, double initial){
	if (x==0 || y==0) {
		col=line=0;mat=NULL;return;
	}
	col=x; line=y;
	mat=(double**)malloc(sizeof(double*)*col);
	size_t sline=sizeof(double)*line;
	for(int i=0;i<col;i++) {
		mat[i]=nspdMalloc(line);
		//mat[i]=(double*)malloc(sline);
		nspdbSet(initial, mat[i], line);
	}
}

DMatrix::~DMatrix(){
	if (mat==NULL) {
		col=line=0;
		return;
	}
		
	for(int i=0;i<col;i++)
		nspFree(mat[i]);
	free(mat);
	col=line=0;
	mat=NULL;
#ifdef DEBUG
	printf("\nDestruction DMatrix...\n");
#endif
}

DMatrix::DMatrix(const DMatrix &m){
	if(m.col==0 || m.line==0 || m.mat==NULL) {
		col=line=0;mat=NULL;
	}
	col=m.col; line=m.line;
	mat=(double**)malloc(sizeof(double*)*col);
	size_t sline=sizeof(double)*line;
	for(int i=0;i<col;i++) {
		mat[i]=nspdMalloc(line);
		//mat[i]=(double*)malloc(sline);
		nspdbCopy(m.mat[i], mat[i], line);
	}
}


void DMatrix::transpose(){
	if (col==0 || line==0) return;
	
	double **m = (double**)malloc(sizeof(double*)*line);
	for(int i=0;i<line;i++)
		m[i]=nspdMalloc(col);

	for(int x=0;x<col;x++) {
		for(int y = 0;y<line;y++)
			m[y][x]=mat[x][y];
	}
	
	for(x=0;x<col;x++) nspFree(mat[x]); 
	free(mat);

	x=col;col=line;line=x;
	mat=m;
}


DMatrix DMatrix::transpose_nip(){
	if (line==0 || col==0) return DMatrix();
	DMatrix m(line, col);

	for(int x=0;x<col;x++) {
		for(int y = 0;y<line;y++)
			m.mat[y][x]=mat[x][y];
	}

	return m;
}

ZMatrix DMatrix::ztranspose_nip(){
	if (line==0 || col==0) return ZMatrix();

	ZMatrix m(line, col);

	for(int x=0;x<col;x++) {
		for(int y = 0;y<line;y++){
			m.mat[y][x].re=mat[x][y];
			m.mat[y][x].im=0.0;
		}
	}
	return m;
}

DMatrix DMatrix::operator+(const DMatrix &rhs){
	if (line==0 || col==0) return DMatrix();
	if(rhs.col!=col || rhs.line!=line) return DMatrix();

	DMatrix m(col, line);
	for(int i=0;i<col;i++){
		nspdbAdd3(mat[i], rhs.mat[i], m.mat[i], line);
	}
	return m;
}

ZMatrix DMatrix::operator+(const ZMatrix &rhs){
	if (line==0 || col==0) return ZMatrix();
	if(rhs.col!=col || rhs.line!=line) return ZMatrix();

	ZMatrix m(col, line);
	double *im=nspdMalloc(line);
	nspdbZero(im, line);
	DCplx *tmp=nspzMalloc(line);
	for(int i=0;i<col;i++){
		nspzb2RealToCplx(mat[i], im, tmp, line);
		nspzbAdd3(tmp, rhs.mat[i], m.mat[i], line);
	}
	nspFree(tmp);
	nspFree(im);
	return m;
}

DMatrix DMatrix::operator-(const DMatrix &rhs){
	if (line==0 || col==0) return DMatrix();
	if(rhs.col!=col || rhs.line!=line) return DMatrix();

	DMatrix m(col, line);
	for(int i=0;i<col;i++){
		nspdbSub3(mat[i], rhs.mat[i], m.mat[i], line);
	}
	return m;
}

ZMatrix DMatrix::operator-(const ZMatrix &rhs){
	if (line==0 || col==0) return ZMatrix();
	if(rhs.col!=col || rhs.line!=line) return ZMatrix();

	ZMatrix m(col, line);
	double *im=nspdMalloc(line);
	nspdbZero(im, line);
	DCplx *tmp=nspzMalloc(line);
	for(int i=0;i<col;i++){
		nspzb2RealToCplx(mat[i], im, tmp, line);
		nspzbSub3(tmp, rhs.mat[i], m.mat[i], line);
	}
	nspFree(tmp);
	nspFree(im);
	return m;
}

DMatrix DMatrix::operator*(const DMatrix &rhs){
	if (line==0 || col==0) return DMatrix();
	if (this->col!=rhs.line) return DMatrix(); // verif size

	DMatrix t = this->transpose_nip();
	DMatrix n(rhs.col, this->line);

	for(int p = 0; p<this->line;p++) {
		for(int q=0;q<rhs.col;q++) {
			n.mat[q][p]=nspdDotProd(t.mat[p], rhs.mat[q], this->col);
		}
	}

	//t.Free();
	return n;
}


ZMatrix DMatrix::operator*(const ZMatrix &rhs){
	if (line==0 || col==0) return ZMatrix();
	if (this->col!=rhs.line) return ZMatrix();

	ZMatrix t = this->ztranspose_nip();
	ZMatrix n(rhs.col, this->line);

	for(int p = 0; p<this->line;p++) {
		for(int q=0;q<rhs.col;q++) {
			n.mat[q][p]=nspzDotProd(t.mat[p], rhs.mat[q], this->col);
		}
	}
	//t.Free();
	return n;
}


DMatrix DMatrix::operator*(const double &rhs){
	if (line==0 || col==0) return DMatrix();
	DMatrix n(*this);

	for(int i=0;i<col;i++)
		nspdbMpy1(rhs, n.mat[i], line);
	return n;
}

ZMatrix DMatrix::operator*(const DCplx &rhs){
	if (line==0 || col==0) return ZMatrix();
	ZMatrix n(*this);
	for(int i=0;i<col;i++)
		nspzbMpy1(rhs, n.mat[i], line);
	return n;

}

DMatrix DMatrix::operator/(const double &rhs){
	if (line==0 || col==0) return DMatrix();
	if (rhs==0) return DMatrix();
	DMatrix n(*this);
	double hs=1/rhs;
	for(int i=0;i<col;i++)
		nspdbMpy1(hs, n.mat[i], line);
	return n;
}

ZMatrix DMatrix::operator/(const DCplx &rhs){
	if (line==0 || col==0) return ZMatrix();
	if(rhs.re==0.0 && rhs.im==0) return ZMatrix();
	ZMatrix n(*this);
	DCplx one={1, 0};
	DCplx hs=nspzDiv(one, rhs);

	for(int i=0;i<col;i++)
		nspzbMpy1(hs, n.mat[i], line);
	return n;
}

bool DMatrix::operator==(const DMatrix &rhs){
	if(rhs.col!=col || rhs.line!=line) return false;
	for(int i=0;i<col;i++){
		for(int j=0;j<line;j++)
			if (rhs.mat[i][j]!=mat[i][j]) return false;
	}
	return true;
}

bool DMatrix::operator!=(const DMatrix &rhs){
	if(rhs.col!=col || rhs.line!=line) return true;
	for(int i=0;i<col;i++){
		for(int j=0;j<line;j++)
			if (rhs.mat[i][j]!=mat[i][j]) return true;
	}
	return false;
}

/*
DMatrix DMatrix::operator=(const DMatrix &rhs){
	return DMatrix(rhs);
}
*/

DMatrix &DMatrix::operator=(const DMatrix &rhs){
	if(rhs.line==0 || rhs.col==0) {

		if (this->mat!=NULL) {
			for(int k=0;k<this->col;k++) 
				nspFree(this->mat[k]);
			free(this->mat);
		}
			
		this->col=this->line=0;
		this->mat=0;
		return *this;
	}
	
	if (this->mat!=NULL) {
		for(int k=0;k<this->col;k++) 
			nspFree(this->mat[k]);
		free(this->mat);
	}

	this->col=rhs.col;
	this->line=rhs.line;

	this->mat = (double**)malloc(sizeof(double*)*col);

	for(int i=0;i<col;i++) {
		this->mat[i] = nspdMalloc(line);
		nspdbCopy(rhs.mat[i], this->mat[i], line);
	}

	return *this;
}

ZMatrix::ZMatrix(int col, int line){
	if (line==0 || col==0){
		this->col=this->line=0;mat=NULL;
		return;
	}
	this->col=col; this->line=line;
	mat=(DCplx**)malloc(sizeof(DCplx*)*col);
	for(int i=0;i<col;i++) {
		mat[i]=nspzMalloc(line);
		nspzbZero(mat[i], line);
	}
}

ZMatrix::ZMatrix (int col, int line, DCplx x){
	if (line==0 || col==0){
		this->col=this->line=0;mat=NULL;
		return;
	}
	
	this->col=col; this->line=line;
	mat=(DCplx**)malloc(sizeof(DCplx*)*col);
	for(int i=0;i<col;i++) {
		mat[i]=nspzMalloc(line);
		nspzbSet(x.re, x.im, mat[i], line);
	}
}

ZMatrix::ZMatrix(const ZMatrix& m){
	col=m.col; line=m.line;
	if (m.mat==NULL) { mat=NULL; return; }
	mat=(DCplx**)malloc(sizeof(DCplx*)*col);
	for(int i=0;i<col;i++) {
		mat[i]=nspzMalloc(line);
		nspzbCopy(m.mat[i], mat[i], line);
	}
}

ZMatrix::ZMatrix(const DMatrix& m){
	col=m.col; line=m.line;
	if (m.mat==NULL) { mat=NULL; return; }
	mat=(DCplx**)malloc(sizeof(DCplx*)*col);

	double *im=nspdMalloc(line);
	nspdbZero(im, line);
	DCplx *tmp=nspzMalloc(line);

	for(int i=0;i<col;i++) {
		mat[i]=nspzMalloc(line);
		nspzb2RealToCplx(m.mat[i], im, tmp, line);
		nspzbCopy(tmp, mat[i], line);
	}
	nspFree(im);
	nspFree(tmp);
}

ZMatrix::~ZMatrix(){
#ifdef DEBUG
	printf("\nDestruction of ZMatrix ...\n");
#endif
	if (mat==NULL) {
		col=line=0;
		return;
	}
	for(int i = 0; i < col; i++) {
		nspFree(mat[i]);
	}
	free(mat);
	col=line=0;
	mat=NULL;
}


/*void ZMatrix::Free(){
	for(int i=0;i<col;i++) nspFree(mat[i]);
	free(mat);
}
*/
void ZMatrix::transpose(){
	if (line==0 || col==0) return;
	DCplx **m = (DCplx**)malloc(sizeof(DCplx*)*line);
	for(int i=0;i<line;i++)
		m[i]=nspzMalloc(col);

	for(int x=0;x<col;x++) {
		for(int y = 0;y<line;y++)
			m[y][x]=mat[x][y];
	}
	
	for(x=0;x<col;x++) nspFree(mat[x]); 
	free(mat);

	x=col;col=line;line=x;
	mat=m;
}


ZMatrix ZMatrix::transpose_nip(){
	if (line==0 || col==0) return ZMatrix();
	ZMatrix m(line, col);

	for(int x=0;x<col;x++) {
		for(int y = 0;y<line;y++)
			m.mat[y][x]=mat[x][y];
	}

	return m;
}

ZMatrix  ZMatrix::operator+(const ZMatrix &rhs){
	if (line==0 || col==0) return ZMatrix();
	if(rhs.col!=col || rhs.line!=line) return ZMatrix();

	ZMatrix m(col, line);
	for(int i=0;i<col;i++){
		nspzbAdd3(mat[i], rhs.mat[i], m.mat[i], line);
	}
	return m;
}

ZMatrix  ZMatrix::operator+(const DMatrix &rhs){
	if (line==0 || col==0) return ZMatrix();
	if(rhs.col!=col || rhs.line!=line) return ZMatrix();

	double *im=nspdMalloc(line);
	nspdbZero(im, line);
	DCplx *tmp=nspzMalloc(line);

	ZMatrix m(col, line);

	for(int i=0;i<col;i++){
		nspzb2RealToCplx(rhs.mat[i], im, tmp, line);
		nspzbAdd3(mat[i], tmp, m.mat[i], line);
	}
	return m;

}

ZMatrix  ZMatrix::operator-(const ZMatrix &rhs){
	if (line==0 || col==0) return ZMatrix();
	if(rhs.col!=col || rhs.line!=line) return ZMatrix();

	ZMatrix m(col, line);
	for(int i=0;i<col;i++){
		nspzbSub3(mat[i], rhs.mat[i], m.mat[i], line);
	}
	return m;
}

ZMatrix  ZMatrix::operator-(const DMatrix &rhs){
	if (line==0 || col==0) return ZMatrix();
	if(rhs.col!=col || rhs.line!=line) return ZMatrix();

	double *im=nspdMalloc(line);
	nspdbZero(im, line);
	DCplx *tmp=nspzMalloc(line);

	ZMatrix m(col, line);

	for(int i=0;i<col;i++){
		nspzb2RealToCplx(rhs.mat[i], im, tmp, line);
		nspzbSub3(mat[i], tmp, m.mat[i], line);
	}
	return m;
}

ZMatrix ZMatrix::operator*(const ZMatrix &rhs){
	if (line==0 || col==0) return ZMatrix();
	if (this->col!=rhs.line) return ZMatrix(); // verif size

	ZMatrix t = this->transpose_nip();
	ZMatrix n(rhs.col, this->line);

	for(int p = 0; p<this->line;p++) {
		for(int q=0;q<rhs.col;q++) {
			n.mat[q][p]=nspzDotProd(t.mat[p], rhs.mat[q], this->col);
		}
	}

	//t.Free();
	return n;
}

ZMatrix ZMatrix::operator*(const DMatrix &rhs){
	if (line==0 || col==0) return ZMatrix();
	if(this->col!=rhs.line) return ZMatrix();

	double *im=nspdMalloc(line);
	nspdbZero(im, line);
	DCplx *tmp=nspzMalloc(line);

	ZMatrix t = this->transpose_nip();
	ZMatrix n(rhs.col, this->line);

	for(int p = 0; p<this->line;p++) {
		for(int q=0;q<rhs.col;q++) {
			nspzb2RealToCplx(rhs.mat[q], im, tmp, this->col);
			n.mat[q][p]=nspzDotProd(t.mat[p], tmp, this->col);
		}
	}

	//t.Free();
	nspFree(im);
	nspFree(tmp);
	return n;
}

ZMatrix ZMatrix::operator*(const double &rhs){
	if (line==0 || col==0) return ZMatrix();
	ZMatrix n(*this);
	DCplx tmp;
	tmp.re=rhs;
	tmp.im=0.0;

	for(int i=0;i<col;i++)
		nspzbMpy1(tmp, n.mat[i], line);
	return n;
}

ZMatrix ZMatrix::operator*(const DCplx &rhs){
	if (line==0 || col==0) return ZMatrix();
	ZMatrix n(*this);
	for(int i=0;i<col;i++)
		nspzbMpy1(rhs, n.mat[i], line);
	return n;
}

ZMatrix ZMatrix::operator/(const double &rhs){
	if (line==0 || col==0) return ZMatrix();
	if(rhs==0) return ZMatrix();

	ZMatrix n(*this);
	DCplx hs;
	hs.re = 1/rhs; hs.im=0.0;

	for(int i=0;i<col;i++)
		nspzbMpy1(hs, n.mat[i], line);
	return n;
}

ZMatrix ZMatrix::operator/(const DCplx &rhs){
	if (line==0 || col==0) return ZMatrix();
	if(rhs.re==0 && rhs.im==0) return ZMatrix();
	ZMatrix n(*this);
	DCplx one={1, 0};
	DCplx hs=nspzDiv(one, rhs);
	for(int i=0;i<col;i++)
		nspzbMpy1(hs, n.mat[i], line);
	return n;
}

ZMatrix &ZMatrix::operator=(const ZMatrix &rhs) {
	if (rhs.col==0 || rhs.line ==0) {
		if (this->mat!=NULL) {
			for(int k=0;k<this->col;k++) 
				nspFree(this->mat[k]);
			free(this->mat);
		}
		this->col = this->line = 0;
		return (*this);
	}
	this->col=rhs.col;
	this->line=rhs.line;
	this->mat=(DCplx**)malloc(sizeof(DCplx*)*col);
	for(int i=0;i<rhs.col;i++) {
		this->mat[i]=nspzMalloc(line);
		nspzbCopy(rhs.mat[i], this->mat[i], line);
	}
	return (*this);
}

ZMatrix &ZMatrix::operator=(const DMatrix &rhs) {
	if (rhs.col==0 || rhs.line ==0) {
		if (this->mat!=NULL) {
			for(int k=0;k<this->col;k++) 
				nspFree(this->mat[k]);
			free(this->mat);
		}
		this->col = this->line = 0;
		return (*this);
	}
	this->col=rhs.col;
	this->line=rhs.line;
	this->mat=(DCplx**)malloc(sizeof(DCplx*)*col);
	double *im=nspdMalloc(line);

	nspdbZero(im, line);
	DCplx *tmp=nspzMalloc(line);


	for(int i=0;i<rhs.col;i++) {
		this->mat[i]=nspzMalloc(line);
		nspzb2RealToCplx(rhs.mat[i], im, tmp, line);
		nspzbCopy(tmp, this->mat[i], line);
	}
	nspFree(im);
	nspFree(tmp);
	return (*this);
}

bool ZMatrix::operator==(const ZMatrix &rhs){
	if(rhs.col!=col || rhs.line!=line) return false;
	for(int i=0;i<col;i++){
		for(int j=0;j<line;j++)
			if (rhs.mat[i][j].re!=mat[i][j].re || rhs.mat[i][j].im!=mat[i][j].im) return false;
	}
	return true;
}

bool ZMatrix::operator!=(const ZMatrix &rhs){
	if(rhs.col!=col || rhs.line!=line) return true;
	for(int i=0;i<col;i++){
		for(int j=0;j<line;j++)
			if (rhs.mat[i][j].re!=mat[i][j].re || rhs.mat[i][j].im!=mat[i][j].im) return true;
	}
	return false;
}

DMatrix ZMatrix::real(){
	if (line==0 || col==0) return DMatrix();
	DMatrix n=DMatrix(col, line);
	for(int i=0; i<col;i++)
		nspzbReal(mat[i], n.mat[i], line);
	printf("\nreal.n : col = %d line=%d mat=%p", n.col, n.line, n.mat);
	return n;
}

DMatrix ZMatrix::imag(){
	if (line==0 || col==0) return DMatrix();
	DMatrix n=DMatrix(col, line);
	for(int i=0; i<col;i++)
		nspzbImag(mat[i], n.mat[i], line);
	return n;
}

DMatrix ZMatrix::abs(){
	if (line==0 || col==0) return DMatrix();
	DMatrix n=DMatrix(col, line);
	for(int i=0; i<col;i++)
		nspzbMag(mat[i], n.mat[i], line);
	return n;
}

DMatrix ZMatrix::arg(){
	if (line==0 || col==0) return DMatrix();
	DMatrix n=DMatrix(col, line);
	for(int i=0; i<col;i++)
		nspzbPhase(mat[i], n.mat[i], line);
	return n;
}		

ZMatrix ZMatrix::conj(){
	if (line==0 || col==0) return DMatrix();
	ZMatrix n=ZMatrix(col, line);
	for(int i=0; i<col;i++)
		nspzbConj2(mat[i], n.mat[i], line);
	return n;
}

ZMatrix operator*(const DCplx &lhs, const ZMatrix &rhs){
	if (rhs.line==0 || rhs.col==0) return ZMatrix();
	ZMatrix n(rhs);
	for(int i=0;i<rhs.col;i++)
		nspzbMpy1(lhs, n.mat[i], rhs.line);
	return n;
}

ZMatrix operator*(const double &lhs, const ZMatrix &rhs){
	if (rhs.line==0 || rhs.col==0) return ZMatrix();
	ZMatrix n(rhs);
	DCplx t;
	t.re=lhs;
	t.im=0.0;

	for(int i=0;i<rhs.col;i++)
		nspzbMpy1(t, n.mat[i], rhs.line);
	return n;
}

DMatrix operator*(const double &lhs, const DMatrix &rhs){
	if (rhs.line==0 || rhs.col==0) return DMatrix();
	DMatrix n(rhs);

	for(int i=0;i<rhs.col;i++)
		nspdbMpy1(lhs, n.mat[i], rhs.line);
	return n;
}

ZMatrix operator*(const DCplx &lhs, const DMatrix &rhs){
	if (rhs.line==0 || rhs.col==0) return ZMatrix();
	ZMatrix n(rhs);

	for(int i=0;i<rhs.col;i++)
		nspzbMpy1(lhs, n.mat[i], rhs.line);
	return n;
}

void print(const char* c, DMatrix &m){
	printf("\n matrix %s with addr 0x%p/%p", c, &m, m.mat);
	if(m.mat==NULL) { printf("\n ZERO MATRIX ...\n"); return;}
	printf("\n[ ");
	for(int j=0;j<m.line;j++){
		if(j!=0) printf("] ");
		printf("\n [ ");
		for(int i=0;i<m.col;i++) 
			printf("%4.4f ", m.mat[i][j]);
	}
	printf("] \n]\n");
}

void print(const char* c, ZMatrix &m){
	printf("\n matrix %s with addr 0x%p/%p", c, &m, m.mat);
	if(m.mat==NULL) { printf("\n ZERO MATRIX ...\n"); return;}
	printf("\n[ ");
	for(int j=0;j<m.line;j++){
		if(j!=0) printf("] ");
		printf("\n [ ");
		for(int i=0;i<m.col;i++) 
			printf("(%4.4f,%4.4f) ", m.mat[i][j].re, m.mat[i][j].im);
	}
	printf("] \n]\n");
}

#define NSPMATRICE_LIBRARY_TEST
#ifdef NSPMATRICE_LIBRARY_TEST

int main(int argc, char* argv[]) {
	//printf("Hello World!\n");
	printf("\n");
	DMatrix d1(2,4);
	for(int i=0;i<2;i++){
		for(int j=0;j<4;j++)
			d1.mat[i][j]=(i+1)*(j+1);
	}
	DCplx c1 = {0.5,0.2};
	ZMatrix z1=d1*c1;
	ZMatrix z2=z1+d1;
	ZMatrix z3=d1+z1;
	printf("\nZ3 : col = %d line=%d mat=%p", z3.col, z3.line, z3.mat);
	DMatrix d2=z3.real();
	printf("\nD2 : col = %d line=%d mat=%p", d2.col, d2.line, d2.mat);
	DMatrix d3=z3.imag();
	ZMatrix z5=d3/c1;
	ZMatrix z6=z5*c1;
	z6.transpose();
	ZMatrix z7=z5*z6;
	ZMatrix z8=z6*z5;

	DMatrix d4(3,2);
	for(i=0;i<3;i++){
		for(int j=0;j<2;j++)
			d4.mat[i][j]=(i+3)*(j-1);
	}
	
	DMatrix d5 = d4.transpose_nip();
	DMatrix d6=d5*d1;
	DMatrix d7=d1*d5;
	DMatrix d8=d1*d4;
	DMatrix d9=d4*d1;
	ZMatrix s1=d8/c1;
	DMatrix e1=d8/2.0;
	ZMatrix s2=s1/2.0;
	ZMatrix s3=e1*s2;
	ZMatrix s4=s2*e1;
	ZMatrix s5=d1*d4/2.0;
	ZMatrix s6;
	s6=s5;
	DMatrix d17;
	d17=d1;
	ZMatrix s7;
	s7=d17;



	print("d1",d1);
	print("z1",z1);
	print("z2",z2);
	print("z3",z3);
	print("d2",d2);
	print("d3",d3);
	print("z5",z5);
	print("z6",z6);
	print("z7",z7);
	print("z8",z8);
	print("d1",d1);
	print("d4",d4);
	print("d8",d8);
	print("s1",s1);
	print("e1",e1);
	print("s2",s2);
	print("s3",s3);
	print("s4",s4);
	print("s5",s5);
	print("s6",s6);
	print("d1",d1);
	print("d17",d17);
	print("s7",s7);

	return 0;
}


#endif

#endif //_NSP_INTEL_