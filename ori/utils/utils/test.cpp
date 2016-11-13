#include <stdio.h>
#include "../utilitis.h"

int main(int argc, char **argv){
// DCPLX TEST 
	DCplx a(3.0,-2.5), b(1.0,-1.0), c(2,1), d,e,f,g;
	double h=3;
	d=a+b;
	e=b+h;
	f=h;
	
	a.print();
	b.print();
	d.print();
	e.print();
	f.print();
	g=d-a;
	g.print();
	printf("\n");
	d=2-a;
	e=a-2;
	d.print();
	e.print();





	return 1;
}