#include <iostream.h>
#include "../utilitis.h"
#include <conio.h>

int dcplxmain(int argc, char **argv){
	DCplx z1,z2(2,3);
	z1=3;
	DCplx z3;
	z3=z1+z2;
	cout << z1 <<'+'<<z2<<"="<<z3 <<endl;
	cout << z1 <<'+'<<z2<<"="<< z1+z2 <<endl;
	cout << z1 <<'-'<<z2<<"="<< z1-z2 <<endl;
	cout << z1 <<'*'<<z2<<"="<< z1*z2 <<endl;
	cout << z1 <<'/'<<z2<<"="<< z1/z2 <<endl;
	cout << z3 <<'/'<<z2 <<"="<< z3/z2 <<endl;
	cout << endl;

	cout << z3 <<'+' << -3 <<"="<< z3+(-3) <<endl;
	cout << z3 <<'-' << -3 <<"="<< z3-(-3) <<endl;
	cout << z3 <<'*' << -3 <<"="<< z3*(-3) <<endl;
	cout << z3 <<'/' << -3 <<"="<< z3/(-3) <<endl;
	cout << endl;

	cout << -5 <<'+' << z3 <<"="<< (-5)+z3 <<endl;
	cout << -5 <<'-' << z3 <<"="<< (-5)-z3 <<endl;
	cout << -5 <<'*' << z3 <<"="<< (-5)*z3 <<endl;
	cout << -5 <<'/' << z3 <<"="<< (-5)/z3 <<endl;
	cout << endl;

	cout << z3 << "==" << z1 << "=" ;
	if (z3==z1)  cout << "true" << endl; else cout << "false" << endl;
	cout << z3 << "==" << z3 << "="  ;
	if (z3==z3)  cout << "true" << endl; else cout << "false" << endl;
	cout << z3 << "==" << 3.0 << "=" ;
	if (z3==3.0)  cout << "true" << endl; else cout << "false" << endl;
	cout << z2 << "==" << 2.0 << "=" ;
	if (z2==2.0)  cout << "true" << endl; else cout << "false" << endl;
	cout << z1 << "==" << 3.0 << "=" ;
	if (z1==3.0)  cout << "true" << endl; else cout << "false" << endl;
	cout << endl;

	cout << z3 << "!=" << z1 << "=" ;
	if (z3!=z1)  cout << "true" << endl; else cout << "false" << endl;
	cout << z3 << "!=" << z3 << "="  ;
	if (z3!=z3)  cout << "true" << endl; else cout << "false" << endl;
	cout << z3 << "!=" << 3.0 << "=" ;
	if (z3!=3.0)  cout << "true" << endl; else cout << "false" << endl;
	cout << z2 << "!=" << 2.0 << "=" ;
	if (z2!=2.0)  cout << "true" << endl; else cout << "false" << endl;
	cout << z1 << "!=" << 3.0 << "=" ;
	if (z1!=3.0)  cout << "true" << endl; else cout << "false" << endl;

	cout << endl;

	cout << z3 << "+= " << 3 << "=";
	z3+=3;
	cout << z3 << endl;

	cout << z3 << "+= " << z2 << "=";
	z3+=z2;
	cout << z3 << endl;

	cout << z3 << "-= " << 3 << "=";
	z3-=3;
	cout << z3 << endl;

	cout << z3 << "-= " << z2 << "=";
	z3-=z2;
	cout << z3 << endl;

	cout << z3 << "*= " << 3 << "=";
	z3*=3.0;
	cout << z3 << endl;

	cout << z3 << "*= " << z2 << "=";
	z3*=z2;
	cout << z3 << endl;


	cout << z3 << ".mag() = " << z3.mag() << endl;
	cout << z3 << ".mod2() = " << z3.mod2() << endl;
	cout << z3 << ".phase() = " << z3.phase() << endl;



	return 1;

}

int main(int argc, char **argv){
	WVector dv1(20), dv2;
	for(register int i=0; i<dv1.taille;i++) dv1.vect[i]=i;
	cout << dv1;
	dv2=dv1+3;
/*	cout << dv1 << "+" << 3 << dv2 << endl;
	cout << dv1 << "+" << 3 << dv1+3 << endl;
	cout << dv1 << "-" << 3 << dv1-3 << endl;
	cout << dv1 << "*" << 3 << dv1*3 << endl;
	cout << dv1 << "/" << 3 << dv1/3 << endl;
	cout << endl;
	cout << 3 << "+" << dv1 << 3+dv1 << endl;
	cout << 3 << "-" << dv1 << 3-dv1 << endl;
	cout << 3 << "*" << dv1 << 3*dv1 << endl;
*/
/*	cout << dv1 << "+" << c << z1 << endl;
	cout << dv1 << "+" << c << dv1+c << endl;
	cout << dv1 << "-" << c << dv1-c << endl;
	cout << dv1 << "*" << c << dv1*c << endl;
	cout << dv1 << "/" << c << dv1/c << endl;
	cout << endl;
	cout << c << "+" << dv1 << c+dv1 << endl;
	cout << c << "-" << dv1 << c-dv1 << endl;
	cout << c << "*" << dv1 << c*dv1 << endl;
	cout << endl;
*/
	DVector n(dv2);
/*	cout << dv1 << "+" << n << "=" << dv1+n << endl;
	cout << dv1 << "-" << n << "=" << dv1-n << endl;
	cout << dv1 << "*" << n << "=" << dv1*n << endl;
	cout << endl;

	cout << n << "+" << dv1 << "=" << n+dv1 << endl;
	cout << n << "-" << dv1 << "=" << n-dv1 << endl;
	cout << n << "*" << dv1 << "=" << n*dv1 << endl;
	cout << endl;

*/	DVector v1(dv1);
/*	cout << v1 << "+" << n << "=" << v1+n << endl;
	cout << v1 << "-" << n << "=" << v1-n << endl;
	cout << v1 << "*" << n << "=" << v1*n << endl;
	cout << endl;
*/
	DCplx c(-1,1);
	ZVector z1=dv1+c;

/*	cout << z1 << "+" << dv2 << "=" << z1+dv2 << endl;
	getch();
	cout << z1 << "-" << dv2 << "=" << z1-dv2 << endl;
	getch();
	cout << z1 << "*" << dv2 << "=" << z1*dv2 << endl;
	cout << endl;
	getch();

	cout << dv2 << "+" << z1 << "=" << dv2+z1 << endl;
	getch();
	cout << dv2 << "-" << z1 << "=" << dv2-z1 << endl;
	getch();
	cout << dv2 << "*" << z1 << "=" << dv2*z1 << endl;
	cout << endl;
	getch();

	cout << z1 << "+" << v1 << "=" << z1+v1 << endl;
	getch();
	cout << z1 << "-" << v1 << "=" << z1-v1 << endl;
	getch();
	cout << z1 << "*" << v1 << "=" << z1*v1 << endl;
	cout << endl;
	getch();

	cout << v1 << "+" << z1 << "=" << v1+z1 << endl;
	getch();
	cout << v1 << "-" << z1 << "=" << v1-z1 << endl;
	getch();
	cout << v1 << "*" << z1 << "=" << v1*z1 << endl;
	cout << endl;
*/
	ZVector z2=dv2*2*c.conj();
/*
	cout << z1 << "+" << z2 << "=" << z1+z2 << endl; getch();
	cout << z1 << "-" << z2 << "=" << z1-z2 << endl; getch();
	cout << z1 << "*" << z2 << "=" << z1*z2 << endl; getch();*/

/*	cout << z2 << ".phase()" << "=" << z2.phase() << endl; getch();
	cout << z2 << ".real()" << "=" << z2.real() << endl; getch();
	cout << z2 << ".imag()" << "=" << z2.imag() << endl; getch();
	cout << z2 << ".mag()" << "=" << z2.mag() << endl; getch();
	cout << z2 << ".mod2()" << "=" << z2.mod2() << endl; getch();
	cout << z2 << ".conj()" << "=" << z2.conj() << endl; getch();
*/


	return 1;

}