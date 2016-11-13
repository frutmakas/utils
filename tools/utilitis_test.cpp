#include <tools/all.h>
#include <iostream.h>

int main(){
	// test DCPlx
	DCplx c1, c2(1.0,2.1);
	c1=10*c2*0.1;
	DCplx c3, c4, c5, c6, c7, c8, c9(c1);
	c3=c1+c2;
	c4=c3.conj();
	c5=c3-c2;
	c6 = c1/0.5;
	c7=c3;
	c8=c1;
	c7+=3.0;
	c8+=c7;


	cout.precision(15);
	cout << "DCplx result test ... " << endl;
	cout << "c1, c2, c5, c9 : 1.0+i2.1 <==> " << c1 << " , " << c2<< " , " << c5 << " , " << c9 << endl;
	cout << "c3, c6 : 2.0+i4.2 <==> " << c3 << " , " << c6 << endl;
	cout << "c4     : 2.0-i4.2 <==> " << c4 << endl;
	cout << "c7     : 5.0+i4.2 <==> " << c7 << endl;
	cout << "c8     : 6.0+i6.3 <==> " << c8 << endl;
	cout << "mag(c8): 8,7 <==> " << c8.mag() << endl;
	cout << "mod2(c8): 75.69 <==> " << c8.mod2() << endl;
	cout << "c8!=c7 : true <==> " << (c8!=c7) << endl;
	cout << "c3!=c6 : false <==> " << (c6!=c3) << endl;
	cout << "c3*c4  : 21.64 <==> " << c3*c4 << endl;
	c1=3;
	c2=c2+c2.conj();
	c3*=c1;
	c4=c3;
	c4*=3;
	c5=c3+3;
	c6 = c5-3;
	c7-=c3;
	c8-=1;
	cout << " c1 : 3 <==> " << c1 << endl;
	cout << " c2 : 2 <==> " << c2 << endl;
	cout << " c3 : 6.0+i12.6 <==> " << c3 << endl;
	cout << " c4 : 18.0+i37.8 <==> " << c4 << endl;
	cout << " c1!=3 : 0<==> " << (c1!=3) << endl;
	cout << " c5 : 9.0+i12.6 <==> " << c5 << endl;
	cout << " c6 : 6.0+i12.6 <==> " << c6 << endl;
	cout << " -c5 : -9.0-i12.6 <==> " << -c5 << endl;
	cout << "c7     : -1.0-i8.4 <==> " << c7 << endl;
	cout << "c8     : 5.0+i6.3 <==> " << c8 << endl;
	cout << " c1==c2 : false <==> " << (c1==c2) << endl;
	cout << " c1==c1 : true <==> " << (c1==c1) << endl;
	cout << " c1==3 : true <==> " << (c1==3) << endl;
	cout << " c1==2 : true <==> " << (c1==3) << endl;
	cout << " c8==5 : false <==> " << (c8==5) << endl;
	c1=DCplx(30,30);
	c2=DCplx(30,-30);
	c3=DCplx(-30,30);
	c4=DCplx(-30,-30);
	cout << " c1.phase() : 45°/0,785398163 <==> " << c1.phase() << endl;
	cout << " c2.phase() : -45°/-0,785398163 <==> " << c2.phase() << endl;
	cout << " c3.phase() : 135°/2,356194490 <==> " << c3.phase() << endl;
	cout << " c4.phase() : -135°/-2,356194490 <==> " << c4.phase() << endl;









return 0;

	// test DMatrix
	const int SIZE_L=11, SIZE_C=10;
	DMatrix a(SIZE_L,SIZE_C) , b(SIZE_L, SIZE_C);
	for(int i=0;i<SIZE_L;i++) {
		for(int j=0;j<SIZE_C;j++) {
			a.mat[i][j]=i*SIZE_C+j;
		}
	}
	b = 1.1*a-5.7*a; 
	DMatrix c = a*a.transpose_nip(), d = a.transpose_nip()*a;
	c.transpose();
	cout << a << b << endl;
	cout << c << d << endl;
	cout << b/(1.1-5.7) - a << endl;
	DMatrix e(1,10,0.55557);
	cout << e << endl;
	a = e;
	DMatrix f(a);
	cout << a << f << endl;
	return 1;

}