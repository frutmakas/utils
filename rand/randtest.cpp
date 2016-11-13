#include <iostream>
#include <tools/all.h>
#include <time.h>
using namespace std;

int main(){
	int starttime=time(NULL);
	const int MAX_SIZE=3000000;
	DVector d(MAX_SIZE);
	dRandGausStatePtr dp;
	dRandGausInit(&dp, 0,1.75);
	
	ZVector z(3000000);
	zRandGausStatePtr zp;
	zRandGausInit(&zp, 0,1.75);
	DVector hist(2000), hupper(2000);
	double hmin=-10.0, hmax=10.0,hstep=0.01,hc=hmin+hstep;
	for(int h=0;h<2000;h++,hc+=hstep) {
	  hupper.vect[h]=hc;
	}
	double s1 = 0 /*d.sum()*/,s2=0;
	DCplx zs1(0.0,0.0), /*d.sum()*/ zs2(0.0,0.0);
	int NBX=50;
	for (int xx=0;xx<NBX;xx++)  {
		cout <<".";cout.flush();
		dbRandGaus(&dp, d);
		zbRandGaus(&zp, z);
		for(int i=0;i<MAX_SIZE;i++) { 
		    if (d.vect[i]>=hmin&&d.vect[i]<hmax) {
		          int pos = (int)floor(1000.0+d.vect[i]*100);
		          while(pos<2000 && hupper.vect[pos]<d.vect[i]) {
                        pos++;
		          }
                  while(pos>1&&hupper.vect[pos-1]>d.vect[i]) pos--;
                  hist.vect[pos]+=1.0;
            }
			s1+=d.vect[i]; 
			s2+=d.vect[i]*d.vect[i];
			zs1+=z.vect[i]; 
			zs2+=z.vect[i]*z.vect[i].conj();
		}
	}
	cout << endl;
	
	cout << "somme = " << s1 << endl;
	cout << "somme2 = " << s2 << endl;
	double dmoyenne=s1/((double)MAX_SIZE*NBX);
	cout << "moyenne = " << dmoyenne<<endl;
	cout << "variance = " << (double)(s2)/(NBX*(double)MAX_SIZE)-dmoyenne*dmoyenne << endl;
	cout << "variance parametré : 1.25" <<endl;
	
	cout << "somme = " << zs1 << endl;
	cout << "somme2 = " << zs2 << endl;
	DCplx zmoyenne=zs1/((double)MAX_SIZE*NBX);
	cout << "moyenne = " << zmoyenne << endl; 
	cout << "variance = " << (zs2)/((double)MAX_SIZE*NBX)-zmoyenne*zmoyenne.conj() << endl << endl;
	cout << "variance parametré : 1.75" <<endl;
	
	cout << "Calculs fait sur " << ((double)MAX_SIZE*NBX) << " echantillons en " << time(NULL)-starttime <<" secondes" <<endl;
	ofstream os("matlab_rand.m");
	matlaboutput(os, "h7", hist, 15);

	os.close();
	system("pause");
	return 0; 
}

int main2() {
    ZVector toto(300);
	zRandGausStatePtr zp;
	zRandGausInit(&zp, 0,1.75);
	cout << zbRandGausApply(&zp, toto) << endl;
	system("pause");
   return 0; 
}
