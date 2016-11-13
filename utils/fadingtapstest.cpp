#include <all.h>
#include <iostream.h>
#include <fstream.h>

int main() {

	dRandUniStatePtr cseed;
	//ofstream os("retestfadingtaps.txt", ios::app);
	for(int i=0;i<2;i++) {
		cout << i << endl;
	
		dRandUniInit(&cseed, 87000, 0, 1);

		fadingtaps(10,1,5,1000,8192,&cseed);
	}
//	os.close();
	return 1;

}