#include <time.h>
#include <iostream.h>
#include <tools/all.h>
#include <rand/randgen.h>

int main(void) {

	int iseed=17;
	TRandBitStatePtr bseed;
	bseed.iseed = iseed;
	ofstream m("mranddebug.m");
	WVector b(200);
	bRandBit1(&bseed, b);
	matlaboutput(m, "bvect", b, 17);
	dRandUniStatePtr dseed;
	dRandUniInit(&dseed, iseed, 0, 1);
	DVector d(2000);
	dbRandUni(&dseed, d);
	matlaboutput(m, "dvect", d, 17);

	dRandGausStatePtr dgseed;
	dRandGausInit(&dgseed, 0, 1);
	DVector dg(2000);
	dbRandGaus(&dgseed, dg);
	matlaboutput(m, "dgvect", dg, 17);

	zRandGausStatePtr zgseed;
	zRandGausInit(&zgseed, 0, 1);
	ZVector zg(2000);
	zbRandGaus(&zgseed, zg);
	matlaboutput(m, "zgvect", zg, 17);

	return 1;
}
	
