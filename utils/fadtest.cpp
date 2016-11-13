#include <time.h>
#include <iostream.h>
#include <tools/all.h>
#include <rand/randgen.h>
#include <fading/fadingtdma.h>

int main(void) {

	int iseed=17;
	TRandBitStatePtr bseed;
	bseed.iseed = iseed;
	ofstream m("mranddebug.m");
	dRandUniStatePtr dseed;
	dRandUniInit(&dseed, iseed, 0, 1);
	
	DVector fd = freqd(1024, 150.0, &dseed);
	matlaboutput(m, "fd", fd, 15);

	DMatrix fg = fgwave(1024, 1, 1.0, &dseed);
	matlaboutput(m, "fg", fg, 15);
	m.close();
	/*
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
*/
    getchar();
	return 1;
}
	
