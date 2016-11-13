#include <iostream.h>
#include <fstream.h>
#include <tools/all.h>
#include <fading/fadingtdma.h>
#include <rand/randgen.h>
#include <time.h>

int main() {
	dRandUniStatePtr tapseed;
	dRandUniInit(&tapseed, time(NULL), 0,1);
	ZMatrix taps = fadingtaps(1500,1,3, 3000, 1024, &tapseed);
	ofstream of("tapseed.txt");
	of << taps;
	of.close();
	return 1;
}