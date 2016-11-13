// qpskawgn.cpp : Defines the entry point for the console application.
//


#include  <math.h>
#include "canal/canal.h"
#include "globaldef.h"
#include "tools/utilitis.h"

ZVector canal_nip(const ZVector &data, const ZVector &coef) {
	return data.conv_nip(coef);
}

DVector canal_nip(const DVector &data, const DVector &coef) {
	return data.conv_nip(coef);
}

ZVector canal_nip(const DVector &data, const ZVector &coef) {
	return data.conv_nip(coef);
}

ZVector canal_nip(const ZVector &data, const DVector &coef) {
	return data.conv_nip(coef);
}

/** coef puissance decroissant */
void canal(ZVector &data, const ZVector &coef) {
	data.fixedconv(coef);
}

/** coef puissance decroissant */
void canal(DVector &data, const DVector &coef) {
	data.fixedconv(coef);
}
