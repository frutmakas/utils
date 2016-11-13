/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/01/20 13:39:52 $
 * $Revision: 1.1.2.1 $
 * $Id: spread.cpp,v 1.1.2.1 2004/01/20 13:39:52 syed Exp $
 ********************************************************************/

#include <tools/utilitis.h>

ZVector SpreadData(const ZVector &data, const ZVector &spreadcode){
	ZVector spreaded(data.taille*spreadcode.taille);

	for(int i=0;i<data.taille;i++) {
		spreaded.insert(spreadcode*data.vect[i], i*spreadcode.taille);
	}
	return spreaded;
}

ZVector SpreadData(const DVector &data, const ZVector &spreadcode){
	ZVector spreaded(data.taille*spreadcode.taille);

	for(int i=0;i<data.taille;i++) {
		spreaded.insert(spreadcode*data.vect[i], i*spreadcode.taille);
	}
	return spreaded;
}

ZVector SpreadData(const ZVector &data, const DVector &spreadcode){
	ZVector spreaded(data.taille*spreadcode.taille);

	for(int i=0;i<data.taille;i++) {
		spreaded.insert(spreadcode*data.vect[i], i*spreadcode.taille);
	}
	return spreaded;
}

DVector SpreadData(const DVector &data, const DVector &spreadcode){
	DVector spreaded(data.taille*spreadcode.taille);

	for(int i=0;i<data.taille;i++) {
		spreaded.insert(spreadcode*data.vect[i], i*spreadcode.taille);
	}
	return spreaded;
}

DVector UnspreadData(const DVector &data, const DVector &spreadcode){
	if (data.taille==0 || spreadcode.taille==0 || data.vect==NULL || spreadcode.vect==NULL) return DVector();
	DVector despread((int)(ceil((double)(data.taille)/(double)(spreadcode.taille))));
	for(int i=0; i<despread.taille;i++) {
		double integral =0.0; int limit=(i+1)*spreadcode.taille;
		for(int j=i*spreadcode.taille, k=0;j<limit;j++,k++) {
			integral+=data.vect[j]*spreadcode.vect[k];
		}
		despread.vect[i]=integral;
	}
	return despread;
}


ZVector UnspreadData(const ZVector &data, const ZVector &spreadcode){
	if (data.taille==0 || spreadcode.taille==0 || data.vect==NULL || spreadcode.vect==NULL) return ZVector();
	ZVector despread((int)(ceil((double)(data.taille)/(double)(spreadcode.taille))));
	for(int i=0; i<despread.taille;i++) {
		DCplx integral; integral  =0.0; int limit=(i+1)*spreadcode.taille;
		for(int j=i*spreadcode.taille, k=0;j<limit;j++,k++) {
			integral+=data.vect[j]*spreadcode.vect[k];
		}
		despread.vect[i]=integral;
	}
	return despread;
}

ZVector UnspreadData(const ZVector &data, const DVector &spreadcode){
	if (data.taille==0 || spreadcode.taille==0 || data.vect==NULL || spreadcode.vect==NULL) return ZVector();
	ZVector despread((int)(ceil((double)(data.taille)/(double)(spreadcode.taille))));
	for(int i=0; i<despread.taille;i++) {
		DCplx integral; integral=0.0; int limit=(i+1)*spreadcode.taille;
		for(int j=i*spreadcode.taille, k=0;j<limit;j++,k++) {
			integral+=data.vect[j]*spreadcode.vect[k];
		}
		despread.vect[i]=integral;
	}
	return despread;
}

ZVector UnspreadData(const DVector &data, const ZVector &spreadcode){
	if (data.taille==0 || spreadcode.taille==0 || data.vect==NULL || spreadcode.vect==NULL) return ZVector();
	ZVector despread((int)(ceil((double)(data.taille)/(double)(spreadcode.taille))));
	for(int i=0; i<despread.taille;i++) {
		DCplx integral; integral=0.0; int limit=(i+1)*spreadcode.taille;
		for(int j=i*spreadcode.taille, k=0;j<limit;j++,k++) {
			integral+=data.vect[j]*spreadcode.vect[k];
		}
		despread.vect[i]=integral;
	}
	return despread;
}

