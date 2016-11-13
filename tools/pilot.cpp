/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: jebon $
 * $Date: 2003/03/31 14:25:59 $
 * $Revision: 1.1.2.3 $
 * $Id: pilot.cpp,v 1.1.2.3 2003/03/31 14:25:59 jebon Exp $
 ********************************************************************/
#include "tools/utilitis.h"
#include <math.h>

ZVector insert_pilot(const ZVector &pilot, const ZVector &data, int symbols_between_pilot) {
	int nb_block = (int)ceil((double)data.taille/symbols_between_pilot);
	int outsize = nb_block*pilot.taille+data.taille;
	ZVector output(outsize);
	for(int sym=0, out=0;sym<data.taille;sym+=symbols_between_pilot) {
		output.insert(pilot, out);
		out+=pilot.taille;
		output.insert(data.copy(sym, sym+symbols_between_pilot), out);
		out+=symbols_between_pilot;
	}
	return output;
}

ZVector delete_pilot(const ZVector data, int pilot_size, int symbols_between_pilot) {
	int nb_block=(int)floor((double)(data.taille+symbols_between_pilot)/(symbols_between_pilot+pilot_size));
	int outsize=data.taille-pilot_size*nb_block;
	ZVector output(outsize);
	for(int i=pilot_size,j=0;i<data.taille;i+=symbols_between_pilot+pilot_size, j+=symbols_between_pilot) {
		output.insert(data.copy(i,i+symbols_between_pilot),j);
	}
	return output;
}

ZVector extract_pilot(const ZVector data, int pilot_size, int symbols_between_pilot) {
	int nb_block=(int)floor((double)(data.taille+symbols_between_pilot)/(symbols_between_pilot+pilot_size));
	int outsize=pilot_size*nb_block;
	ZVector output(outsize);
	for(int i=0,j=0;i<data.taille;i+=symbols_between_pilot+pilot_size, j+=pilot_size) {
		output.insert(data.copy(i,i+pilot_size),j);
	}
	return output;
}

int separate_pilot(const ZVector data, int pilot_size, int symbols_between_pilot, ZVector &pilot, ZVector &raw_data) {
	int nb_block=(int)floor((double)(data.taille+symbols_between_pilot)/(symbols_between_pilot+pilot_size));
	int pilot_outsize=pilot_size*nb_block;
	int outsize=data.taille-pilot_size*nb_block;
	raw_data = ZVector(outsize);
	pilot = ZVector(pilot_outsize);
	for(int i=0,j=0, k=0;i<data.taille; j+=pilot_size, k+=symbols_between_pilot) {
		pilot.insert(data.copy(i,i+pilot_size),j);
		i+=pilot_size;
		raw_data.insert(data.copy(i,i+symbols_between_pilot),k);
		i+=symbols_between_pilot;
	}
	return nb_block;
}

