#ifndef CODING_H
#define CODING_H
#include "tools/utilitis.h"

long_vector codeur_rsc_2(long_vector x, long_vector g, long k, long param);
void backward_treillis_rsc
	(long_vector g, long k, long_matrice &last_state, long_matrice &last_out);
void forward_treillis_rsc
		(long_vector g,long k, long_matrice &next_state, long_matrice &next_out);
vector sova_paquet(vector sig_rec, long info_len, long n, long k, vector L_a,
			 long_matrice last_state, long_matrice last_out, long ind_dec);
vector matrix_interleaver(vector a,long lin, long col,long param);
long_vector matrix_interleaver(long_vector a,long lin, long col, long param);

#endif

