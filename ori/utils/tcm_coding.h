#ifndef TCM_CODING_H
#define TCM_CODING_H
double       log_add(double a,double b);
long_matrice back_ward_treillis(long_matrice FW_treillis);
long_vector  calcul_tail(long etat_fin, long_matrice trans);
long_vector  tcm_encoder(long_vector &x, long_matrice trans, long_matrice output, long fermeture);
cmplx_vect   mapper_8PSK(long_vector x);
cmplx_vect   mapper(long_vector x, cmplx_vect scat);
cmplx_vect   scat_8PSK();
cmplx_vect   scat_Shift8PSK();
long_vector  symbol_interleaver(long_vector in, long_vector intrl,long );
void         sym_interleaver(matrice in, long_vector intrl);
long_vector  symbol_deinterleaver(long_vector in,long_vector intrl,long );
cmplx_vect   symbol_deinterleaver(cmplx_vect in,long_vector intrl,long );
void         sym_deinterleaver(matrice in, long_vector intrl);
void         swap_selector(cmplx_vect a, cmplx_vect b, cmplx_vect out);
void         swap_combiner(cmplx_vect r, cmplx_vect a, cmplx_vect b);
void         TCM_decoder_SIHO(cmplx_vect x, long_vector y, long_matrice trans, long_matrice output,
                   cmplx_vect scat, long close_switch);
void         TCM_log_MAP_decoder(cmplx_vect x, matrice app, matrice LLR, long_matrice trans,
                   long_matrice BW_trans,long_matrice output,cmplx_vect scat, double sigma2,long fermeture);
void         Turbo_TCM_I(cmplx_vect r, matrice app, long_vector bit_dec, long_vector intrl,long_matrice trans,
				long_matrice BW_trans, long_matrice output, cmplx_vect scat, double sigma2);
#endif