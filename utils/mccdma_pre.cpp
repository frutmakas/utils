/********************************************************************
 *                        File Information                          *
 ********************************************************************
 * $Name:  $
 * $Author: syed $
 * $Date: 2004/01/19 23:57:50 $
 * $Revision: 1.1.2.11 $
 * $Id: mccdma_pre.cpp,v 1.1.2.11 2004/01/19 23:57:50 syed Exp $
 ********************************************************************/
#include <iostream.h>
#include <all.h>
#include <gold/gold.h>
#include <time.h>
#include <conio.h>
int main() {

	int InitialValue[2]={31,31};
	int Polynomial[2]={18, 27};

	int nb_user_max = 3, nb_symbol=200, nb_tx=2, nb_rx=1/*, nb_taps=3*/, nb_carrier=31;
	cout << "NBuser max = "; cout.flush();
	cin >> nb_user_max;
	cout << "nb symbol= "; cout.flush();
	cin >> nb_symbol;

	dRandUniStatePtr cseed;
	int seed = time(NULL);
//	seed= 1074336888;
	dRandUniInit(&cseed, seed, 0, 1);

	char filename[2500];
	sprintf(filename, "%ld_mmsecdma_res.txt", seed);
	ofstream xxx(filename);
	xxx << " seed = " << seed << endl;
	ZMatrix theeye(nb_carrier,nb_carrier);
	theeye.eye();
	for(int nb_user=1;nb_user<=nb_user_max; nb_user++) {	
		xxx << "-----------------------------------          " <<endl<< nb_user << "-----------------------------------           " <<endl;
		cout << "-----------------------------------          " <<endl << nb_user << "-----------------------------------           " <<endl;
		int code_size=0, spreadsize=0;
		UVector<DVector> gold(nb_user);
		TRandBitStatePtr bseed;
		UVector<WVector> userbit(nb_user);
		UVector<ZVector> usersymbol(nb_user);
		UVector<ZVector> userspreaded(nb_user);


		typedef UMatrix<ZMatrix> TCanalSTBC; //(nb_tx, nb_rx);
		typedef UVector<ZMatrix>TTimedCanal; // mot stbc
		typedef UMatrix<TTimedCanal> TTimedCanalSTBC;
		
		UVector<TCanalSTBC> canal(nb_user);
		UVector<TTimedCanalSTBC> usertimedchannel(nb_user);
 
		int GoldIDX =(int)ceil(log2(nb_carrier));
		ZVector tmp = dGoldCode(InitialValue, Polynomial, GoldIDX, GoldIDX , 0);//, &code_size);
		nb_carrier=tmp.taille;
		cout << " nb carrier : " << nb_carrier << endl;
		//construction de matrice Q
		ZMatrix Q(nb_carrier, nb_carrier);
		UVector<ZMatrix> r2p(nb_user), r2p1(nb_user);
		typedef UMatrix<DVector> TChannelProfile;
		UVector<TChannelProfile> channelprofile(nb_user);

		double qnorm = 1.0; // /sqrt(nb_carrier);
		int i;
		for( i =0;i<nb_carrier;i++) {
			for(int j=0;j<nb_carrier;j++) {
				Q.mat[i][j] = qnorm*DCplx::exp2RI(1, -M_2_PI*i*j/nb_carrier);
			}
		}

		long max_doppler=10, filter_length=5, nb_rays=8192;

		double  symbol_time=1;
		ZMatrix alluser_2p, alluser_2p1;
		ZMatrix theeye = ZMatrix::eye(nb_carrier);
		for ( i=0;i<nb_user;i++) {
			channelprofile.vect[i] = TChannelProfile(nb_tx, nb_rx);
			gold.vect[i]=dGoldCode(InitialValue, Polynomial, GoldIDX, GoldIDX, i);//, &code_size);
			nb_carrier=gold.vect[i].taille;
			userbit.vect[i] = WVector(3*nb_symbol*nb_tx);
			bRandBit1(&bseed, userbit.vect[i]);
			usersymbol.vect[i]=(*psk8.modulate_function)(userbit.vect[i]);
			userspreaded.vect[i] = SpreadData(usersymbol.vect[i], gold.vect[i]);
			// on a les sequences de l'utilisateur i a transmettre .. 
			// on va generer le canal vue de chaque antenne de l'utilisateur i
			// on aura nb_symbol*nb_tx bloc a traiter .. 
			// on aura nb_symbols mot stbc

			canal.vect[i]=TCanalSTBC(nb_tx, nb_rx);
			usertimedchannel.vect[i] = TTimedCanalSTBC(nb_tx, nb_rx);
			for(int t=0;t<nb_tx;t++) {
				for(int r=0;r<nb_rx;r++) {
					channelprofile.vect[i].mat[t][r] = DVector(filter_length);
					dbRandUni(&cseed, channelprofile.vect[i].mat[t][r]);
					channelprofile.vect[i].mat[t][r].normalizepower();
					channelprofile.vect[i].mat[t][r]*=0.5;
					canal.vect[i].mat[t][r]=fadingtaps(max_doppler, symbol_time, filter_length, nb_symbol, nb_rays, &cseed);
					MultEachColinMatrixWithY(canal.vect[i].mat[t][r], channelprofile.vect[i].mat[t][r]);
					usertimedchannel.vect[i].mat[t][r] = TTimedCanal(nb_symbol);
				}
			}
			// on a le canal, on a les symboles .. il faut maintenant le transmettre
			// on travaille sur des blocs de nb_carrier
			// le canal reste invariant pdt un mot stbc
			r2p.vect[i] = ZMatrix(nb_symbol*nb_carrier,1);
			r2p1.vect[i] = ZMatrix(nb_symbol*nb_carrier,1);
			for(int ht=0, st=0;ht<nb_symbol; ht++) {
				// extraire LAMDA1 et LAMDA2
				ZMatrix xlamba1(nb_carrier, 1), xlamba2(nb_carrier, 1);
				for(int l=0;l<filter_length;l++) {
					xlamba1.mat[l][0] = canal.vect[i].mat[0][0].mat[ht][l];
					xlamba2.mat[l][0] = canal.vect[i].mat[1][0].mat[ht][l];
				}
#define useeye
#ifdef useeye
				usertimedchannel.vect[i].mat[0][0].vect[ht]=theeye;//ColXinMatrixAsDiagonal(Q*xlamba1, 0);
				usertimedchannel.vect[i].mat[1][0].vect[ht]=theeye;//ColXinMatrixAsDiagonal(Q*xlamba2, 0);
#else 
				usertimedchannel.vect[i].mat[0][0].vect[ht]=ColXinMatrixAsDiagonal(Q*xlamba1, 0);
				usertimedchannel.vect[i].mat[1][0].vect[ht]=ColXinMatrixAsDiagonal(Q*xlamba2, 0);
#endif
				// on a les lamdas .. 
				// cherchons X1 et X2 ... 
				ZMatrix X1 = VectorInColMatrixRepresentation(userspreaded.vect[i].copy(st, st+nb_carrier));
				st+=nb_carrier;
				ZMatrix X2 = VectorInColMatrixRepresentation(userspreaded.vect[i].copy(st, st+nb_carrier));
				st+=nb_carrier;

				// j'ai mon canal + mes mots stbc .. 
				// je vais obtenir le signal recu pour chaque mot mccdma
				InsertSubMatrixIntoMatrix(r2p.vect[i],  usertimedchannel.vect[i].mat[0][0].vect[ht]*X1 - usertimedchannel.vect[i].mat[1][0].vect[ht]*X2.conj(), ht*nb_carrier, 0);
				InsertSubMatrixIntoMatrix(r2p1.vect[i], usertimedchannel.vect[i].mat[0][0].vect[ht]*X2 + usertimedchannel.vect[i].mat[1][0].vect[ht]*X1.conj(), ht*nb_carrier, 0);
			}
			// on a traiter tout les symboles .. 
			if(i==0) { 
				alluser_2p = r2p.vect[i];
				
				alluser_2p1 = r2p1.vect[i].conj();
			}else {
				alluser_2p += r2p.vect[i];
				alluser_2p1 += r2p1.vect[i].conj();
			}
			// j'ai ma superposition des symboles .. 
		}

		canal.~UVector();
		r2p.~UVector();
		r2p1.~UVector();

		// maintenant appliquons le filtre de wiener
		// d=w.hermit*r
		// w = inv(E[rr*])*rho
		// rho = cµ * (sigma2_2p*LAMDA1µ + sigma2_2p1*LAMDA2µ)

		// canal bouge lentement ==> E[rr*] = rr*
		// decodage utilisateur µ
		int mu = 0;
		double sigmaX2p =1.0;
		double sigmaX2p1 =1.0;
		WVector dec(userbit.vect[mu].taille);
		int decptr=0;
		for(int t=0;t<nb_symbol;t++) {
			// traitement de temp 2p
			ZMatrix r2p = SubMatrix(alluser_2p, t*nb_carrier, 0, nb_carrier, 1);
			ZMatrix Err_2p = (r2p*r2p.hermit()).inv_nip();
			cout << r2p;
			cout << Err_2p;
			cout << sigmaX2p*usertimedchannel.vect[mu].mat[0][0].vect[t]+sigmaX2p1*usertimedchannel.vect[mu].mat[1][0].vect[t];
			ZVector rho2p = gold.vect[mu]*(sigmaX2p*usertimedchannel.vect[mu].mat[0][0].vect[t]+sigmaX2p1*usertimedchannel.vect[mu].mat[1][0].vect[t]);
			ZVector w2p = Err_2p*rho2p;
			cout << w2p;

			ZVector beta2p = w2p.conj()*r2p;

			// traitement de temp 2p+1
			ZMatrix r2p1 = SubMatrix(alluser_2p1, t*nb_carrier, 0, nb_carrier, 1);
			ZMatrix Err_2p1 = (r2p1*r2p1.hermit()).inv_nip();
			ZVector rho2p1 = gold.vect[mu]*(sigmaX2p*usertimedchannel.vect[mu].mat[0][0].vect[t]+sigmaX2p1*usertimedchannel.vect[mu].mat[1][0].vect[t]);
			ZVector w2p1 = Err_2p1*rho2p1;

			ZVector beta2p1 = w2p1.conj()*r2p1;

			// extraire symbol
			ZVector X2p = 0.5*(beta2p+beta2p1);
			dec.insert((*psk8.demodulate_function)(X2p), decptr);
			decptr+=3;

			ZVector X2p1= 0.5*(beta2p1-beta2p);
			dec.insert((*psk8.demodulate_function)(X2p1.conj()), decptr);
			decptr+=3;
			// done !
		}

		cout << userbit.vect[mu] << dec << userbit.vect[mu] - dec <<endl;
	}

	return 0;
}

