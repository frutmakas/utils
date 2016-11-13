 /********************************************************************
  *                        File Information                          *
  ********************************************************************
  * $Name:  $
  * $Author: syed $
  * $Date: 2003/11/12 14:36:08 $
  * $Revision: 1.1.2.13 $
  * $Id: canal_ofdm.cpp,v 1.1.2.13 2003/11/12 14:36:08 syed Exp $
  ********************************************************************/

#include "tools/all.h"
#include "tools/superclass.h"
#include "tools/ext_utilitis.h"
#include "tools/ext_matrix.h"
#include "canal/canal_ofdm.h"
#include "rand/randgen.h"
#include "modulation/psk.h"
#include "stbc/conversion.h"
#include "fading/fadingtdma.h"

ZMatrix ofdm_time_frequency_channel_estimation(const ZMatrix &C, const ZMatrix &D, int training_length, int Lf,  double sigma, double fd /*max_doppler*/, ZMatrix &nC) {
	ZMatrix yob; //((training_length+1)*Lf,1);
	ZMatrix A(training_length*Lf, Lf);
	UVector<ZVector> uD = SplitMatrixByCol(D);
	UVector<ZVector> uC =  SplitMatrixByCol(C);

	register int i, j, k, c, l;
	yob = ConcatUVectorAsColumnMatrix(uC);
	for(k=0;k<training_length;k++){
		InsertSubMatrixIntoMatrix(A, diag(uD.vect[k]), k*Lf, 0);
	}

	ZMatrix ztmp=A.hermit()*A+sigma*ZMatrix::eye(Lf);
	ztmp.inv();

	ZMatrix h = ztmp*A.hermit()*yob;

	//% Calcul de la fonction d'autocorrélation fréquentielle, on la calcule à l'aide de la fonction xcorr
 	ZVector Sf = MatrixInVectorRepresentation(h).xcorr_nip();

	//% Construction de la matrice de Toeplitz de corrélation fréquentielle du
	//% canal

	ZVector Sf_i(Lf);
	for(j=0;j<Lf;j++) 
		Sf_i.vect[j] = Sf.vect[Lf-1-j];

	ZMatrix Af = gen_toeplitz(Sf_i, Sf.copy(Lf-1, 2*Lf-1));

	ZMatrix U,  V;
	DMatrix S;
	zsvd(Af, U, S, V);

	//% Calcul des coefficients de la matrice diagonale Phi(w)
	//% On utilise la distribution de Jakes pour le spectre Doppler
	DVector a(Lf), b(Lf), gg(Lf), f(Lf);
	for(int m=0;m<Lf;m++) {
		a.vect[m]=2*S.mat[m][m]/(M_2_PI*fd*sigma);
		if(a.vect[m]>=1){
			double temp = sqrt(a.vect[m]*a.vect[m]-1);
			b.vect[m]=M_PI_2*a.vect[m]-temp*(M_PI_2-asin(1.0/a.vect[m]));
		} else {
			double temp = sqrt(1-a.vect[m]*a.vect[m]);
			b.vect[m]=M_PI_2*a.vect[m]+temp*log((1+temp)/a.vect[m]);
		}
		gg.vect[m]=exp(-2*fd*(b.vect[m]+log(a.vect[m]/2)));
		f.vect[m]=1-gg.vect[m];
	}

	//% Construction de la matrice de l'estimateur fréquentiel 
	nC = diag(f);
	ZMatrix nD=U*nC*V.hermit()*h;
	ZMatrix CF(Lf,1);
	for(i=0;i<Lf;i++) {
		CF.mat[i][0]=(nD.mat[i][0]).conj()/nD.mat[i][0].mag();
	}

	return CF;
}

// 
// function  [TEB,SER] = MIMO1OFDM
DMatrix  MIMO1OFDM(DMatrix &TEB ) {
// % On considère une modulation MDP-4 transmise en OFDM %
// % sur deux antennes d'émission en utilisant le schéma d'Alamouti %
// % M=((s0, s1),(-s1* s0*)) %
// n = 100;
// L = 8;
// % taille d'un paquet en nombre de symboles MDP-4 %
// M = 2; 
// Nbpacket = 100;
// Lf = 8;
// % taille de la FFT: Lf doit être un multiple du nombre d'antennes d'émission %
// Nt = 2;
// % nombre d'antennes de transmission %
// training_length=5;
	int n=200, L=8, M=2, Nbpacket=1000, Lf=8, Nt=2, training_length=5;
// % frequence Doppler maximale dans le specte fréquentiel du canal %
// fd = 50;
	double fd=50.0;
// coding_rate = 1;
// nb_step = 10;
	int coding_rate=1, nb_step=10;
// Gn=[];
// Gn1=[];
// Gn2=[];
//#error "Please define type of Gn?"
	ZMatrix Gn(Lf*L,2*Lf),Gn1(Lf*L,Lf),Gn2(Lf*L,Lf);

// SNR = [1,2,3,4,5,6];
	int NB_SNR=6;
// for xx=1:length(SNR)
//     for jj=1:Nbpacket
//         num_error(xx,jj)=0;  
//     end;
// end;
	WMatrix num_error(NB_SNR, Nbpacket,0);

// for dd=1:length(SNR)
//     SER(dd,1)=0;
// end;
//	DVector SER(NB_SNR,0.0);
	// generating binary signals .. 
	WVector seq(n*M*Lf*Nt*Nbpacket);
	TRandBitStatePtr bptr;
	bRandBit1(&bptr, seq);
#define DEBUG__
#ifdef DEBUG__
	cout << " Debug : Signal modulation " << endl;
#endif
	ZVector vn = z4PSK_mod(seq);
#ifdef DEBUG__
	cout << " Debug : STBC Signal splitting " << endl;
#endif
	UVector<ZVector> sig = SplitMatrixByCol(G2(vn, vn.taille));

#ifdef DEBUG__
	cout << " Debug : Antenna data stockage reformat" << endl;
#endif
	ZMatrix ant1 = ColumnPackedVectorAsMatrix(sig.vect[0], Lf, n*Nt*Nbpacket);
	ZMatrix ant2 = ColumnPackedVectorAsMatrix(sig.vect[1], Lf, n*Nt*Nbpacket);
	ZMatrix dcn =  ColumnPackedVectorAsMatrix(vn, Lf, n*Nt*Nbpacket);
	sig.~UVector<ZVector>();
	vn.~ZVector();
#ifdef DEBUG__
	cout << " Debug : creating channel model " << endl;
#endif

	dRandUniStatePtr path1seed, path2seed;
	dRandUniInit(&path1seed, 11224567, 0,1);
	dRandUniInit(&path2seed, 78852248, 0,1);
	ZMatrix path1=fadingtaps(fd,1,1,Nbpacket*n*Nt,512, &path1seed);
	ZMatrix path2=fadingtaps(fd,1,1,Nbpacket*n*Nt,512, &path2seed);
	//ZMatrix ofdmpath1(Lf, Nbpacket*n*Nt);
	//ZMatrix ofdmpath2(Lf, Nbpacket*n*Nt);
	ZMatrix ofdmsig1(Lf, Nbpacket*n*Nt);
	ZMatrix ofdmsig2(Lf, Nbpacket*n*Nt);
	ZMatrix resultingsig(Lf, Nbpacket*n*Nt);
	zRandGausStatePtr noiseptr;
	ZMatrix noise(Lf, Nbpacket*n*Nt);

	ZVector angle(Lf);
#ifdef DEBUG__
	cout << " Debug : Preparing ofdm signal " << endl;
#endif

	for(int k=0;k<Lf;k++) {
		angle.vect[k]=DCplx::exp2RI(1, -M_2_PI*k/Lf);
		for(int ii=0;ii<n*Nt*Nbpacket;ii++) {
			//ofdmpath1.mat[k][ii]=path1.mat[ii][0]*angle.vect[k];
			//ofdmpath2.mat[k][ii]=path2.mat[ii][0]*angle.vect[k];
			ofdmsig1.mat[k][ii]=path1.mat[ii][0]*angle.vect[k]*ant1.mat[k][ii];
			ofdmsig2.mat[k][ii]=path2.mat[ii][0]*angle.vect[k]*ant2.mat[k][ii];
		}
	}

	path2.~ZMatrix();
	path1.~ZMatrix();

//   % Création du signal X1,n.H1,n %
//   for i=1:n*Nt
//        dd1(:,i)=diag(w1(:,i))*zy1(:,i);
//        dd2(:,i)=diag(w2(:,i))*zy2(:,i);
//        [d(:,i),std]=awgn(coding_rate,Lf,SNR(da),dd1(:,i)+dd2(:,i));
//    end;

// 
	const int SNR_MAX = 6;
	TEB=DMatrix(SNR_MAX, Nbpacket);
// for da=1:length(SNR)
	for(int snr_sweep=0;snr_sweep<SNR_MAX; snr_sweep++) {
		zRandGausInit(&noiseptr ,0, convertsnrtosigma2(coding_rate, snr_sweep));
		
		zmRandGaus(&noiseptr, noise);
		resultingsig = ofdmsig1+ofdmsig2+noise;

//	    for num_packet=1:Nbpacket
		for(int num_packet=0;num_packet<Nbpacket;num_packet++) {
			int packet_offset = n*Nt*num_packet;
			cout << " Debug process : SNR = " << snr_sweep << " and num_packet = " << num_packet << endl;

//         feet=n*M*Lf*Nt;
//         seq(:,num_packet)=(rand(1,feet)>0.5)';
//         mseq=seq(:,num_packet);
//         
//         for i=1:(n*Lf*Nt)
//             sn=mseq(M*(i-1)+1:M*i);
//             if(sn==[1 1]')
//                 vn(i)=exp(j*pi/4);
//             elseif(sn==[0 1]')
//                 vn(i)=exp(j*(pi/4+pi/2));
//             elseif(sn==[0 0]')
//                 vn(i)=exp(j*(pi/4+pi));
//             else
//                 vn(i)=exp(j*(pi/4+3*pi/2));
//             end;
//         end;
//         
//         % création du codage spatio-temporel % 
//         for jj=1:(n*Lf)
//             wn(Nt*(jj-1)+1:Nt*jj)=vn(Nt*(jj-1)+1:Nt*jj);
//             tn(Nt*(jj-1)+1)=-conj(vn(Nt*jj));
//             tn(Nt*jj)=conj(vn(Nt*(jj-1)+1));
//             s1(Nt*(jj-1)+1)=vn(Nt*(jj-1)+1);  
//             s1(Nt*jj)=tn(Nt*(jj-1)+1);
//             s2(Nt*(jj-1)+1)=vn(Nt*jj);
//             s2(Nt*jj)=tn(Nt*jj);
//         end;
//         
//         % creation de la modulation multiporteuse OFDM %
//         for i=1:n*Nt
//             w1(:,i)=conj((s1(Lf*(i-1)+1:Lf*i))'); // ant1
//             w2(:,i)=conj((s2(Lf*(i-1)+1:Lf*i))'); // ant2
//             dcn(:,i)=conj((wn(Lf*(i-1)+1:Lf*i))');
//         end;
//         
//         % On a à créer deux canaux de propagation sur chacun des deux flux %
//         % de sortie des deux antennes %
//         
//         % première antenne %
//         x1=fadingtaps(fd,1,1,n*Nt,512);
//         % deuxième antenne %
//         x2=fadingtaps(fd,1,1,n*Nt,512);
//         
//         % On construit les répliques fréquentielles décalées de ce trajet multiple
//         % On construit Lf colonnes chacune ayant Ns (nombre symboles) lignes 
//         % Ainsi, on obtient une matrice temps fréquence (colonnes => fréquences et 
//         % lignes => temps)
//         
//         for k=1:Lf
//             for ii=1:n*Nt
//                 y1(k,ii)=x1(ii);
//                 y2(k,ii)=x2(ii);  
//             end;
//             ww(k,1)=exp(-j*2*pi*((k-1)/Lf));
//             zy1(k,:)=y1(k,:).*ww(k,1); // ofdmpath1
//             zy2(k,:)=y2(k,:).*ww(k,1); //  // ofdmpath2
//         end; 
//         
//         % Création du signal X1,n.H1,n %
//         for i=1:n*Nt
//             dd1(:,i)=diag(w1(:,i))*zy1(:,i);
//             dd2(:,i)=diag(w2(:,i))*zy2(:,i);
//             [d(:,i),std]=awgn(coding_rate,Lf,SNR(da),dd1(:,i)+dd2(:,i));
//         end;

//         
//         for i=1:L
//             Gn=[Gn;diag(w1(:,i)),diag(w2(:,i))];
//             rn((i-1)*Lf+1:i*Lf,1)=d(:,i);
//             Gn1 = [Gn1;diag(w1(:,i))];   
//             Gn2 = [Gn2;diag(w2(:,i))]; 
//         end;
			ZMatrix rn(L*Lf,1);
			for(int i=0;i<L;i++) {
				for(int l=0;l<Lf;l++) {
					rn.mat[i*Lf+l][0]=resultingsig.mat[l][packet_offset+i];
					Gn1.mat[i*Lf+l][l]=Gn.mat[i*Lf+l][l] = ant1.mat[l][packet_offset+i];
					Gn2.mat[i*Lf+l][l]=Gn.mat[i*Lf+l][Lf+l] = ant2.mat[l][packet_offset+i];
				}
			}
			ofstream c2m("m:\\work\\cances.old\\ctom.m");
			matlaboutput(c2m, "rn", rn, 15);
			matlaboutput(c2m, "Gn", Gn, 15);
			matlaboutput(c2m, "Gn1", Gn1, 15);
			matlaboutput(c2m, "Gn2", Gn2, 15);
			c2m.close();

//         hest=pinv(Gn)*rn;
		   ZMatrix hest = pinv(Gn)*rn;

//         
//         % Test de l'EM et de l'algorithme SAGE %
//         % On initialise : algorithme EM %
//         
//         hest1EM = inv(Gn1'*Gn1)*Gn1'*rn;
//         hest2EM = inv(Gn2'*Gn2)*Gn2'*rn;
		   ZMatrix hest1EM = (Gn1.hermit()*Gn1).inv_nip()*Gn1.hermit()*rn;
		   ZMatrix hest2EM = (Gn2.hermit()*Gn2).inv_nip()*Gn2.hermit()*rn;
//         
		   c2m.open("m:\\work\\cances.old\\mmdebux.m");
		   matlaboutput(c2m, "xhest", hest,15);
		   matlaboutput(c2m, "xhest1EM", hest1EM,15);
		   matlaboutput(c2m, "xhest2EM", hest2EM,15);
		   c2m.close();

//         for ii = 1:nb_step
		   int ii;
		   for(ii=0;ii<nb_step;ii++) {
//             % On commence par calculer Zi,n (cf article Georghiadhes IEEE Trans. on %
//             % Commun. janvier 2003) %
//             Z1 = Gn1*hest1EM;
//             Z2 = Gn2*hest2EM;
			   ZMatrix Z1=Gn1*hest1EM;
			   ZMatrix Z2=Gn2*hest2EM;
				c2m.open("m:\\work\\cances.old\\mmdebux2.m");
				matlaboutput(c2m, "xZ1", Z1,15);
				matlaboutput(c2m, "xZ2", Z2,15);
			   c2m.close();

			   //             
//             Y1 = Z1 + 0.5*(rn-Z1-Z2);
//             Y2 = Z2 + 0.5*(rn-Z1-Z2);
			   ZMatrix Y1 = Z1+0.5*(rn-Z1-Z2);
			   ZMatrix Y2 = Z2+0.5*(rn-Z1-Z2);
				c2m.open("m:\\work\\cances.old\\mmdebux3.m");
				matlaboutput(c2m, "xY1", Y1,15);
				matlaboutput(c2m, "xY2", Y2,15);
			   c2m.close();
//             
//             hest1EM = inv(Gn1'*Gn1)*Gn1'*Y1;
//             hest2EM = inv(Gn2'*Gn2)*Gn2'*Y2;
			   hest1EM = (Gn1.hermit()*Gn1).inv_nip()*Gn1.hermit()*Y1;
			   hest2EM = (Gn2.hermit()*Gn2).inv_nip()*Gn2.hermit()*Y2;
//             ii = ii+1;
				c2m.open("m:\\work\\cances.old\\mmdebux4.m");
				matlaboutput(c2m, "xhest1EMi", hest1EM,15);
				matlaboutput(c2m, "xhest2EMi", hest2EM,15);
			   c2m.close();
//         end;
			}
//         % algorithme SAGE %


		   ZMatrix hest1sage = (Gn1.hermit()*Gn1).inv_nip()*Gn1.hermit()*rn;
		   ZMatrix hest2sage = (Gn2.hermit()*Gn2).inv_nip()*Gn2.hermit()*rn;

//         hest1sage = inv(Gn1'*Gn1)*Gn1'*rn;
//         hest2sage = inv(Gn2'*Gn2)*Gn2'*rn;

//         Z1 = Gn1*hest1sage;
//         Z2 = Gn2*hest2sage;
		   ZMatrix Z1 = Gn1*hest1sage;
		   ZMatrix Z2 = Gn2*hest2sage;
//         
//         for nn = 1:nb_step
		   for(int nn=0;nn<nb_step;nn++) {
//             % On commence par calculer Zi,n (cf article Georghiadhes IEEE Trans. on %
//             % Commun. janvier 2003) %
//             Y1 = rn-Z2;
			   ZMatrix Y1 = rn - Z2;
//             hest1sage = inv(Gn1'*Gn1)*Gn1'*Y1;
			   hest1sage = (Gn1.hermit()*Gn1).inv_nip()*Gn1.hermit()*Y1;
//             K1 = Gn1*hest1sage;
			   ZMatrix K1 = Gn1*hest1sage;
//             
//             Y2 = rn-Z1;
			   ZMatrix Y2 = rn-Z1;
//             hest2sage = inv(Gn2'*Gn2)*Gn2'*Y2;
			   hest2sage = (Gn2.hermit()*Gn2).inv_nip()*Gn2.hermit()*Y2;
//             K2 = Gn2*hest2sage;
			   ZMatrix K2 = Gn2*hest2sage;
//             Z1=K1;
			   Z1=K1;
//             Z2=K2;
			   Z2=K2;
//             ii = ii+1;
//         end;
		   }
//         
			c2m.open("m:\\work\\cances.old\\mmdebux5.m");
			matlaboutput(c2m, "xhest1sage", hest1sage,15);
			matlaboutput(c2m, "xhest2sage", hest2sage,15);
			c2m.close();

		   //         Gn=[];
//         Gn1=[];
//         Gn2=[];
		   DCplx ztmp(0,0);
		   Gn.reset(ztmp);
		   Gn1.reset(ztmp);
		   Gn2.reset(ztmp);
//         
//         % A partir de maintenant l'estimateur opère en mode DD: Decision Directed %
//         ddn(:,1:L)=dcn(:,1:L);
		   ZMatrix ddn(Lf, n*Nt);
           for(int al=0, bl=0;al<Lf;al++,bl++) {
                   for(int ac=0,bc=packet_offset;ac<L;ac++,bc++) {
                           ddn.mat[al][ac]=dcn.mat[bl][bc];
                   }
           }
//                   SubMatrix(dcn, 0, packet_offset, Lf, L);
		   c2m.open("m:\\work\\cances.old\\ctom2.m");
		   matlaboutput(c2m, "ddn", ddn, 15);
		   matlaboutput(c2m, "d", SubMatrix(resultingsig, 0, packet_offset, Lf, n*Nt), 15);
		   c2m.close();
			c2m.open("m:\\work\\cances.old\\mmdebux6.m");
			matlaboutput(c2m, "xddn", ddn,15);
			matlaboutput(c2m, "xresultingsig", resultingsig,15);
			c2m.close();

		   ZMatrix x(Lf,1);
		   ZMatrix w1(Lf, n*Nt), w2(Lf, n*Nt);
//         for ii=L+1:n*Nt
		   for(ii=L;ii<n*Nt;ii++) {
//             for k=1:Lf/Nt
			   ZMatrix H(2,2);
			   for(int k=0;k<Lf/Nt;k++) {
					cout << "Second ii loop : ii = " << ii << "  ,  k = " << k<< endl;

//                 H(1,:)=[hest(Nt*(k-1)+1),hest(Nt*(k-1)+1+Lf)];
//				   cout << hest << endl;
				   H.mat[0][0]=hest.mat[Nt*k][0];
				   H.mat[0][1]=hest.mat[Nt*k+Lf][0];
//                 H(2,:)=[conj(hest(Nt*k+Lf)),-conj(hest(Nt*k))];
				   H.mat[1][0]=hest.mat[Nt*(1+k)+Lf-1][0].conj();
				   H.mat[1][1]=-hest.mat[Nt*(k+1)-1][0].conj();
//                 A=conj([d((k-1)*Nt+1,ii),conj(d(k*Nt,ii))]');
				   ZMatrix A(2,1);
				   A.mat[0][0]=resultingsig.mat[k*Nt][packet_offset+ii];
				   A.mat[1][0]=resultingsig.mat[(k+1)*Nt-1][packet_offset+ii].conj();

//                 %x((k-1)*Nt+1:k*Nt)=pinv(H)*A;
//                 x((k-1)*Nt+1)=(H(1,2)*A(2)-H(2,2)*A(1))/(H(2,1)*H(1,2)-H(1,1)*H(2,2));
				   x.mat[k*Nt][0]=(H.mat[0][1]*A.mat[1][0]-H.mat[1][1]*A.mat[0][0])/(H.mat[1][0]*H.mat[0][1]-H.mat[0][0]*H.mat[1][1]);
//                 x(k*Nt)=(H(2,1)*A(1)-H(1,1)*A(2))/(H(1,2)*H(2,1)-H(1,1)*H(2,2));
				   x.mat[(k+1)*Nt-1][0]=(H.mat[1][0]*A.mat[0][0]-H.mat[0][0]*A.mat[1][0])/(H.mat[0][1]*H.mat[1][0]-H.mat[0][0]*H.mat[1][1]);
//                 ddn((k-1)*Nt+1:k*Nt,ii)=conj([cos(pi/4)*(sign(real(x((k-1)*Nt+1)))+j*sign(imag((x((k-1)*Nt+1))))),cos(pi/4)*(sign(real(x(k*Nt)))+j*sign(imag(x(k*Nt))))]');

				   ddn.mat[k*Nt][ii]=cos(M_2_PI/4.0)*DCplx(sign(x.mat[k*Nt][0].re),sign(x.mat[k*Nt][0].im));
//				                           cos(pi/4)*(sign(real(x(k*Nt)))+j*sign(imag(x(k*Nt))))]');
				   ddn.mat[(k+1)*Nt-1][ii]=cos(M_2_PI/4.0)*DCplx(sign(x.mat[(k+1)*Nt-1][0].re),sign(x.mat[(k+1)*Nt-1][0].im));
//                 w1((k-1)*Nt+1:k*Nt,ii)=conj([ddn((k-1)*Nt+1,ii),-conj(ddn(k*Nt,ii))]');
				   w1.mat[k*Nt][ii]=ddn.mat[k*Nt][ii];
				   w1.mat[(k+1)*Nt-1][ii]=-ddn.mat[(k+1)*Nt-1][ii].conj();
//                 w2((k-1)*Nt+1:k*Nt,ii)=conj([ddn(k*Nt,ii),conj(ddn((k-1)*Nt+1,ii))]');
				   w2.mat[k*Nt][ii]=ddn.mat[(k+1)*Nt-1][ii];
				   w2.mat[(k+1)*Nt-1][ii]=ddn.mat[k*Nt][ii].conj();
//             end; 
			   }

//             % Une fois que l'on a déterminé le vecteur des symboles décidés, on réestime le canal % 
//             % On recrée d'abord le codage STBC %
//             for pp=1:L

			   for(int pp=0;pp<L;pp++) {
				   for(int l=0;l<Lf;l++) {
					   cout << "Second ii loop : ii = " << ii << " ,  pp = " << pp  << "  , l = " << l << endl;
					   rn.mat[pp*Lf+l][0]=resultingsig.mat[l][packet_offset+pp+ii-L];
//                 rn((pp-1)*Lf+1:pp*Lf,1)=d(:,pp+ii-L);
//                 Gn=[Gn;diag(w1(:,pp+ii-L)),diag(w2(:,pp+ii-L))];
						Gn.mat[pp*Lf+l][l]=w1.mat[l][pp+ii-L];
						Gn.mat[pp*Lf+l][Lf+l]=w2.mat[l][pp+ii-L];
				   }
//             end;
			   }

			   c2m.open("m:\\work\\cances.old\\ctom3.m");
			   matlaboutput(c2m, "xrn", rn, 15);
			   matlaboutput(c2m, "xGn", Gn, 15);
			   c2m.close();

//             hest=pinv(Gn)*rn;
			   //cout << Gn << endl;
			   ofstream mfile("gndeb.m");
			   matlaboutput(mfile, "xGn", Gn, 15);
			   cout << "Second ii loop : ii = " << ii << " doing hest " << endl;
			   matlaboutput(mfile, "pinvGn", pinv(Gn), 15);

			   cout << " setp 2 ... " << endl;

    		   hest=pinv(Gn)*rn;
			   cout << "Second ii loop : ii = " << ii << " : Hest done ... "<< endl; 

			   matlaboutput(mfile, "hest", hest, 15);
			   mfile.close();

	//		   cout << hest << endl;
//             Gn=[];
			   Gn.reset(ztmp);
//         end;
		   }
		   ZMatrix err = SubMatrix(dcn, 0, packet_offset, Lf, n*Nt)-ddn;
//         % On compte les erreurs dans le paquet transmis %
//         for i=L+1:n*Nt
//             err(:,i)=fix(dcn(:,i)-ddn(:,i));
//         end;
//         for i=L+1:n*Nt

		   for(int i=L;i<n*Nt;i++) {
//             for kk=1:Lf
			   for(int kk=0;kk<Lf;kk++) {
//                 if abs(err(kk,i))>0.1
				   if(err.mat[kk][i].mag()>0.1) {
//                     num_error(da,num_packet)=num_error(da,num_packet)+1;
					   num_error.mat[snr_sweep][num_packet]++;
//                 end;
				   }
//             end;
			   }
//         end;
		   }
//         TEB(da,num_packet)=num_error(da,num_packet)/((n*Nt-L)*Lf);

		   TEB.mat[snr_sweep][num_packet]=num_error.mat[snr_sweep][num_packet]/(double)((n*Nt-L)*Lf);
//     end;
		}
//     SER(da) = sum(TEB(da,:),2)/Nbpacket;
// end;
	}
	cout << "TEB / SNR " << endl;
        DMatrix dtemp(TEB.linesum()/(double)Nbpacket);
	cout << TEB << "    :    " << dtemp;
	return TEB.linesum()/(double)Nbpacket;
// grid,x = 1.:1.:6.;
// grid on ;
// grid,semilogy(x,SER(:),'b');
// grid on ;
// grid,title('SER versus Es/N0');
//
// return;
}
// 
// 
// 
// 








