#include "tools/all.h"
#include "modulation/psk.h"
#include "rand/randgen.h"
#include "fading/fadingtdma.h"
#include "tools/clapack_tools.h"
#include "tools/ext_matrix.h"
#include "coding/block/ldpc.h"
#include "stbc/conversion.h"

//function [BER,BER1,BER2,BER3,BER4,TEB,TEB1,TEB2,TEB3,TEB4]= mimowangEM
int mimowangEM(UVector<DMatrix> &BER, UVector<DMatrix> &TEB) {
	int BER_TEB_CNT=5;
//% On considère une modulation MDP-4 transmise en OFDM %
//% sur deux antennes d'émission en utilisant un codage LDPC %
//% conforme à l'article de Wang (IEEE Trans. Commun. Jan. 2002 %
	int ncheck=512, block_size=1024, taille=block_size-ncheck, L=10;
//ncheck=512;
//block_size=1024;
//taille=block_size-ncheck;
//L=10;
//% taille de la séquence d'apprentissage du canal: soit L*Lf symboles MDP-4 %
	int M=2, Nbpacket=50, Lf=8; 
//M=2; 
//Nbpacket=50;
//Lf=8;
//% taille de la FFT: Lf doit être un multiple du nombre d'antennes d'émission %
//Nt=2;
	int Nt=2;
//% nombre d'antennes de transmission %
//Lh=4;
	int Lh=4;
//% étalement des trajets multiples maximum %
//fd=100;
	int fd=100;
//% frequence Doppler maximale dans le specte fréquentiel du canal %
	double coding_rate =0.5;
//coding_rate=0.5;
	
//% code LDPC 512/1024 %
//nb_step=5;
	int nb_step=5;
//% nombre d'itérations pour l'algorithme SAGE d'apprentissage du canal %
//EM_iter=3;
	int EM_iter=3;
//% nombre d'itérations pour l'algorithme EM de démodulation %

//s1=exp(j*pi/4);
//s2=exp(j*(pi/4+pi/2));
//s3=exp(j*(pi/4+pi));
//s4=exp(j*(pi/4+3*pi/2));
//% construction de la matrice de test pour le maximum de vraisemblance %
//MV(1,1)=s3; MV(1,2)=s3; 
//MV(2,1)=s3; MV(2,2)=s2; 
//MV(3,1)=s3; MV(3,2)=s4; 
//MV(4,1)=s3; MV(4,2)=s1; 
//MV(5,1)=s2; MV(5,2)=s3; 
//MV(6,1)=s2; MV(6,2)=s2; 
//MV(7,1)=s2; MV(7,2)=s4; 
//MV(8,1)=s2; MV(8,2)=s1; 
//MV(9,1)=s4; MV(9,2)=s3; 
//MV(10,1)=s4; MV(10,2)=s2; 
//MV(11,1)=s4; MV(11,2)=s4; 
//MV(12,1)=s4; MV(12,2)=s1; 
//MV(13,1)=s1; MV(13,2)=s3; 
//MV(14,1)=s1; MV(14,2)=s2; 
//MV(15,1)=s1; MV(15,2)=s4; 
//MV(16,1)=s1; MV(16,2)=s1; 
//
		dRandUniStatePtr tapseed;

		dRandUniInit(&tapseed, 557578, 0,1);

		int BLOCKSIZE_M_Nt_Lf = block_size/(M*Nt*Lf);



	//% on construit la matrice des coefficients de IFFT %

	ZMatrix ww(Lf, Lh);
//for k=1:Lf
	double isqrt_lf=1.0/sqrt(Lf), UTILITIS_M_2_PI_Lf=UTILITIS_M_2_PI/(double)Lf;
	for(int k=0;k<Lf;k++) {
//    for i=1:Lh
		for(int i=0;i<Lh;i++) {
//        ww(k,i)=(1/sqrt(Lf))*exp(-j*(2*pi*(k-1)*(i-1))/Lf);
			ww.mat[k][i]=isqrt_lf*DCplx::exp2RI(1, UTILITIS_M_2_PI_Lf*k*i);
//    end;
		}
//end;
	}
//
//Gn=[];
	ZMatrix Gn(Lf*L, Lf*2);
//W=[ww,zeros(Lf,Lh);zeros(Lf,Lh),ww];   
	ZMatrix W;
	RepeatMatrixDiagonally(ww, W, 2, DCplx(0,0));
//
//SNR=[12,13,14,15];
	int tsnr[]={12,13,14,15};
	WVector SNR(tsnr, 5);
	BER = UVector<DMatrix>(SNR.taille);
	TEB = UVector<DMatrix>(SNR.taille);

//for xx=1:length(SNR)
//    for jj=1:Nbpacket
//        BER(xx,jj)=0;
//        BER1(xx,jj)=0;
//        BER2(xx,jj)=0;
//        BER3(xx,jj)=0;
//        BER4(xx,jj)=0;
//    end;
//end;
//for dd=1:length(SNR)
//    TEB(dd,1)=0;
//    TEB1(dd,1)=0;
//    TEB2(dd,1)=0;
//    TEB3(dd,1)=0;
//    TEB4(dd,1)=0;
//end;
	DMatrix dtemp = DMatrix(SNR.taille, Nbpacket,0), dtemp2= DMatrix(SNR.taille,1);
	for(int xx=0;xx<BER_TEB_CNT;xx++) {
		BER.vect[xx]= dtemp;
		TEB.vect[xx]= dtemp2;
	}
	dtemp.~DMatrix();
	dtemp2.~DMatrix();

//
//H=ldpc_generate(ncheck,block_size,3,2,456);
	WLDPCParityCheck ldpc(ncheck, block_size);
//[H,G]=Ldpc_h2g(H);
	ldpc.make_dense_mixed();
//
//for da=1:length(SNR)
	TRandBitStatePtr seed;
	for(int da=0;da<SNR.taille;da++) {
//    for num_packet=1:Nbpacket
		UVector<WVector> seq(Nbpacket);
		WVector gseq;

		for(int num_packet=0;num_packet<Nbpacket;num_packet++) {
//        % On code une séquence aléatoire à l'aide d'un code LDPC de rendement 1/2 %
//        % 512/1024 %
//        seq(:,num_packet)=((sign(randn(1,size(G,1)))+1)/2)';
			gseq = seq.vect[num_packet]=WVector(block_size-ncheck);
			bRandBit1(&seed, seq.vect[num_packet]);
//        mseq=seq(:,num_packet);
//        gseq=mseq';
//        y=mod(mseq'*G,2);  % coding 
			WVector y=ldpc.encode(seq.vect[num_packet]);
//        
//        for i=1:(block_size/M)
//            ssn=y(M*(i-1)+1:M*i);
//            if(ssn==[1 1])
//                vn(i)=s1;
//            elseif(ssn==[0 1])
//                vn(i)=s2;
//            elseif(ssn==[0 0])
//                vn(i)=s3;
//            else
//                vn(i)=s4;
//            end;
//        end;
			ZVector vn = zPSK4Gray_mod(y);
//        
//        % création du codage spatio-temporel % 
//        for jj=1:(block_size/(M*Nt))
//            ss1(jj)=vn(Nt*(jj-1)+1);
//            ss2(jj)=vn(Nt*jj);
//        end;
			UVector<ZVector> ss = SplitMatrixByCol(G2(vn, vn.taille));
//        
//        % creation de la modulation multiporteuse OFDM %
//        for i=1:(BLOCKSIZE_M_Nt_Lf)
//            w1(:,i)=conj((ss1(Lf*(i-1)+1:Lf*i))');
//            w2(:,i)=conj((ss2(Lf*(i-1)+1:Lf*i))');
//            dcn(:,i)=conj((vn(Nt*Lf*(i-1)+1:Nt*Lf*i))');
//        end;
			ZMatrix w1 = ColumnPackedVectorAsMatrix(ss.vect[0], Lf, BLOCKSIZE_M_Nt_Lf);
			ZMatrix w2 = ColumnPackedVectorAsMatrix(ss.vect[1], Lf, BLOCKSIZE_M_Nt_Lf);
			ZMatrix dcn = ColumnPackedVectorAsMatrix(vn, Nt*Lf, BLOCKSIZE_M_Nt_Lf);

//        % création de la séquence d'apprentissage : à cause du format des %
//        % paquets (utilisation d'un code LDPC), on choisit de la rajouter %
//        % en plus des données %
//        
//        feet=L*M*Lf*Nt;
			int feet=L*M*Lf*Nt;
//        train=(rand(1,feet)>0.5)';
			WVector train(feet);
			bRandBit1(&seed, train);
//        for i=1:feet/M
//            dsn=train(M*(i-1)+1:M*i);
//            if (dsn==[1 1]')
//                vvn(i)=s1;
//            elseif(dsn==[0 1]')
//                vvn(i)=s2;
//            elseif(dsn==[0 0]')
//                vvn(i)=s3;
//            else
//                vvn(i)=s4;
//            end;
//        end;
			ZVector vvn=zPSK4Gray_mod(train);
//        
//        for jj=1:L*Lf
//            vvn1(jj)=vvn(Nt*(jj-1)+1);
//            vvn2(jj)=vvn(Nt*jj);
//        end;
			UVector<ZVector> svvn;
			StreamVector(vvn, svvn, 2);
//        
//        for i=1:L
//            training1(:,i)=conj((vvn1(Lf*(i-1)+1:Lf*i))');
//            training2(:,i)=conj((vvn2(Lf*(i-1)+1:Lf*i))');
//        end;
			UVector<ZMatrix> training(2);
			for(int i=0;i<2;i++) {
				training.vect[i]=ColumnPackedVectorAsMatrix(svvn.vect[i], Lf, L);
			}
//        
//        % On a à créer deux canaux de propagation sur chacun des deux flux %
//        % de sortie des deux antennes: chacun étant sélectif en fréquence: %
//        % le premier comporte deux coefficients et le deuxième trois %
//        
//        % première antenne %
//        z1=0.574*fadingtaps(fd,1,1,(BLOCKSIZE_M_Nt_Lf)+L,512);
//        z2=0.574*fadingtaps(fd,1,1,(BLOCKSIZE_M_Nt_Lf)+L,512);
//        z3=0.574*fadingtaps(fd,1,1,(BLOCKSIZE_M_Nt_Lf)+L,512);
//        z4=0.408*fadingtaps(fd,1,1,(BLOCKSIZE_M_Nt_Lf)+L,512);
//        %z5=0.408*fadingtaps(fd,1,1,(BLOCKSIZE_M_Nt_Lf)+L,512);
//        %z6=0.408*fadingtaps(fd,1,1,(BLOCKSIZE_M_Nt_Lf)+L,512);



			ZMatrix zfad = fadingtaps(fd, 1, 4, (BLOCKSIZE_M_Nt_Lf)+L, 512, &tapseed);
			double coef[]={0.574,0.574,0.574,0.408};
			
			MultEachColinMatrixWithY(zfad, DVector(coef,4));
//        
//        % deuxième antenne %
//        %z7=0.408*fadingtaps(fd,1,1,(BLOCKSIZE_M_Nt_Lf)+L,512);
//        %z8=0.408*fadingtaps(fd,1,1,(BLOCKSIZE_M_Nt_Lf)+L,512);
//        %z9=0.408*fadingtaps(fd,1,1,(BLOCKSIZE_M_Nt_Lf)+L,512);
//        %z10=0.707*fadingtaps(fd,1,1,(BLOCKSIZE_M_Nt_Lf)+L,512);
//        %z11=0.0707*fadingtaps(fd,1,1,(BLOCKSIZE_M_Nt_Lf)+L,512);
//        %z12=0.0.707*fadingtaps(fd,1,1,(BLOCKSIZE_M_Nt_Lf)+L,512);
//        
//        % on construit la matrice des coefficients de Rayleigh du canal % 
//        for time=1:L+(BLOCKSIZE_M_Nt_Lf) 
//            HH(:,time)=conj([z1(time),0,z2(time),0,z3(time),0,z4(time),0]');
//        end;
			ZMatrix HH= InterleaveTwoMatrixLineWise(zfad.transpose_nip(), ZMatrix(zfad.col, zfad.line));
//        HH
//        
//        % construction du signal émis %
//        for jj=1:(BLOCKSIZE_M_Nt_Lf) 
			UVector<ZMatrix> X(BLOCKSIZE_M_Nt_Lf);
			UVector<ZMatrix> d(BLOCKSIZE_M_Nt_Lf);
			DVector std(BLOCKSIZE_M_Nt_Lf);
			for(int jj=0;jj<BLOCKSIZE_M_Nt_Lf;jj++) {
//            X(:,jj)=[diag(w1(:,jj)),diag(w2(:,jj))]*W*HH(:,jj+L);
				X.vect[jj] = col(ColXinMatrixAsDiagonal(w1,jj), ColXinMatrixAsDiagonal(w2,jj))*W*SubMatrix(HH, 0, jj+L, HH.line, 1);
//            [d(:,jj),std(jj)]=awgn(coding_rate,Lf,SNR(da),X(:,jj));
				d.vect[jj]=awgn(coding_rate, Lf, SNR.vect[da], X.vect[jj], std.vect[jj]);
//        end;
			}

//        for kk=1:L
			UVector<ZMatrix> X0(L);
			UVector<ZMatrix> d0(L);
			DVector std3(L);
			for(int kk=0;kk<L;kk++) {
//            %X01(:,kk)=[diag(training1(:,kk)),zeros(Lf)]*W*HH(:,kk);
//            %[d01(:,kk),std1(kk)]=awgn(1,Lf,SNR(da),X01(:,kk));
//            %X02(:,kk)=[zeros(Lf),diag(training2(:,kk))]*W*HH(:,kk);
//            %[d02(:,kk),std2(kk)]=awgn(1,Lf,SNR(da),X02(:,kk));
//            X0(:,kk)=[diag(training1(:,kk)),diag(training2(:,kk))]*W*HH(:,kk);
				X0.vect[kk]=col(ColXinMatrixAsDiagonal(training.vect[0], kk), ColXinMatrixAsDiagonal(training.vect[1],kk))*W*SubMatrix(HH, 0,kk,HH.line,1);
//            [d0(:,kk),std3(kk)]=awgn(1,Lf,SNR(da),X0(:,kk));
				d0.vect[kk]=awgn(1,Lf, SNR.vect[da], X0.vect[kk], std3.vect[kk]);
//        end;
			}
//        
//        for i=1:L
			for(i=0;i<L;i++) {
//            Gn=[Gn;[diag(training1(:,i)),diag(training2(:,i))]*W];
				InsertSubMatrixIntoMatrix(Gn, col(ColXinMatrixAsDiagonal(training.vect[0],i), ColXinMatrixAsDiagonal(training.vect[1],i))*W,i*training.vect[0].line ,0);
//            rn((i-1)*Lf+1:i*Lf,1)=d0(:,i);
//        end;
			}
			ZMatrix rn=line(d0);
//        hest=pinv(Gn)*rn
			ZMatrix hest = pinv(Gn*rn);
//        Gn=[];
			Gn=ZMatrix(Lf*L,2*L);
//        
//        % algorithme SAGE %
//        % valeur initiale des coefficients %
//        
//        %for mm=1:L
//        % valeur initiale des coefficients %
//        %hest1sage(:,mm)=ww'*inv(diag(training1(:,mm)))*d01(:,mm);
//        %hest2sage(:,mm)=ww'*inv(diag(training2(:,mm)))*d02(:,mm);
//        %Z1(:,mm)=diag(training1(:,mm))*ww*hest1sage(:,mm);
//        %Z2(:,mm)=diag(training2(:,mm))*ww*hest2sage(:,mm);
//        
//        %for ii=1:nb_step 
//        %Y1(:,mm)=d0(:,mm)-Z2(:,mm);
//        %hest1sage(:,mm)=ww'*inv(diag(training1(:,mm)))*Y1(:,mm);
//        %K1(:,mm)=diag(training1(:,mm))*ww*hest1sage(:,mm);
//        %Y2(:,mm)=d0(:,mm)-Z1(:,mm);
//        %hest2sage(:,mm)=ww'*inv(diag(training2(:,mm)))*Y2(:,mm);
//        %K2(:,mm)=diag(training2(:,mm))*ww*hest2sage(:,mm);
//        %Z1(:,mm)=K1(:,mm);
//        %Z2(:,mm)=K2(:,mm);
//        %ii=ii+1;
//        %end;
//        %end;
//        %GG=hest1sage(:,1:L);
//        %KK=hest2sage(:,1:L);
//        %hest1sge=(1/L)*sum(GG,2);
//        %hest2sge=(1/L)*sum(KK,2);
//        %hest(1:Lh,1)=hest1sge;
//        %hest(Lh+1:Nt*Lh,1)=hest2sge;
//        %hest;

//        sigmah=diag(hest)*(diag(hest))';
			ZMatrix diag_hest=ColXinMatrixAsDiagonal(hest,0);
			ZMatrix sigmah = diag_hest*diag_hest.conj();
//        hestim=W*hest;
			ZMatrix hestim=W*hest;
//        
//        for ii=1:(BLOCKSIZE_M_Nt_Lf)
			//int BLOCKSIZE_M_Nt_Lf=(BLOCKSIZE_M_Nt_Lf;
			ZMatrix ddn(Nt*Lf, BLOCKSIZE_M_Nt_Lf);
			ZMatrix zx1(Lf, BLOCKSIZE_M_Nt_Lf);
			ZMatrix zx2(Lf, BLOCKSIZE_M_Nt_Lf);

			UVector<ZMatrix> htil(BLOCKSIZE_M_Nt_Lf);
			DMatrix delta1(Lf, BLOCKSIZE_M_Nt_Lf), delta2(Lf, BLOCKSIZE_M_Nt_Lf),delta3(Lf, BLOCKSIZE_M_Nt_Lf),delta4(Lf, BLOCKSIZE_M_Nt_Lf);
			DMatrix zeta(4*Lf, BLOCKSIZE_M_Nt_Lf);


			for(int ii=0;ii<BLOCKSIZE_M_Nt_Lf;ii++) {
//            % On commence par une détection MVP classique %
//            for n=1:Lf
				DMatrix Q(16,Lf);
				for(int n=0;n<Lf;n++) {
//                Rmini=10^5;
					double Rmini=TINY;
					WVector minpos(2,0);

//                for m=1:16
//                    Q(m,n)=abs(d(n,ii)-hestim(n)*MV(m,1)-hestim(n+Lf)*MV(m,2))^2;               
//                    if Q(m,n)<Rmini
//                        min=m;
//                        Rmini=Q(m,n);
//                    end;
//                end;
					for(int p1=0;p1<4;p1++) {
						for(int p2=0;p2<4;p2++) {
							Q.mat[p1*4+p2][n]=(d.vect[ii].mat[n][0]-hestim.mat[n][0]*psk4gray.sym[p1]-hestim.mat[n+Lf][0]*psk4gray.sym[p2]).mod2();
							if(Q.mat[p1*4+p2][n]<Rmini) {
								Rmini=Q.mat[p1*4+p2][n];
								minpos.vect[0]=p1;minpos.vect[1]=p2;
							}
						}
					}
//                ddn(Nt*(n-1)+1,ii)=MV(min,1);
//                zx1(n,ii)=ddn(Nt*(n-1)+1,ii);
					zx1.mat[n][ii]=ddn.mat[Nt*n][ii]=psk4gray.sym[minpos.vect[0]];
//                ddn(n*Nt,ii)=MV(min,2);
//                zx2(n,ii)=ddn(n*Nt,ii);
					zx2.mat[n][ii]=ddn.mat[Nt*(n+1)-1][ii]=psk4gray.sym[minpos.vect[1]];
//            end;
				}
//            T=[diag(zx1(:,ii)),diag(zx2(:,ii))];
				ZMatrix T = col(ColXinMatrixAsDiagonal(zx1,ii), ColXinMatrixAsDiagonal(zx2,ii));
//            % Calcul des matrices htilde et sigmatildeh %
//            Q=inv((1/std(ii))*W'*T'*T*W+inv(sigmah))*W'*T'*(1/std(ii));

#pragma message("possible confusion on Q")
				Q = (((1/std.vect[ii])*W.hermit()*T.hermit()*T*W+sigmah.inv_nip() ).inv_nip()*W.hermit()*T.hermit()*(1/std.vect[ii])).real();
//            htil(:,ii)=Q*d(:,ii);
				//InsertSubMatrixIntoMatrix(htil, Q*d.vect[ii], 0, ii);
				htil.vect[ii]=Q*d.vect[ii];
//            sigmatil=sigmah-Q*T*W*sigmah;
				ZMatrix sigmatil=sigmah-Q*T*W*sigmah;
//            dsigma=W*sigmatil*W';
				ZMatrix dsigma=W*sigmatil*W.hermit();
//            % Première itération de l'algorithme E.M. (cf calculs Wang)%
//            % On reconstitue la matrice des symboles décidés %
//            for ll=1:EM_iter
				ZMatrix wwdiag;
				RepeatMatrixDiagonally(ExtractLineXFromMatrix(ww,n), wwdiag, 2, DCplx(0,0));

				for(int ll=0;ll<EM_ITER;ll++) {
//                % On calcule les métriques de l'E.M. %
//                for n=1:Lf
					DMatrix R(16, Lf);
					for(int n=0;n<Lf;n++) {
//                    Rmini=10^5;
						double Rmini=TINY;
//                    for m=1:16
//                        R(m,n)=(1/(Lf*std(ii)))*(abs(d(n,ii)-[MV(m,1),MV(m,2)]* [ww(n,:),zeros(1,Lh);zeros(1,Lh),ww(n,:)]*htil(:,ii))^2 + ...
//                            real([MV(m,1),MV(m,2)]*[dsigma(n,n),dsigma(n+Lf,n);dsigma(n,n+Lf),dsigma(n+Lf,n+Lf)]*[MV(m,1),MV(m,2)]'))-log(1/16);
//                        if R(m,n)<Rmini
//                            min=m;
//                            Rmini=R(m,n);
//                        end;
//                    end;
						// [dsigma(n,n),dsigma(n+Lf,n);dsigma(n,n+Lf),dsigma(n+Lf,n+Lf)]
						ZMatrix dsigma_mat(2,2);
						dsigma_mat.mat[0][0] = dsigma.mat[n][n];
						dsigma_mat.mat[0][1] = dsigma.mat[n+Lf][n];
						dsigma_mat.mat[1][0] = dsigma.mat[n][n+Lf];
						dsigma_mat.mat[1][1] = dsigma.mat[n+Lf][n+Lf];
						ZMatrix xtest(2,1)
						for(int p1=0;p1<4;p1++) {
							for(int p2=0;p2<4;p2++) {
								ZMatrix MV(1,2), 
								MV.mat[0][0]=psk4gray.sym[p1];MV.mat[0][1]=psk4gray.sym[p2];
								R.mat[p1*4+p2][n]=(1.0/(Lf*std.vect[ii])) * ( mod2(d.vect[ii].mat[n]-MV * wwdiag*htil.vect[ii] ) + ( MV*dsigma_mat*MV.hermit() ).re - log(1.0/16.0);
								if(R.mat[p1*4+p2][n] < Rmini) {
									xmin=MV;
									Rmini = R.mat[p1*4+p2][n];
								}
							}
						}
//                    mini(n,ii)=min;
//                    ddn(Nt*(n-1)+1,ii)=MV(min,1);
//                    xx1(n,ii)=ddn(Nt*(n-1)+1,ii);
						xx1.mat[n][ii]=ddn.mat[Nt*n][ii]=xmin.mat[0][0];
//                    ddn(n*Nt,ii)=MV(min,2);
//                    xx2(n,ii)=ddn(n*Nt,ii); 
						xx2.mat[n][ii]=ddn.mat[Nt*(1+n)-1][ii]=xmin.mat[0][1];
//                end;
					}
//                T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
					T = col(ColXinMatrixAsDiagonal(xx1, ii), ColXinMatrixAsDiagonal(xx2, ii));
//                % Calcul des matrices htilde et sigmatildeh %
//                Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
					Q = (1.0/std.vect[ii])*W.hermit()*T.hermit()*T*W*pinv(sigmah)).inv_nip()*W.hermit()*T.hermit()*(1.0/std.vect[ii]);
//                htil(:,ii)=Q*d(:,ii);
					htil.vect[ii] = Q*d.vect[ii];
//                sigmatil=sigmah-Q*T*W*sigmah;
					sigmatil =sigmah-Q*T*W*sigmah;
//                dsigma=W*sigmatil*W';
					dsigma=W*sigmatil*W.hermit();
//                ll=ll+1;
//            end;
				}
//            % Fin de la boucle EM et calcul des LLR's %
//            for n=1:Lf

				for(n=0;n<Lf;n++) {
//                delta1(n,ii)=log(exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
//                    log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n)));
					delta1.mat[n][ii]=log(exp(-R.mat[8][n]) + exp(-R.mat[9][n])+exp(-R.mat[10][n])+exp(-R.mat[11][n])+exp(-R.mat[12][n])+exp(-R.mat[13][n])+exp(-R.mat[14][n])+exp(-R.mat[15][n]))-...
				                      log(exp(-R.mat[4][n])+exp(-R.mat[5][n])+exp(-R.mat[6][n])+exp(-R.mat[7][n])+exp(-R.mat[0][n])+exp(-R.mat[1][n])+exp(-R.mat[2][n])+exp(-R.mat[3][n]));
//                delta2(n,ii)=log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
//                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n)));
					delta2.mat[n][ii]=log(exp(-R.mat[4][n]) + exp(-R.mat[5][n])+exp(-R.mat[6][n])+exp(-R.mat[7][n])+exp(-R.mat[12][n])+exp(-R.mat[13][n])+exp(-R.mat[14][n])+exp(-R.mat[15][n]))-...
				                      log(exp(-R.mat[0][n])+exp(-R.mat[1][n])+exp(-R.mat[2][n])+exp(-R.mat[3][n])+exp(-R.mat[8][n])+exp(-R.mat[9][n])+exp(-R.mat[10][n])+exp(-R.mat[11][n]));
//                delta3(n,ii)=log(exp(-R(3,n))+exp(-R(4,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(15,n))+exp(-R(16,n)))-...
//                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(5,n))+exp(-R(6,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(13,n))+exp(-R(14,n)));
					delta3.mat[n][ii]=log(exp(-R.mat[2][n]) + exp(-R.mat[3][n])+exp(-R.mat[6][n])+exp(-R.mat[7][n])+exp(-R.mat[10][n])+exp(-R.mat[11][n])+exp(-R.mat[14][n])+exp(-R.mat[15][n]))-...
				                      log(exp(-R.mat[0][n])+exp(-R.mat[1][n])+exp(-R.mat[4][n])+exp(-R.mat[5][n])+exp(-R.mat[8][n])+exp(-R.mat[9][n])+exp(-R.mat[12][n])+exp(-R.mat[13][n]));
//                delta4(n,ii)=log(exp(-R(2,n))+exp(-R(4,n))+exp(-R(6,n))+exp(-R(8,n))+exp(-R(10,n))+exp(-R(12,n))+exp(-R(14,n))+exp(-R(16,n)))-...
//                    log(exp(-R(1,n))+exp(-R(3,n))+exp(-R(5,n))+exp(-R(7,n))+exp(-R(9,n))+exp(-R(11,n))+exp(-R(13,n))+exp(-R(15,n)));
					delta4.mat[n][ii]=log(exp(-R.mat[1][n]) + exp(-R.mat[3][n])+exp(-R.mat[5][n])+exp(-R.mat[7][n])+exp(-R.mat[9][n])+exp(-R.mat[11][n])+exp(-R.mat[13][n])+exp(-R.mat[15][n]))-...
				                      log(exp(-R.mat[0][n])+exp(-R.mat[2][n])+exp(-R.mat[4][n])+exp(-R.mat[6][n])+exp(-R.mat[8][n])+exp(-R.mat[10][n])+exp(-R.mat[12][n])+exp(-R.mat[14][n]));
//                zeta((n-1)*4+1,ii)=delta1(n,ii);
					zeta.mat[n*4][ii]=delta1.mat[n][ii];
//                zeta((n-1)*4+2,ii)=delta2(n,ii);
					zeta.mat[n*4+1][ii]=delta2.mat[n][ii];
//                zeta((n-1)*4+3,ii)=delta3(n,ii);
					zeta.mat[n*4+2][ii]=delta3.mat[n][ii];
//                zeta((n-1)*4+4,ii)=delta4(n,ii);
					zeta.mat[n*4+3][ii]=delta4.mat[n][ii];
//            end;
				}
//        end;  
			}
//        mini;
//        delta=reshape(zeta,block_size,1);
			DMatrix delta = ReshapeColumnwise(zeta, block_size,1);
//        % Calcul des probabilités P0 et P1 %
//        P0=1./(1+exp(delta));
			DMatrix  P0 = 1.0/(1.0+exp(delta));

//        for i=1:block_size
//            % On évite les problèmes de log 0 %
//            if P0(i)>=0.99999
//                P0(i)=0.99999;
//            elseif P0(i)<=0.00001
//                P0(i)=0.00001;
//            end;
//        end;
			Bound(P0, 0.00001, 0.99999);
//        P1=1-P0;
			DMatrix P1=1.0-P0;
			WVector z_hat = ldpc.decode(P0, APP);
//        [z_hat,success,k,APP]=ldpc_decode(P0,P1,H);
//        k;
//        success;
			// ldpc.prev_decode_status;
//        extrinsic=APP-delta;
			DVector extrinsic = APP-delta;
//        P0=1./(1+exp(extrinsic));
			P0=1.0/(1+exp(extrinsic));
//        for i=1:block_size
//            % On évite les problèmes de log 0 %
//            if P0(i)>=0.99999
//                P0(i)=0.99999;
//            elseif P0(i)<=0.00001
//                P0(i)=0.00001;
//            end;
//        end;
			Bound(P0, 0.00001, 0.99999);
//        P1=1-P0;
			P1=1-P0;
//        P1;
//        x_hat = z_hat(size(G,2)+1-size(G,1):size(G,2));
//        x_hat = x_hat';
			WVector x_hat = VectorInLineMatrixRepresentation(z_hat.copy(ldpc.GCol()-ldpc.GLine(), ldpc.GCol()));
			
//        o = xor(x_hat,gseq);
//        numerror = sum(o);
			int numerror=x_hat.diff_count(gseq);
//        BER(da,num_packet) = numerror/(taille); 
			BER.vect[0]=numerror/x_hat.taille;
//        % fin de la première turbo-itération: on commence la seconde itération turbo %
//        % on utilise les loglikelihood ratios au niveau des bits en sortie %
//        % du LDPC pour les convertir en rapports de vraisemblance au niveau %
//        % des symboles %
//        % On commence par calculer les nouvelles probabilités à priori sur les symboles %
//        % à partir des APP extrinsèques du LDPC %
//        
			WMatrix sP0 = ReshapeColumnwise(P0, M*Nt, block_size/(M*Nt)); 
			WMatrix sP1 = ReshapeColumnwise(P1, M*Nt, block_size/(M*Nt));
			DMatrix P(16,block_size/(M*Nt));
//        for i=1:(block_size/(M*Nt))
//            sP0(:,i)=P0(M*Nt*(i-1)+1:M*Nt*i,1);
//            sP1(:,i)=P1(M*Nt*(i-1)+1:M*Nt*i,1); 
//            P(1,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP0(4,i));
//            P(2,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP1(4,i));
//            P(3,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP0(4,i));
//            P(4,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP1(4,i));
//            P(5,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP0(4,i));
//            P(6,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP1(4,i));
//            P(7,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP0(4,i));
//            P(8,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP1(4,i));
//            P(9,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP0(4,i));
//            P(10,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP1(4,i));
//            P(11,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP0(4,i));
//            P(12,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP1(4,i));
//            P(13,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP0(4,i));
//            P(14,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP1(4,i));
//            P(15,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP0(4,i));
//            P(16,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP1(4,i));
//            sumP=sum(exp(P(:,i)),1);
//            P(1:16,i)=exp(P(1:16,i))/sumP;
//        end;
			double sumP=0.0;
			for(i=0;i<block_size/(double)(M*Nt);i++) {
				double sumP=0.0;
				sumP+=P.mat[0][i]=log(sP0.mat[0][i])+log(sP0.mat[1][i])+log(sP0.mat[2][i])+log(sP0.mat[3][i]);
				sumP+=P.mat[1][i]=log(sP0.mat[0][i])+log(sP0.mat[1][i])+log(sP0.mat[2][i])+log(sP1.mat[3][i]);
				sumP+=P.mat[2][i]=log(sP0.mat[0][i])+log(sP0.mat[1][i])+log(sP1.mat[2][i])+log(sP0.mat[3][i]);
				sumP+=P.mat[3][i]=log(sP0.mat[0][i])+log(sP0.mat[1][i])+log(sP1.mat[2][i])+log(sP1.mat[3][i]);
				sumP+=P.mat[4][i]=log(sP0.mat[0][i])+log(sP1.mat[1][i])+log(sP0.mat[2][i])+log(sP0.mat[3][i]);
				sumP+=P.mat[5][i]=log(sP0.mat[0][i])+log(sP1.mat[1][i])+log(sP0.mat[2][i])+log(sP1.mat[3][i]);
				sumP+=P.mat[6][i]=log(sP0.mat[0][i])+log(sP1.mat[1][i])+log(sP1.mat[2][i])+log(sP0.mat[3][i]);
				sumP+=P.mat[7][i]=log(sP0.mat[0][i])+log(sP1.mat[1][i])+log(sP1.mat[2][i])+log(sP1.mat[3][i]);
				sumP+=P.mat[8][i]=log(sP1.mat[0][i])+log(sP0.mat[1][i])+log(sP0.mat[2][i])+log(sP0.mat[3][i]);
				sumP+=P.mat[9][i]=log(sP1.mat[0][i])+log(sP0.mat[1][i])+log(sP0.mat[2][i])+log(sP1.mat[3][i]);
				sumP+=P.mat[10][i]=log(sP1.mat[0][i])+log(sP0.mat[1][i])+log(sP1.mat[2][i])+log(sP0.mat[3][i]);
				sumP+=P.mat[11][i]=log(sP1.mat[0][i])+log(sP0.mat[1][i])+log(sP1.mat[2][i])+log(sP1.mat[3][i]);
				sumP+=P.mat[12][i]=log(sP1.mat[0][i])+log(sP1.mat[1][i])+log(sP0.mat[2][i])+log(sP0.mat[3][i]);
				sumP+=P.mat[13][i]=log(sP1.mat[0][i])+log(sP1.mat[1][i])+log(sP0.mat[2][i])+log(sP1.mat[3][i]);
				sumP+=P.mat[14][i]=log(sP1.mat[0][i])+log(sP1.mat[1][i])+log(sP1.mat[2][i])+log(sP0.mat[3][i]);
				sumP+=P.mat[15][i]=log(sP1.mat[0][i])+log(sP1.mat[1][i])+log(sP1.mat[2][i])+log(sP1.mat[3][i]);
				MultColXinMatrixWithY(P,i, 1.0/sumP);
//            sumP=sum(exp(P(:,i)),1);
//            P(1:16,i)=exp(P(1:16,i))/sumP;
			}
//        
//        PP=reshape(P,16*Lf,BLOCKSIZE_M_Nt_Lf);
			DMatrix PP=ReshapeColumnwise(P, 16*Lf, BLOCKSIZE_M_Nt_Lf);
//        PQ=reshape(PP,16,Lf*BLOCKSIZE_M_Nt_Lf);
			int BLOCKSIZE_M_Nt = block_size/(M*Nt);
			DMatrix PQ=ReshapeColumnwise(PP, 16, BLOCKSIZE_M_Nt);
//        for i=1:Lf*BLOCKSIZE_M_Nt_Lf
//            [ZX(i),m]=max(PQ(:,i));

//            N(i)=m;
//        end;
			WMatrix ZXPos;
			DMatrix ZX=ColMax(PQ, ZXPos);
//        ZX=reshape(N,Lf,BLOCKSIZE_M_Nt_Lf);
//#pragma warning("ZX reused before being used !")
			ZX=ReshapeColumnwise(N,Lf, BLOCKSIZE_M_Nt_Lf);
//        for ii=1:(BLOCKSIZE_M_Nt_Lf)
//            for nn=1:Lf
//                xx1(nn,ii)=MV(ZX(nn,ii),1);
//                xx2(nn,ii)=MV(ZX(nn,ii),2);
//            end;
//        end;

#pragma warning("complicate conversion")
			for(ii=0;ii<BLOCKSIZE_M_Nt_Lf;i++) {
				for(int nn<0;nn<Lf;nn++) {
					xx1.mat[nn][ii]=psk8gray.sym[n];
					xx2.mat[nn][ii]=psk8gray.sym[n];
				}
			}
//        
//        % On recommence l'algorithme EM %
//        for ii=1:(BLOCKSIZE_M_Nt_Lf)
//            T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
//            % Calcul des matrices htilde et sigmatildeh %
//            Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
//            sigmatil=sigmah-Q*T*W*sigmah;
//            dsigma=W*sigmatil*W';
//            for ll=1:EM_iter
//                % On calcule les métriques de l'E.M. %
//                for n=1:Lf
//                    Rmini=10^5;
//                    for m=1:16
//                        R(m,n)=(1/(Lf*std(ii)))*(abs(d(n,ii)-[MV(m,1),MV(m,2)]* [ww(n,:),zeros(1,Lh);zeros(1,Lh),ww(n,:)]*htil(:,ii))^2 + ...
//                            real([MV(m,1),MV(m,2)]*[dsigma(n,n),dsigma(n+Lf,n);dsigma(n,n+Lf),dsigma(n+Lf,n+Lf)]*[MV(m,1),MV(m,2)]'))-log(PP(m+16*(n-1),ii));
//                        if R(m,n)<Rmini
//                            min=m;
//                            Rmini=R(m,n);
//                        end;
//                    end;
//                    ddn(Nt*(n-1)+1,ii)=MV(min,1);
//                    xx1(n,ii)=ddn(Nt*(n-1)+1,ii);
//                    ddn(n*Nt,ii)=MV(min,2);
//                    xx2(n,ii)=ddn(n*Nt,ii); 
//                end;
//                T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
//                % Calcul des matrices htilde et sigmatildeh %
//                Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
//                htil(:,ii)=Q*d(:,ii);
//                sigmatil=sigmah-Q*T*W*sigmah;
//                dsigma=W*sigmatil*W';
//                ll=ll+1;
//            end;
//            % Fin de la boucle EM, on recalcule les rapports de vraisemblance %
//            for n=1:Lf
//                delta1(n,ii)=log(exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
//                    log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n)));
//                delta2(n,ii)=log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
//                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n)));
//                delta3(n,ii)=log(exp(-R(3,n))+exp(-R(4,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(15,n))+exp(-R(16,n)))-...
//                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(5,n))+exp(-R(6,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(13,n))+exp(-R(14,n)));
//                delta4(n,ii)=log(exp(-R(2,n))+exp(-R(4,n))+exp(-R(6,n))+exp(-R(8,n))+exp(-R(10,n))+exp(-R(12,n))+exp(-R(14,n))+exp(-R(16,n)))-...
//                    log(exp(-R(1,n))+exp(-R(3,n))+exp(-R(5,n))+exp(-R(7,n))+exp(-R(9,n))+exp(-R(11,n))+exp(-R(13,n))+exp(-R(15,n)));
//                zeta((n-1)*4+1,ii)=delta1(n,ii);
//                zeta((n-1)*4+2,ii)=delta2(n,ii);
//                zeta((n-1)*4+3,ii)=delta3(n,ii);
//                zeta((n-1)*4+4,ii)=delta4(n,ii);
//            end;
//        end;

			for(int ii=0;ii<BLOCKSIZE_M_Nt_Lf;ii++) {
//            % On commence par une détection MVP classique %
				DMatrix Q(16,Lf);
				for(int n=0;n<Lf;n++) {
					double Rmini=TINY;
					WVector minpos(2,0);
					for(int p1=0;p1<4;p1++) {
						for(int p2=0;p2<4;p2++) {
							Q.mat[p1*4+p2][n]=(d.vect[ii].mat[n,0]-hestim.mat[n][0]*psk4gray.sym[p1]-hestim.mat[n+Lf].[0]*psk4gray.sym[p2]).mod2();
							if(Q[p1*4+p2][n]<Rmini) {
								Rmini=Q[p1*4+p2][n];
								minpos.vect[0]=p1;minpos.vect[0]=p2;
							}
						}
					}
					zx1.mat[n][ii]=ddn.mat[Nt*n][ii]=psk4gray.sym[minpos.vect[p1]];
					zx2.mat[n][ii]=ddn.mat[Nt*(n+1)-1][ii]=psk4gray.sym[minpos.vect[p2]];
				}
				ZMatrix T = col(ColXinMatrixAsDiagonal(zx1,ii), ColXinMatrixAsDiagonal(zx2,ii));
//            % Calcul des matrices htilde et sigmatildeh %
				Q = ((1/std.vect[ii])*W.hermit()*T.hermit()*T*W+inv(sigmah))*W.hermit()*T.hermit()*(1/std.vect[ii])).inv();
				htil.vect[ii]=Q*d.vect[ii];
				ZMatrix sigmatil=sigmah-Q*T*W*sigmh;
				ZMatrix dsigma=W*sigmatil*W.hermit();
//            % Première itération de l'algorithme E.M. (cf calculs Wang)%
//            % On reconstitue la matrice des symboles décidés %
				ZMatrix wwdiag;
				RepeatMatrixDiagonally(ExtractLineXFromMatrix(ww,n), wwdiag, 2, DCplx(0,0));
				for(int ll=0;ll<EM_ITER;ll++) {
//                % On calcule les métriques de l'E.M. %
					DMatrix R(16, Lf);
					for(int n=0;n<Lf;n++) {
						double Rmini=TINY;
						ZMatrix dsigma_mat(2,2);
						dsigma_mat.mat[0][0] = dsigma.mat[n][n];
						dsigma_mat.mat[0][1] = dsigma.mat[n+Lf][n];
						dsigma_mat.mat[1][0] = dsigma.mat[n][n+Lf];
						dsigma_mat.mat[1][1] = dsigma.mat[n+Lf][n+Lf];
						ZMatrix xtest(2,1)
						for(int p1=0;p1<4;p1++) {
							for(int p2=0;p2<4;p2++) {
								ZMatrix MV(1,2), 
								MV.mat[0][0]=psk4gray.sym[p1];MV.mat[0][1]=psk4gray.sym[p2];
								R.mat[p1*4+p2][n]=(1.0/(Lf*std.vect[ii])) * ( mod2(d.vect[ii].mat[n]-MV * wwdiag*htil.vect[ii] ) + ( MV*dsigma_mat*MV.hermit() ).re - log(1.0/16.0);
								if(R.mat[p1*4+p2][n] < Rmini) {
									xmin=MV;
									Rmini = R.mat[p1*4+p2][n];
								}
							}
						}
						xx1.mat[n][ii]=ddn.mat[Nt*n][ii]=xmin.mat[0][0];
						xx2.mat[n][ii]=ddn.mat[Nt*(1+n)-1][ii]=xmin.mat[0][1];
					}
					T = col(ColXinMatrixAsDiagonal(xx1, ii), ColXinMatrixAsDiagonal(xx2, ii));
//                % Calcul des matrices htilde et sigmatildeh %
					Q = (1.0/std.vect[ii])*W.hermit()*T.hermit()*T*W*pinv(sigmah)).inv_nip()*W.hermit()*T.hermit()*(1.0/std.vect[ii]);
					htil.vect[ii] = Q*d.vect[ii];
					sigmatil =sigmah-Q*T*W*sigmah;
					dsigma=W*sigmatil*W.hermit();
				}
//            % Fin de la boucle EM et calcul des LLR's %
				for(n=0;n<Lf;n++) {
					delta1.mat[n][ii]=log(exp(-R.mat[8][n]) + exp(-R.mat[9][n])+exp(-R.mat[10][n])+exp(-R.mat[11][n])+exp(-R.mat[12][n])+exp(-R.mat[13][n])+exp(-R.mat[14][n])+exp(-R.mat[15][n]))-...
				                      log(exp(-R.mat[4][n])+exp(-R.mat[5][n])+exp(-R.mat[6][n])+exp(-R.mat[7][n])+exp(-R.mat[0][n])+exp(-R.mat[1][n])+exp(-R.mat[2][n])+exp(-R.mat[3][n]));
					delta2.mat[n][ii]=log(exp(-R.mat[4][n]) + exp(-R.mat[5][n])+exp(-R.mat[6][n])+exp(-R.mat[7][n])+exp(-R.mat[12][n])+exp(-R.mat[13][n])+exp(-R.mat[14][n])+exp(-R.mat[15][n]))-...
				                      log(exp(-R.mat[0][n])+exp(-R.mat[1][n])+exp(-R.mat[2][n])+exp(-R.mat[3][n])+exp(-R.mat[8][n])+exp(-R.mat[9][n])+exp(-R.mat[10][n])+exp(-R.mat[11][n]));
					delta3.mat[n][ii]=log(exp(-R.mat[2][n]) + exp(-R.mat[3][n])+exp(-R.mat[6][n])+exp(-R.mat[7][n])+exp(-R.mat[10][n])+exp(-R.mat[11][n])+exp(-R.mat[14][n])+exp(-R.mat[15][n]))-...
				                      log(exp(-R.mat[0][n])+exp(-R.mat[1][n])+exp(-R.mat[4][n])+exp(-R.mat[5][n])+exp(-R.mat[8][n])+exp(-R.mat[9][n])+exp(-R.mat[12][n])+exp(-R.mat[13][n]));
					delta4.mat[n][ii]=log(exp(-R.mat[1][n]) + exp(-R.mat[3][n])+exp(-R.mat[5][n])+exp(-R.mat[7][n])+exp(-R.mat[9][n])+exp(-R.mat[11][n])+exp(-R.mat[13][n])+exp(-R.mat[15][n]))-...
				                      log(exp(-R.mat[0][n])+exp(-R.mat[2][n])+exp(-R.mat[4][n])+exp(-R.mat[6][n])+exp(-R.mat[8][n])+exp(-R.mat[10][n])+exp(-R.mat[12][n])+exp(-R.mat[14][n]));
					zeta.mat[n*4][ii]=delta1.mat[n][ii];
					zeta.mat[n*4+1][ii]=delta2.mat[n][ii];
					zeta.mat[n*4+2][ii]=delta3.mat[n][ii];
					zeta.mat[n*4+3][ii]=delta4.mat[n][ii];
				}
			}

#pragma message("em2 done ...")
//        delta=reshape(zeta,block_size,1);
//        % Calcul des informations extrinsèques %
//        delta=delta-extrinsic;
//        % Calcul des probabilités P0 et P1 %
//        P0=1./(1+exp(delta));
//        for i=1:block_size
//            % On évite les problèmes de log 0 %
//            if P0(i)>=0.99999
//                P0(i)=0.99999;
//            elseif P0(i)<=0.00001
//                P0(i)=0.00001;
//            end;
//        end;
//        P1=1-P0;
//        P1;
//        [z_hat,success,k,APP]=ldpc_decode(P0,P1,H);
//        success;
//        k;
//        extrinsic=APP-delta;
//        P0=1./(1+exp(extrinsic));
//        for i=1:block_size
//            % On évite les problèmes de log 0 %
//            if P0(i)>=0.99999
//                P0(i)=0.99999;
//            elseif P0(i)<=0.00001
//                P0(i)=0.00001;
//            end;
//        end;
//        P1=1-P0;
//        x_hat = z_hat(size(G,2)+1-size(G,1):size(G,2));
//        x_hat = x_hat';
//        o = xor(x_hat,gseq);
//        numerror1 = sum(o);
//        BER1(da,num_packet) = numerror1/(taille); 
//        % fin de la seconde turbo-itération: on commence la troisième itération turbo %
//        % on utilise les loglikelihood ratios au niveau des bits en sortie %
//        % du LDPC pour les convertir en rapports de vraisemblance au niveau %
//        % des symboles %
//        % On commence par calculer les nouvelles probabilités à priori sur les symboles %
//        % à partir des APP extrinsèques du LDPC %
//        
//        
//        for i=1:(block_size/(M*Nt))
//            sP0(:,i)=P0(M*Nt*(i-1)+1:M*Nt*i,1);
//            sP1(:,i)=P1(M*Nt*(i-1)+1:M*Nt*i,1); 
//            P(1,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP0(4,i));
//            P(2,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP1(4,i));
//            P(3,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP0(4,i));
//            P(4,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP1(4,i));
//            P(5,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP0(4,i));
//            P(6,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP1(4,i));
//            P(7,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP0(4,i));
//            P(8,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP1(4,i));
//            P(9,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP0(4,i));
//            P(10,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP1(4,i));
//            P(11,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP0(4,i));
//            P(12,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP1(4,i));
//            P(13,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP0(4,i));
//            P(14,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP1(4,i));
//            P(15,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP0(4,i));
//            P(16,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP1(4,i));
//            sumP=sum(exp(P(:,i)),1);
//            P(1:16,i)=exp(P(1:16,i))/sumP;
//        end;
//        PP=reshape(P,16*Lf,BLOCKSIZE_M_Nt_Lf);
//        PQ=reshape(PP,16,Lf*BLOCKSIZE_M_Nt_Lf);
//        for i=1:Lf*BLOCKSIZE_M_Nt_Lf
//            [ZX(i),m]=max(PQ(:,i));
//            N(i)=m;
//        end;
//        ZX=reshape(N,Lf,BLOCKSIZE_M_Nt_Lf);
//        for ii=1:(BLOCKSIZE_M_Nt_Lf)
//            for nn=1:Lf
//                xx1(nn,ii)=MV(ZX(nn,ii),1);
//                xx2(nn,ii)=MV(ZX(nn,ii),2);
//            end;
//        end;
//        
//        % On recommence l'algorithme EM %
//        for ii=1:(BLOCKSIZE_M_Nt_Lf)
//            T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
//            % Calcul des matrices htilde et sigmatildeh %
//            Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
//            sigmatil=sigmah-Q*T*W*sigmah;
//            dsigma=W*sigmatil*W';
//            for ll=1:EM_iter
//                % On calcule les métriques de l'E.M. %
//                for n=1:Lf
//                    Rmini=10^5;
//                    for m=1:16
//                        R(m,n)=(1/(Lf*std(ii)))*(abs(d(n,ii)-[MV(m,1),MV(m,2)]* [ww(n,:),zeros(1,Lh);zeros(1,Lh),ww(n,:)]*htil(:,ii))^2 + ...
//                            real([MV(m,1),MV(m,2)]*[dsigma(n,n),dsigma(n+Lf,n);dsigma(n,n+Lf),dsigma(n+Lf,n+Lf)]*[MV(m,1),MV(m,2)]'))-log(PP(m+16*(n-1),ii));
//                        if R(m,n)<Rmini
//                            min=m;
//                            Rmini=R(m,n);
//                        end;
//                    end;
//                    ddn(Nt*(n-1)+1,ii)=MV(min,1);
//                    xx1(n,ii)=ddn(Nt*(n-1)+1,ii);
//                    ddn(n*Nt,ii)=MV(min,2);
//                    xx2(n,ii)=ddn(n*Nt,ii); 
//                end;
//                T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
//                % Calcul des matrices htilde et sigmatildeh %
//                Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
//                htil(:,ii)=Q*d(:,ii);
//                sigmatil=sigmah-Q*T*W*sigmah;
//                dsigma=W*sigmatil*W';
//                ll=ll+1;
//            end;
//            % Calcul des LLR's %
//            for n=1:Lf
//                delta1(n,ii)=log(exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
//                    log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n)));
//                delta2(n,ii)=log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
//                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n)));
//                delta3(n,ii)=log(exp(-R(3,n))+exp(-R(4,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(15,n))+exp(-R(16,n)))-...
//                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(5,n))+exp(-R(6,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(13,n))+exp(-R(14,n)));
//                delta4(n,ii)=log(exp(-R(2,n))+exp(-R(4,n))+exp(-R(6,n))+exp(-R(8,n))+exp(-R(10,n))+exp(-R(12,n))+exp(-R(14,n))+exp(-R(16,n)))-...
//                    log(exp(-R(1,n))+exp(-R(3,n))+exp(-R(5,n))+exp(-R(7,n))+exp(-R(9,n))+exp(-R(11,n))+exp(-R(13,n))+exp(-R(15,n)));
//                zeta((n-1)*4+1,ii)=delta1(n,ii);
//                zeta((n-1)*4+2,ii)=delta2(n,ii);
//                zeta((n-1)*4+3,ii)=delta3(n,ii);
//                zeta((n-1)*4+4,ii)=delta4(n,ii);
//            end;
//        end;
//        delta=reshape(zeta,block_size,1);
//        % Calcul des informations extrinsèques %
//        delta=delta-extrinsic;
//        % Calcul des probabilités P0 et P1 %
//        P0=1./(1+exp(delta));
//        for i=1:block_size
//            % On évite les problèmes de log 0 %
//            if P0(i)>=0.99999
//                P0(i)=0.99999;
//            elseif P0(i)<=0.00001
//                P0(i)=0.00001;
//            end;
//        end;
//        P1=1-P0;
//        P1;
//        [z_hat,success,k,APP]=ldpc_decode(P0,P1,H);
//        extrinsic=APP-delta;
//        P0=1./(1+exp(extrinsic));
//        for i=1:block_size
//            % On évite les problèmes de log 0 %
//            if P0(i)>=0.99999
//                P0(i)=0.99999;
//            elseif P0(i)<=0.00001
//                P0(i)=0.00001;
//            end;
//        end;
//        P1=1-P0;
//        x_hat = z_hat(size(G,2)+1-size(G,1):size(G,2));
//        x_hat = x_hat';
//        o = xor(x_hat,gseq);
//        numerror2 = sum(o);
//        BER2(da,num_packet) = numerror2/(taille); 
//        
//        % fin de la troisième turbo-itération: on commence la quatrième itération turbo %
//        % on utilise les loglikelihood ratios au niveau des bits en sortie %
//        % du LDPC pour les convertir en rapports de vraisemblance au niveau %
//        % des symboles %
//        % On commence par calculer les nouvelles probabilités à priori sur les symboles %
//        % à partir des APP extrinsèques du LDPC %
//        
//        
//        for i=1:(block_size/(M*Nt))
//            sP0(:,i)=P0(M*Nt*(i-1)+1:M*Nt*i,1);
//            sP1(:,i)=P1(M*Nt*(i-1)+1:M*Nt*i,1); 
//            P(1,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP0(4,i));
//            P(2,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP1(4,i));
//            P(3,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP0(4,i));
//            P(4,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP1(4,i));
//            P(5,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP0(4,i));
//            P(6,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP1(4,i));
//            P(7,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP0(4,i));
//            P(8,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP1(4,i));
//            P(9,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP0(4,i));
//            P(10,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP1(4,i));
//            P(11,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP0(4,i));
//            P(12,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP1(4,i));
//            P(13,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP0(4,i));
//            P(14,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP1(4,i));
//            P(15,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP0(4,i));
//            P(16,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP1(4,i));
//            sumP=sum(exp(P(:,i)),1);
//            P(1:16,i)=exp(P(1:16,i))/sumP;
//        end;
//        PP=reshape(P,16*Lf,BLOCKSIZE_M_Nt_Lf);
//        PQ=reshape(PP,16,Lf*BLOCKSIZE_M_Nt_Lf);
//        for i=1:Lf*BLOCKSIZE_M_Nt_Lf
//            [ZX(i),m]=max(PQ(:,i));
//            N(i)=m;
//        end;
//        ZX=reshape(N,Lf,BLOCKSIZE_M_Nt_Lf);
//        for ii=1:(BLOCKSIZE_M_Nt_Lf)
//            for nn=1:Lf
//                xx1(nn,ii)=MV(ZX(nn,ii),1);
//                xx2(nn,ii)=MV(ZX(nn,ii),2);
//            end;
//        end;
//        
//        % On recommence l'algorithme EM %
//        for ii=1:(BLOCKSIZE_M_Nt_Lf)
//            T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
//            % Calcul des matrices htilde et sigmatildeh %
//            Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
//            sigmatil=sigmah-Q*T*W*sigmah;
//            dsigma=W*sigmatil*W';
//            for ll=1:EM_iter
//                % On calcule les métriques de l'E.M. %
//                for n=1:Lf
//                    Rmini=10^5;
//                    for m=1:16
//                        R(m,n)=(1/(Lf*std(ii)))*(abs(d(n,ii)-[MV(m,1),MV(m,2)]* [ww(n,:),zeros(1,Lh);zeros(1,Lh),ww(n,:)]*htil(:,ii))^2 + ...
//                            real([MV(m,1),MV(m,2)]*[dsigma(n,n),dsigma(n+Lf,n);dsigma(n,n+Lf),dsigma(n+Lf,n+Lf)]*[MV(m,1),MV(m,2)]'))-log(PP(m+16*(n-1),ii));
//                        if R(m,n)<Rmini
//                            min=m;
//                            Rmini=R(m,n);
//                        end;
//                    end;
//                    ddn(Nt*(n-1)+1,ii)=MV(min,1);
//                    xx1(n,ii)=ddn(Nt*(n-1)+1,ii);
//                    ddn(n*Nt,ii)=MV(min,2);
//                    xx2(n,ii)=ddn(n*Nt,ii); 
//                end;
//                T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
//                % Calcul des matrices htilde et sigmatildeh %
//                Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
//                htil(:,ii)=Q*d(:,ii);
//                sigmatil=sigmah-Q*T*W*sigmah;
//                dsigma=W*sigmatil*W';
//                ll=ll+1;
//            end;
//            % Calcul des LLR's %
//            for n=1:Lf
//                delta1(n,ii)=log(exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
//                    log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n)));
//                delta2(n,ii)=log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
//                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n)));
//                delta3(n,ii)=log(exp(-R(3,n))+exp(-R(4,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(15,n))+exp(-R(16,n)))-...
//                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(5,n))+exp(-R(6,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(13,n))+exp(-R(14,n)));
//                delta4(n,ii)=log(exp(-R(2,n))+exp(-R(4,n))+exp(-R(6,n))+exp(-R(8,n))+exp(-R(10,n))+exp(-R(12,n))+exp(-R(14,n))+exp(-R(16,n)))-...
//                    log(exp(-R(1,n))+exp(-R(3,n))+exp(-R(5,n))+exp(-R(7,n))+exp(-R(9,n))+exp(-R(11,n))+exp(-R(13,n))+exp(-R(15,n)));
//                zeta((n-1)*4+1,ii)=delta1(n,ii);
//                zeta((n-1)*4+2,ii)=delta2(n,ii);
//                zeta((n-1)*4+3,ii)=delta3(n,ii);
//                zeta((n-1)*4+4,ii)=delta4(n,ii);
//            end;
//        end;
//        delta=reshape(zeta,block_size,1);
//        % Calcul des informations extrinsèques %
//        delta=delta-extrinsic;
//        % Calcul des probabilités P0 et P1 %
//        P0=1./(1+exp(delta));
//        for i=1:block_size
//            % On évite les problèmes de log 0 %
//            if P0(i)>=0.99999
//                P0(i)=0.99999;
//            elseif P0(i)<=0.00001
//                P0(i)=0.00001;
//            end;
//        end;
//        P1=1-P0;
//        P1;
//        [z_hat,success,k,APP]=ldpc_decode(P0,P1,H);
//        extrinsic=APP-delta;
//        P0=1./(1+exp(extrinsic));
//        for i=1:block_size
//            % On évite les problèmes de log 0 %
//            if P0(i)>=0.99999
//                P0(i)=0.99999;
//            elseif P0(i)<=0.00001
//                P0(i)=0.00001;
//            end;
//        end;
//        P1=1-P0;
//        x_hat = z_hat(size(G,2)+1-size(G,1):size(G,2));
//        x_hat = x_hat';
//        o = xor(x_hat,gseq);
//        numerror3 = sum(o);
//        BER3(da,num_packet) = numerror3/(taille);     
//        
//        % fin de la quatrième turbo-itération: on commence la quatrième itération turbo %
//        % on utilise les loglikelihood ratios au niveau des bits en sortie %
//        % du LDPC pour les convertir en rapports de vraisemblance au niveau %
//        % des symboles %
//        % On commence par calculer les nouvelles probabilités à priori sur les symboles %
//        % à partir des APP extrinsèques du LDPC %
//        
//        
//        for i=1:(block_size/(M*Nt))
//            sP0(:,i)=P0(M*Nt*(i-1)+1:M*Nt*i,1);
//            sP1(:,i)=P1(M*Nt*(i-1)+1:M*Nt*i,1); 
//            P(1,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP0(4,i));
//            P(2,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP1(4,i));
//            P(3,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP0(4,i));
//            P(4,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP1(4,i));
//            P(5,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP0(4,i));
//            P(6,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP1(4,i));
//            P(7,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP0(4,i));
//            P(8,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP1(4,i));
//            P(9,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP0(4,i));
//            P(10,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP1(4,i));
//            P(11,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP0(4,i));
//            P(12,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP1(4,i));
//            P(13,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP0(4,i));
//            P(14,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP1(4,i));
//            P(15,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP0(4,i));
//            P(16,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP1(4,i));
//            sumP=sum(exp(P(:,i)),1);
//            P(1:16,i)=exp(P(1:16,i))/sumP;
//        end;
//        PP=reshape(P,16*Lf,BLOCKSIZE_M_Nt_Lf);
//        PQ=reshape(PP,16,Lf*BLOCKSIZE_M_Nt_Lf);
//        for i=1:Lf*BLOCKSIZE_M_Nt_Lf
//            [ZX(i),m]=max(PQ(:,i));
//            N(i)=m;
//        end;
//        ZX=reshape(N,Lf,BLOCKSIZE_M_Nt_Lf);
//        for ii=1:(BLOCKSIZE_M_Nt_Lf)
//            for nn=1:Lf
//                xx1(nn,ii)=MV(ZX(nn,ii),1);
//                xx2(nn,ii)=MV(ZX(nn,ii),2);
//            end;
//        end;
//        
//        % On recommence l'algorithme EM %
//        for ii=1:(BLOCKSIZE_M_Nt_Lf)
//            T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
//            % Calcul des matrices htilde et sigmatildeh %
//            Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
//            sigmatil=sigmah-Q*T*W*sigmah;
//            dsigma=W*sigmatil*W';
//            for ll=1:EM_iter
//                % On calcule les métriques de l'E.M. %
//                for n=1:Lf
//                    Rmini=10^5;
//                    for m=1:16
//                        R(m,n)=(1/(Lf*std(ii)))*(abs(d(n,ii)-[MV(m,1),MV(m,2)]* [ww(n,:),zeros(1,Lh);zeros(1,Lh),ww(n,:)]*htil(:,ii))^2 + ...
//                            real([MV(m,1),MV(m,2)]*[dsigma(n,n),dsigma(n+Lf,n);dsigma(n,n+Lf),dsigma(n+Lf,n+Lf)]*[MV(m,1),MV(m,2)]'))-log(PP(m+16*(n-1),ii));
//                        if R(m,n)<Rmini
//                            min=m;
//                            Rmini=R(m,n);
//                        end;
//                    end;
//                    ddn(Nt*(n-1)+1,ii)=MV(min,1);
//                    xx1(n,ii)=ddn(Nt*(n-1)+1,ii);
//                    ddn(n*Nt,ii)=MV(min,2);
//                    xx2(n,ii)=ddn(n*Nt,ii); 
//                end;
//                T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
//                % Calcul des matrices htilde et sigmatildeh %
//                Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
//                htil(:,ii)=Q*d(:,ii);
//                sigmatil=sigmah-Q*T*W*sigmah;
//                dsigma=W*sigmatil*W';
//                ll=ll+1;
//            end;
//            % Calcul des LLR's %
//            for n=1:Lf
//                delta1(n,ii)=log(exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
//                    log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n)));
//                delta2(n,ii)=log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
//                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n)));
//                delta3(n,ii)=log(exp(-R(3,n))+exp(-R(4,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(15,n))+exp(-R(16,n)))-...
//                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(5,n))+exp(-R(6,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(13,n))+exp(-R(14,n)));
//                delta4(n,ii)=log(exp(-R(2,n))+exp(-R(4,n))+exp(-R(6,n))+exp(-R(8,n))+exp(-R(10,n))+exp(-R(12,n))+exp(-R(14,n))+exp(-R(16,n)))-...
//                    log(exp(-R(1,n))+exp(-R(3,n))+exp(-R(5,n))+exp(-R(7,n))+exp(-R(9,n))+exp(-R(11,n))+exp(-R(13,n))+exp(-R(15,n)));
//                zeta((n-1)*4+1,ii)=delta1(n,ii);
//                zeta((n-1)*4+2,ii)=delta2(n,ii);
//                zeta((n-1)*4+3,ii)=delta3(n,ii);
//                zeta((n-1)*4+4,ii)=delta4(n,ii);
//            end;
//        end;
//        delta=reshape(zeta,block_size,1);
//        % Calcul des informations extrinsèques %
//        delta=delta-extrinsic;
//        % Calcul des probabilités P0 et P1 %
//        P0=1./(1+exp(delta));
//        for i=1:block_size
//            % On évite les problèmes de log 0 %
//            if P0(i)>=0.99999
//                P0(i)=0.99999;
//            elseif P0(i)<=0.00001
//                P0(i)=0.00001;
//            end;
//        end;
//        P1=1-P0;
//        P1;
//        [z_hat,success,k,APP]=ldpc_decode(P0,P1,H);
//        extrinsic=APP-delta;
//        P0=1./(1+exp(extrinsic));
//        for i=1:block_size
//            % On évite les problèmes de log 0 %
//            if P0(i)>=0.99999
//                P0(i)=0.99999;
//            elseif P0(i)<=0.00001
//                P0(i)=0.00001;
//            end;
//        end;
//        P1=1-P0;
//        x_hat = z_hat(size(G,2)+1-size(G,1):size(G,2));
//        x_hat = x_hat';
//        o = xor(x_hat,gseq);
//        numerror4 = sum(o);
//        BER4(da,num_packet) = numerror4/(taille);    
//    end;  
//    TEB(da)=sum(BER(da,:),2)/Nbpacket;
//    TEB1(da)=sum(BER1(da,:),2)/Nbpacket;
//    TEB2(da)=sum(BER2(da,:),2)/Nbpacket;
//    TEB3(da)=sum(BER3(da,:),2)/Nbpacket;
//    TEB4(da)=sum(BER4(da,:),2)/Nbpacket;
//end;
//
//grid,x = 13:1.:14.;
//grid on ;
//grid,semilogy(x,TEB(:),'b', x, TEB1(:),'b*', x, TEB2(:),'r', x, TEB3(:),'go', x, TEB4(:),'b--');
//
//return;
//
//
//
//
//
//
//
//
//




