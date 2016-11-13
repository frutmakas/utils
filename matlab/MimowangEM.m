function [BER,BER1,BER2,BER3,BER4,TEB,TEB1,TEB2,TEB3,TEB4]= mimowangEM
% On considère une modulation MDP-4 transmise en OFDM %
% sur deux antennes d'émission en utilisant un codage LDPC %
% conforme à l'article de Wang (IEEE Trans. Commun. Jan. 2002 %
ncheck=512;
block_size=1024;
taille=block_size-ncheck;
L=10;
% taille de la séquence d'apprentissage du canal: soit L*Lf symboles MDP-4 %
M=2; 
Nbpacket=50;
Lf=8;
% taille de la FFT: Lf doit être un multiple du nombre d'antennes d'émission %
Nt=2;
% nombre d'antennes de transmission %
Lh=4;
% étalement des trajets multiples maximum %
fd=100;
% frequence Doppler maximale dans le specte fréquentiel du canal %
coding_rate=0.5;
% code LDPC 512/1024 %
nb_step=5;
% nombre d'itérations pour l'algorithme SAGE d'apprentissage du canal %
EM_iter=3;
% nombre d'itérations pour l'algorithme EM de démodulation %
s1=exp(j*pi/4);
s2=exp(j*(pi/4+pi/2));
s3=exp(j*(pi/4+pi));
s4=exp(j*(pi/4+3*pi/2));
% construction de la matrice de test pour le maximum de vraisemblance %
MV(1,1)=s3; MV(1,2)=s3; 
MV(2,1)=s3; MV(2,2)=s2; 
MV(3,1)=s3; MV(3,2)=s4; 
MV(4,1)=s3; MV(4,2)=s1; 
MV(5,1)=s2; MV(5,2)=s3; 
MV(6,1)=s2; MV(6,2)=s2; 
MV(7,1)=s2; MV(7,2)=s4; 
MV(8,1)=s2; MV(8,2)=s1; 
MV(9,1)=s4; MV(9,2)=s3; 
MV(10,1)=s4; MV(10,2)=s2; 
MV(11,1)=s4; MV(11,2)=s4; 
MV(12,1)=s4; MV(12,2)=s1; 
MV(13,1)=s1; MV(13,2)=s3; 
MV(14,1)=s1; MV(14,2)=s2; 
MV(15,1)=s1; MV(15,2)=s4; 
MV(16,1)=s1; MV(16,2)=s1; 

% on construit la matrice des coefficients de IFFT %
for k=1:Lf
    for i=1:Lh
        ww(k,i)=(1/sqrt(Lf))*exp(-j*(2*pi*(k-1)*(i-1))/Lf);
    end;
end;

Gn=[];
W=[ww,zeros(Lf,Lh);zeros(Lf,Lh),ww];   

SNR=[12,13,14,15];
for xx=1:length(SNR)
    for jj=1:Nbpacket
        BER(xx,jj)=0;
        BER1(xx,jj)=0;
        BER2(xx,jj)=0;
        BER3(xx,jj)=0;
        BER4(xx,jj)=0;
    end;
end;
for dd=1:length(SNR)
    TEB(dd,1)=0;
    TEB1(dd,1)=0;
    TEB2(dd,1)=0;
    TEB3(dd,1)=0;
    TEB4(dd,1)=0;
end;

H=ldpc_generate(ncheck,block_size,3,2,456);
[H,G]=Ldpc_h2g(H);

for da=1:length(SNR)
    for num_packet=1:Nbpacket
        % On code une séquence aléatoire à l'aide d'un code LDPC de rendement 1/2 %
        % 512/1024 %
        seq(:,num_packet)=((sign(randn(1,size(G,1)))+1)/2)';
        mseq=seq(:,num_packet);
        gseq=mseq';
        y=mod(mseq'*G,2);  % coding 
        
        for i=1:(block_size/M)
            ssn=y(M*(i-1)+1:M*i);
            if(ssn==[1 1])
                vn(i)=s1;
            elseif(ssn==[0 1])
                vn(i)=s2;
            elseif(ssn==[0 0])
                vn(i)=s3;
            else
                vn(i)=s4;
            end;
        end;
        
        % création du codage spatio-temporel % 
        for jj=1:(block_size/(M*Nt))
            ss1(jj)=vn(Nt*(jj-1)+1);
            ss2(jj)=vn(Nt*jj);
        end;
        
        % creation de la modulation multiporteuse OFDM %
        for i=1:(block_size/(M*Nt*Lf))
            w1(:,i)=conj((ss1(Lf*(i-1)+1:Lf*i))');
            w2(:,i)=conj((ss2(Lf*(i-1)+1:Lf*i))');
            dcn(:,i)=conj((vn(Nt*Lf*(i-1)+1:Nt*Lf*i))');
        end;
        % création de la séquence d'apprentissage : à cause du format des %
        % paquets (utilisation d'un code LDPC), on choisit de la rajouter %
        % en plus des données %
        
        feet=L*M*Lf*Nt;
        train=(rand(1,feet)>0.5)';
        for i=1:feet/M
            dsn=train(M*(i-1)+1:M*i);
            if (dsn==[1 1]')
                vvn(i)=s1;
            elseif(dsn==[0 1]')
                vvn(i)=s2;
            elseif(dsn==[0 0]')
                vvn(i)=s3;
            else
                vvn(i)=s4;
            end;
        end;
        
        for jj=1:L*Lf
            vvn1(jj)=vvn(Nt*(jj-1)+1);
            vvn2(jj)=vvn(Nt*jj);
        end;
        
        for i=1:L
            training1(:,i)=conj((vvn1(Lf*(i-1)+1:Lf*i))');
            training2(:,i)=conj((vvn2(Lf*(i-1)+1:Lf*i))');
        end;
        
        % On a à créer deux canaux de propagation sur chacun des deux flux %
        % de sortie des deux antennes: chacun étant sélectif en fréquence: %
        % le premier comporte deux coefficients et le deuxième trois %
        
        % première antenne %
        z1=0.574*fadingtaps(fd,1,1,(block_size/(M*Nt*Lf))+L,512);
        z2=0.574*fadingtaps(fd,1,1,(block_size/(M*Nt*Lf))+L,512);
        z3=0.574*fadingtaps(fd,1,1,(block_size/(M*Nt*Lf))+L,512);
        z4=0.408*fadingtaps(fd,1,1,(block_size/(M*Nt*Lf))+L,512);
        %z5=0.408*fadingtaps(fd,1,1,(block_size/(M*Nt*Lf))+L,512);
        %z6=0.408*fadingtaps(fd,1,1,(block_size/(M*Nt*Lf))+L,512);
        
        % deuxième antenne %
        %z7=0.408*fadingtaps(fd,1,1,(block_size/(M*Nt*Lf))+L,512);
        %z8=0.408*fadingtaps(fd,1,1,(block_size/(M*Nt*Lf))+L,512);
        %z9=0.408*fadingtaps(fd,1,1,(block_size/(M*Nt*Lf))+L,512);
        %z10=0.707*fadingtaps(fd,1,1,(block_size/(M*Nt*Lf))+L,512);
        %z11=0.0707*fadingtaps(fd,1,1,(block_size/(M*Nt*Lf))+L,512);
        %z12=0.0.707*fadingtaps(fd,1,1,(block_size/(M*Nt*Lf))+L,512);
        
        % on construit la matrice des coefficients de Rayleigh du canal % 
        for time=1:L+(block_size/(M*Nt*Lf)) 
            HH(:,time)=conj([z1(time),0,z2(time),0,z3(time),0,z4(time),0]');
        end;
        HH
        
        % construction du signal émis %
        for jj=1:(block_size/(M*Nt*Lf)) 
            X(:,jj)=[diag(w1(:,jj)),diag(w2(:,jj))]*W*HH(:,jj+L);
            [d(:,jj),std(jj)]=awgn(coding_rate,Lf,SNR(da),X(:,jj));
        end;
        for kk=1:L
            %X01(:,kk)=[diag(training1(:,kk)),zeros(Lf)]*W*HH(:,kk);
            %[d01(:,kk),std1(kk)]=awgn(1,Lf,SNR(da),X01(:,kk));
            %X02(:,kk)=[zeros(Lf),diag(training2(:,kk))]*W*HH(:,kk);
            %[d02(:,kk),std2(kk)]=awgn(1,Lf,SNR(da),X02(:,kk));
            X0(:,kk)=[diag(training1(:,kk)),diag(training2(:,kk))]*W*HH(:,kk);
            [d0(:,kk),std3(kk)]=awgn(1,Lf,SNR(da),X0(:,kk));
        end;
        
        for i=1:L
            Gn=[Gn;[diag(training1(:,i)),diag(training2(:,i))]*W];
            rn((i-1)*Lf+1:i*Lf,1)=d0(:,i);
        end;
        hest=pinv(Gn)*rn
        Gn=[];
        
        % algorithme SAGE %
        % valeur initiale des coefficients %
        
        %for mm=1:L
            % valeur initiale des coefficients %
            %hest1sage(:,mm)=ww'*inv(diag(training1(:,mm)))*d01(:,mm);
            %hest2sage(:,mm)=ww'*inv(diag(training2(:,mm)))*d02(:,mm);
            %Z1(:,mm)=diag(training1(:,mm))*ww*hest1sage(:,mm);
            %Z2(:,mm)=diag(training2(:,mm))*ww*hest2sage(:,mm);
            
            %for ii=1:nb_step 
                %Y1(:,mm)=d0(:,mm)-Z2(:,mm);
                %hest1sage(:,mm)=ww'*inv(diag(training1(:,mm)))*Y1(:,mm);
                %K1(:,mm)=diag(training1(:,mm))*ww*hest1sage(:,mm);
                %Y2(:,mm)=d0(:,mm)-Z1(:,mm);
                %hest2sage(:,mm)=ww'*inv(diag(training2(:,mm)))*Y2(:,mm);
                %K2(:,mm)=diag(training2(:,mm))*ww*hest2sage(:,mm);
                %Z1(:,mm)=K1(:,mm);
                %Z2(:,mm)=K2(:,mm);
                %ii=ii+1;
                %end;
                %end;
        %GG=hest1sage(:,1:L);
        %KK=hest2sage(:,1:L);
        %hest1sge=(1/L)*sum(GG,2);
        %hest2sge=(1/L)*sum(KK,2);
        %hest(1:Lh,1)=hest1sge;
        %hest(Lh+1:Nt*Lh,1)=hest2sge;
        %hest;
        sigmah=diag(hest)*(diag(hest))';
        hestim=W*hest;
        
        for ii=1:(block_size/(M*Nt*Lf))
            % On commence par une détection MVP classique %
            for n=1:Lf
                Rmini=10^5;
                for m=1:16
                    Q(m,n)=abs(d(n,ii)-hestim(n)*MV(m,1)-hestim(n+Lf)*MV(m,2))^2;               
                    if Q(m,n)<Rmini
                        min=m;
                        Rmini=Q(m,n);
                    end;
                end;
                ddn(Nt*(n-1)+1,ii)=MV(min,1);
                zx1(n,ii)=ddn(Nt*(n-1)+1,ii);
                ddn(n*Nt,ii)=MV(min,2);
                zx2(n,ii)=ddn(n*Nt,ii);
            end;
            T=[diag(zx1(:,ii)),diag(zx2(:,ii))];
            % Calcul des matrices htilde et sigmatildeh %
            Q=inv((1/std(ii))*W'*T'*T*W+inv(sigmah))*W'*T'*(1/std(ii));
            htil(:,ii)=Q*d(:,ii);
            sigmatil=sigmah-Q*T*W*sigmah;
            dsigma=W*sigmatil*W';
            % Première itération de l'algorithme E.M. (cf calculs Wang)%
            % On reconstitue la matrice des symboles décidés %
            for ll=1:EM_iter
                % On calcule les métriques de l'E.M. %
                for n=1:Lf
                    Rmini=10^5;
                    for m=1:16
                        R(m,n)=(1/(Lf*std(ii)))*(abs(d(n,ii)-[MV(m,1),MV(m,2)]* [ww(n,:),zeros(1,Lh);zeros(1,Lh),ww(n,:)]*htil(:,ii))^2 + ...
                            real([MV(m,1),MV(m,2)]*[dsigma(n,n),dsigma(n+Lf,n);dsigma(n,n+Lf),dsigma(n+Lf,n+Lf)]*[MV(m,1),MV(m,2)]'))-log(1/16);
                        if R(m,n)<Rmini
                            min=m;
                            Rmini=R(m,n);
                        end;
                    end;
                    mini(n,ii)=min;
                    ddn(Nt*(n-1)+1,ii)=MV(min,1);
                    xx1(n,ii)=ddn(Nt*(n-1)+1,ii);
                    ddn(n*Nt,ii)=MV(min,2);
                    xx2(n,ii)=ddn(n*Nt,ii); 
                end;
                T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
                % Calcul des matrices htilde et sigmatildeh %
                Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
                htil(:,ii)=Q*d(:,ii);
                sigmatil=sigmah-Q*T*W*sigmah;
                dsigma=W*sigmatil*W';
                ll=ll+1;
            end;
            % Fin de la boucle EM et calcul des LLR's %
            for n=1:Lf
                delta1(n,ii)=log(exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
                    log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n)));
                delta2(n,ii)=log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n)));
                delta3(n,ii)=log(exp(-R(3,n))+exp(-R(4,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(15,n))+exp(-R(16,n)))-...
                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(5,n))+exp(-R(6,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(13,n))+exp(-R(14,n)));
                delta4(n,ii)=log(exp(-R(2,n))+exp(-R(4,n))+exp(-R(6,n))+exp(-R(8,n))+exp(-R(10,n))+exp(-R(12,n))+exp(-R(14,n))+exp(-R(16,n)))-...
                    log(exp(-R(1,n))+exp(-R(3,n))+exp(-R(5,n))+exp(-R(7,n))+exp(-R(9,n))+exp(-R(11,n))+exp(-R(13,n))+exp(-R(15,n)));
                zeta((n-1)*4+1,ii)=delta1(n,ii);
                zeta((n-1)*4+2,ii)=delta2(n,ii);
                zeta((n-1)*4+3,ii)=delta3(n,ii);
                zeta((n-1)*4+4,ii)=delta4(n,ii);
            end;
        end;  
        mini;
        delta=reshape(zeta,block_size,1);
        % Calcul des probabilités P0 et P1 %
        P0=1./(1+exp(delta));
        for i=1:block_size
            % On évite les problèmes de log 0 %
            if P0(i)>=0.99999
                P0(i)=0.99999;
            elseif P0(i)<=0.00001
                P0(i)=0.00001;
            end;
        end;
        P1=1-P0;
        [z_hat,success,k,APP]=ldpc_decode(P0,P1,H);
        k;
        success;
        extrinsic=APP-delta;
        P0=1./(1+exp(extrinsic));
        for i=1:block_size
            % On évite les problèmes de log 0 %
            if P0(i)>=0.99999
                P0(i)=0.99999;
            elseif P0(i)<=0.00001
                P0(i)=0.00001;
            end;
        end;
        P1=1-P0;
        P1;
        x_hat = z_hat(size(G,2)+1-size(G,1):size(G,2));
        x_hat = x_hat';
        o = xor(x_hat,gseq);
        numerror = sum(o);
        BER(da,num_packet) = numerror/(taille); 
        % fin de la première turbo-itération: on commence la seconde itération turbo %
        % on utilise les loglikelihood ratios au niveau des bits en sortie %
        % du LDPC pour les convertir en rapports de vraisemblance au niveau %
        % des symboles %
        % On commence par calculer les nouvelles probabilités à priori sur les symboles %
        % à partir des APP extrinsèques du LDPC %
        
        for i=1:(block_size/(M*Nt))
            sP0(:,i)=P0(M*Nt*(i-1)+1:M*Nt*i,1);
            sP1(:,i)=P1(M*Nt*(i-1)+1:M*Nt*i,1); 
            P(1,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP0(4,i));
            P(2,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP1(4,i));
            P(3,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP0(4,i));
            P(4,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP1(4,i));
            P(5,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP0(4,i));
            P(6,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP1(4,i));
            P(7,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP0(4,i));
            P(8,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP1(4,i));
            P(9,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP0(4,i));
            P(10,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP1(4,i));
            P(11,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP0(4,i));
            P(12,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP1(4,i));
            P(13,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP0(4,i));
            P(14,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP1(4,i));
            P(15,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP0(4,i));
            P(16,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP1(4,i));
            sumP=sum(exp(P(:,i)),1);
            P(1:16,i)=exp(P(1:16,i))/sumP;
        end;
        
        PP=reshape(P,16*Lf,block_size/(M*Nt*Lf));
        PQ=reshape(PP,16,Lf*block_size/(M*Nt*Lf));
        for i=1:Lf*block_size/(M*Nt*Lf)
            [ZX(i),m]=max(PQ(:,i));
            N(i)=m;
        end;
        ZX=reshape(N,Lf,block_size/(M*Nt*Lf));
        for ii=1:(block_size/(M*Nt*Lf))
            for nn=1:Lf
                xx1(nn,ii)=MV(ZX(nn,ii),1);
                xx2(nn,ii)=MV(ZX(nn,ii),2);
            end;
        end;
        
        % On recommence l'algorithme EM %
        for ii=1:(block_size/(M*Nt*Lf))
            T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
            % Calcul des matrices htilde et sigmatildeh %
            Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
            sigmatil=sigmah-Q*T*W*sigmah;
            dsigma=W*sigmatil*W';
            for ll=1:EM_iter
                % On calcule les métriques de l'E.M. %
                for n=1:Lf
                    Rmini=10^5;
                    for m=1:16
                        R(m,n)=(1/(Lf*std(ii)))*(abs(d(n,ii)-[MV(m,1),MV(m,2)]* [ww(n,:),zeros(1,Lh);zeros(1,Lh),ww(n,:)]*htil(:,ii))^2 + ...
                            real([MV(m,1),MV(m,2)]*[dsigma(n,n),dsigma(n+Lf,n);dsigma(n,n+Lf),dsigma(n+Lf,n+Lf)]*[MV(m,1),MV(m,2)]'))-log(PP(m+16*(n-1),ii));
                        if R(m,n)<Rmini
                            min=m;
                            Rmini=R(m,n);
                        end;
                    end;
                    ddn(Nt*(n-1)+1,ii)=MV(min,1);
                    xx1(n,ii)=ddn(Nt*(n-1)+1,ii);
                    ddn(n*Nt,ii)=MV(min,2);
                    xx2(n,ii)=ddn(n*Nt,ii); 
                end;
                T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
                % Calcul des matrices htilde et sigmatildeh %
                Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
                htil(:,ii)=Q*d(:,ii);
                sigmatil=sigmah-Q*T*W*sigmah;
                dsigma=W*sigmatil*W';
                ll=ll+1;
            end;
            % Fin de la boucle EM, on recalcule les rapports de vraisemblance %
            for n=1:Lf
                delta1(n,ii)=log(exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
                    log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n)));
                delta2(n,ii)=log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n)));
                delta3(n,ii)=log(exp(-R(3,n))+exp(-R(4,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(15,n))+exp(-R(16,n)))-...
                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(5,n))+exp(-R(6,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(13,n))+exp(-R(14,n)));
                delta4(n,ii)=log(exp(-R(2,n))+exp(-R(4,n))+exp(-R(6,n))+exp(-R(8,n))+exp(-R(10,n))+exp(-R(12,n))+exp(-R(14,n))+exp(-R(16,n)))-...
                    log(exp(-R(1,n))+exp(-R(3,n))+exp(-R(5,n))+exp(-R(7,n))+exp(-R(9,n))+exp(-R(11,n))+exp(-R(13,n))+exp(-R(15,n)));
                zeta((n-1)*4+1,ii)=delta1(n,ii);
                zeta((n-1)*4+2,ii)=delta2(n,ii);
                zeta((n-1)*4+3,ii)=delta3(n,ii);
                zeta((n-1)*4+4,ii)=delta4(n,ii);
            end;
        end;
        delta=reshape(zeta,block_size,1);
        % Calcul des informations extrinsèques %
        delta=delta-extrinsic;
        % Calcul des probabilités P0 et P1 %
        P0=1./(1+exp(delta));
        for i=1:block_size
            % On évite les problèmes de log 0 %
            if P0(i)>=0.99999
                P0(i)=0.99999;
            elseif P0(i)<=0.00001
                P0(i)=0.00001;
            end;
        end;
        P1=1-P0;
        P1;
        [z_hat,success,k,APP]=ldpc_decode(P0,P1,H);
        success;
        k;
        extrinsic=APP-delta;
        P0=1./(1+exp(extrinsic));
        for i=1:block_size
            % On évite les problèmes de log 0 %
            if P0(i)>=0.99999
                P0(i)=0.99999;
            elseif P0(i)<=0.00001
                P0(i)=0.00001;
            end;
        end;
        P1=1-P0;
        x_hat = z_hat(size(G,2)+1-size(G,1):size(G,2));
        x_hat = x_hat';
        o = xor(x_hat,gseq);
        numerror1 = sum(o);
        BER1(da,num_packet) = numerror1/(taille); 
        % fin de la seconde turbo-itération: on commence la troisième itération turbo %
        % on utilise les loglikelihood ratios au niveau des bits en sortie %
        % du LDPC pour les convertir en rapports de vraisemblance au niveau %
        % des symboles %
        % On commence par calculer les nouvelles probabilités à priori sur les symboles %
        % à partir des APP extrinsèques du LDPC %
        
        
        for i=1:(block_size/(M*Nt))
            sP0(:,i)=P0(M*Nt*(i-1)+1:M*Nt*i,1);
            sP1(:,i)=P1(M*Nt*(i-1)+1:M*Nt*i,1); 
            P(1,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP0(4,i));
            P(2,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP1(4,i));
            P(3,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP0(4,i));
            P(4,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP1(4,i));
            P(5,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP0(4,i));
            P(6,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP1(4,i));
            P(7,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP0(4,i));
            P(8,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP1(4,i));
            P(9,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP0(4,i));
            P(10,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP1(4,i));
            P(11,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP0(4,i));
            P(12,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP1(4,i));
            P(13,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP0(4,i));
            P(14,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP1(4,i));
            P(15,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP0(4,i));
            P(16,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP1(4,i));
            sumP=sum(exp(P(:,i)),1);
            P(1:16,i)=exp(P(1:16,i))/sumP;
        end;
        PP=reshape(P,16*Lf,block_size/(M*Nt*Lf));
        PQ=reshape(PP,16,Lf*block_size/(M*Nt*Lf));
        for i=1:Lf*block_size/(M*Nt*Lf)
            [ZX(i),m]=max(PQ(:,i));
            N(i)=m;
        end;
        ZX=reshape(N,Lf,block_size/(M*Nt*Lf));
        for ii=1:(block_size/(M*Nt*Lf))
            for nn=1:Lf
                xx1(nn,ii)=MV(ZX(nn,ii),1);
                xx2(nn,ii)=MV(ZX(nn,ii),2);
            end;
        end;
        
        % On recommence l'algorithme EM %
        for ii=1:(block_size/(M*Nt*Lf))
            T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
            % Calcul des matrices htilde et sigmatildeh %
            Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
            sigmatil=sigmah-Q*T*W*sigmah;
            dsigma=W*sigmatil*W';
            for ll=1:EM_iter
                % On calcule les métriques de l'E.M. %
                for n=1:Lf
                    Rmini=10^5;
                    for m=1:16
                        R(m,n)=(1/(Lf*std(ii)))*(abs(d(n,ii)-[MV(m,1),MV(m,2)]* [ww(n,:),zeros(1,Lh);zeros(1,Lh),ww(n,:)]*htil(:,ii))^2 + ...
                            real([MV(m,1),MV(m,2)]*[dsigma(n,n),dsigma(n+Lf,n);dsigma(n,n+Lf),dsigma(n+Lf,n+Lf)]*[MV(m,1),MV(m,2)]'))-log(PP(m+16*(n-1),ii));
                        if R(m,n)<Rmini
                            min=m;
                            Rmini=R(m,n);
                        end;
                    end;
                    ddn(Nt*(n-1)+1,ii)=MV(min,1);
                    xx1(n,ii)=ddn(Nt*(n-1)+1,ii);
                    ddn(n*Nt,ii)=MV(min,2);
                    xx2(n,ii)=ddn(n*Nt,ii); 
                end;
                T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
                % Calcul des matrices htilde et sigmatildeh %
                Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
                htil(:,ii)=Q*d(:,ii);
                sigmatil=sigmah-Q*T*W*sigmah;
                dsigma=W*sigmatil*W';
                ll=ll+1;
            end;
            % Calcul des LLR's %
            for n=1:Lf
                delta1(n,ii)=log(exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
                    log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n)));
                delta2(n,ii)=log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n)));
                delta3(n,ii)=log(exp(-R(3,n))+exp(-R(4,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(15,n))+exp(-R(16,n)))-...
                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(5,n))+exp(-R(6,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(13,n))+exp(-R(14,n)));
                delta4(n,ii)=log(exp(-R(2,n))+exp(-R(4,n))+exp(-R(6,n))+exp(-R(8,n))+exp(-R(10,n))+exp(-R(12,n))+exp(-R(14,n))+exp(-R(16,n)))-...
                    log(exp(-R(1,n))+exp(-R(3,n))+exp(-R(5,n))+exp(-R(7,n))+exp(-R(9,n))+exp(-R(11,n))+exp(-R(13,n))+exp(-R(15,n)));
                zeta((n-1)*4+1,ii)=delta1(n,ii);
                zeta((n-1)*4+2,ii)=delta2(n,ii);
                zeta((n-1)*4+3,ii)=delta3(n,ii);
                zeta((n-1)*4+4,ii)=delta4(n,ii);
            end;
        end;
        delta=reshape(zeta,block_size,1);
        % Calcul des informations extrinsèques %
        delta=delta-extrinsic;
        % Calcul des probabilités P0 et P1 %
        P0=1./(1+exp(delta));
        for i=1:block_size
            % On évite les problèmes de log 0 %
            if P0(i)>=0.99999
                P0(i)=0.99999;
            elseif P0(i)<=0.00001
                P0(i)=0.00001;
            end;
        end;
        P1=1-P0;
        P1;
        [z_hat,success,k,APP]=ldpc_decode(P0,P1,H);
        extrinsic=APP-delta;
        P0=1./(1+exp(extrinsic));
        for i=1:block_size
            % On évite les problèmes de log 0 %
            if P0(i)>=0.99999
                P0(i)=0.99999;
            elseif P0(i)<=0.00001
                P0(i)=0.00001;
            end;
        end;
        P1=1-P0;
        x_hat = z_hat(size(G,2)+1-size(G,1):size(G,2));
        x_hat = x_hat';
        o = xor(x_hat,gseq);
        numerror2 = sum(o);
        BER2(da,num_packet) = numerror2/(taille); 
        
        % fin de la troisième turbo-itération: on commence la quatrième itération turbo %
        % on utilise les loglikelihood ratios au niveau des bits en sortie %
        % du LDPC pour les convertir en rapports de vraisemblance au niveau %
        % des symboles %
        % On commence par calculer les nouvelles probabilités à priori sur les symboles %
        % à partir des APP extrinsèques du LDPC %
        
        
        for i=1:(block_size/(M*Nt))
            sP0(:,i)=P0(M*Nt*(i-1)+1:M*Nt*i,1);
            sP1(:,i)=P1(M*Nt*(i-1)+1:M*Nt*i,1); 
            P(1,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP0(4,i));
            P(2,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP1(4,i));
            P(3,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP0(4,i));
            P(4,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP1(4,i));
            P(5,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP0(4,i));
            P(6,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP1(4,i));
            P(7,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP0(4,i));
            P(8,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP1(4,i));
            P(9,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP0(4,i));
            P(10,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP1(4,i));
            P(11,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP0(4,i));
            P(12,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP1(4,i));
            P(13,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP0(4,i));
            P(14,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP1(4,i));
            P(15,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP0(4,i));
            P(16,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP1(4,i));
            sumP=sum(exp(P(:,i)),1);
            P(1:16,i)=exp(P(1:16,i))/sumP;
        end;
        PP=reshape(P,16*Lf,block_size/(M*Nt*Lf));
        PQ=reshape(PP,16,Lf*block_size/(M*Nt*Lf));
        for i=1:Lf*block_size/(M*Nt*Lf)
            [ZX(i),m]=max(PQ(:,i));
            N(i)=m;
        end;
        ZX=reshape(N,Lf,block_size/(M*Nt*Lf));
        for ii=1:(block_size/(M*Nt*Lf))
            for nn=1:Lf
                xx1(nn,ii)=MV(ZX(nn,ii),1);
                xx2(nn,ii)=MV(ZX(nn,ii),2);
            end;
        end;
        
        % On recommence l'algorithme EM %
        for ii=1:(block_size/(M*Nt*Lf))
            T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
            % Calcul des matrices htilde et sigmatildeh %
            Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
            sigmatil=sigmah-Q*T*W*sigmah;
            dsigma=W*sigmatil*W';
            for ll=1:EM_iter
                % On calcule les métriques de l'E.M. %
                for n=1:Lf
                    Rmini=10^5;
                    for m=1:16
                        R(m,n)=(1/(Lf*std(ii)))*(abs(d(n,ii)-[MV(m,1),MV(m,2)]* [ww(n,:),zeros(1,Lh);zeros(1,Lh),ww(n,:)]*htil(:,ii))^2 + ...
                            real([MV(m,1),MV(m,2)]*[dsigma(n,n),dsigma(n+Lf,n);dsigma(n,n+Lf),dsigma(n+Lf,n+Lf)]*[MV(m,1),MV(m,2)]'))-log(PP(m+16*(n-1),ii));
                        if R(m,n)<Rmini
                            min=m;
                            Rmini=R(m,n);
                        end;
                    end;
                    ddn(Nt*(n-1)+1,ii)=MV(min,1);
                    xx1(n,ii)=ddn(Nt*(n-1)+1,ii);
                    ddn(n*Nt,ii)=MV(min,2);
                    xx2(n,ii)=ddn(n*Nt,ii); 
                end;
                T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
                % Calcul des matrices htilde et sigmatildeh %
                Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
                htil(:,ii)=Q*d(:,ii);
                sigmatil=sigmah-Q*T*W*sigmah;
                dsigma=W*sigmatil*W';
                ll=ll+1;
            end;
            % Calcul des LLR's %
            for n=1:Lf
                delta1(n,ii)=log(exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
                    log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n)));
                delta2(n,ii)=log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n)));
                delta3(n,ii)=log(exp(-R(3,n))+exp(-R(4,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(15,n))+exp(-R(16,n)))-...
                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(5,n))+exp(-R(6,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(13,n))+exp(-R(14,n)));
                delta4(n,ii)=log(exp(-R(2,n))+exp(-R(4,n))+exp(-R(6,n))+exp(-R(8,n))+exp(-R(10,n))+exp(-R(12,n))+exp(-R(14,n))+exp(-R(16,n)))-...
                    log(exp(-R(1,n))+exp(-R(3,n))+exp(-R(5,n))+exp(-R(7,n))+exp(-R(9,n))+exp(-R(11,n))+exp(-R(13,n))+exp(-R(15,n)));
                zeta((n-1)*4+1,ii)=delta1(n,ii);
                zeta((n-1)*4+2,ii)=delta2(n,ii);
                zeta((n-1)*4+3,ii)=delta3(n,ii);
                zeta((n-1)*4+4,ii)=delta4(n,ii);
            end;
        end;
        delta=reshape(zeta,block_size,1);
        % Calcul des informations extrinsèques %
        delta=delta-extrinsic;
        % Calcul des probabilités P0 et P1 %
        P0=1./(1+exp(delta));
        for i=1:block_size
            % On évite les problèmes de log 0 %
            if P0(i)>=0.99999
                P0(i)=0.99999;
            elseif P0(i)<=0.00001
                P0(i)=0.00001;
            end;
        end;
        P1=1-P0;
        P1;
        [z_hat,success,k,APP]=ldpc_decode(P0,P1,H);
        extrinsic=APP-delta;
        P0=1./(1+exp(extrinsic));
        for i=1:block_size
            % On évite les problèmes de log 0 %
            if P0(i)>=0.99999
                P0(i)=0.99999;
            elseif P0(i)<=0.00001
                P0(i)=0.00001;
            end;
        end;
        P1=1-P0;
        x_hat = z_hat(size(G,2)+1-size(G,1):size(G,2));
        x_hat = x_hat';
        o = xor(x_hat,gseq);
        numerror3 = sum(o);
        BER3(da,num_packet) = numerror3/(taille);     
        
        % fin de la quatrième turbo-itération: on commence la quatrième itération turbo %
        % on utilise les loglikelihood ratios au niveau des bits en sortie %
        % du LDPC pour les convertir en rapports de vraisemblance au niveau %
        % des symboles %
        % On commence par calculer les nouvelles probabilités à priori sur les symboles %
        % à partir des APP extrinsèques du LDPC %
        
        
        for i=1:(block_size/(M*Nt))
            sP0(:,i)=P0(M*Nt*(i-1)+1:M*Nt*i,1);
            sP1(:,i)=P1(M*Nt*(i-1)+1:M*Nt*i,1); 
            P(1,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP0(4,i));
            P(2,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP1(4,i));
            P(3,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP0(4,i));
            P(4,i)=log(sP0(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP1(4,i));
            P(5,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP0(4,i));
            P(6,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP1(4,i));
            P(7,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP0(4,i));
            P(8,i)=log(sP0(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP1(4,i));
            P(9,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP0(4,i));
            P(10,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP0(3,i))+log(sP1(4,i));
            P(11,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP0(4,i));
            P(12,i)=log(sP1(1,i))+log(sP0(2,i))+log(sP1(3,i))+log(sP1(4,i));
            P(13,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP0(4,i));
            P(14,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP0(3,i))+log(sP1(4,i));
            P(15,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP0(4,i));
            P(16,i)=log(sP1(1,i))+log(sP1(2,i))+log(sP1(3,i))+log(sP1(4,i));
            sumP=sum(exp(P(:,i)),1);
            P(1:16,i)=exp(P(1:16,i))/sumP;
        end;
        PP=reshape(P,16*Lf,block_size/(M*Nt*Lf));
        PQ=reshape(PP,16,Lf*block_size/(M*Nt*Lf));
        for i=1:Lf*block_size/(M*Nt*Lf)
            [ZX(i),m]=max(PQ(:,i));
            N(i)=m;
        end;
        ZX=reshape(N,Lf,block_size/(M*Nt*Lf));
        for ii=1:(block_size/(M*Nt*Lf))
            for nn=1:Lf
                xx1(nn,ii)=MV(ZX(nn,ii),1);
                xx2(nn,ii)=MV(ZX(nn,ii),2);
            end;
        end;
        
        % On recommence l'algorithme EM %
        for ii=1:(block_size/(M*Nt*Lf))
            T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
            % Calcul des matrices htilde et sigmatildeh %
            Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
            sigmatil=sigmah-Q*T*W*sigmah;
            dsigma=W*sigmatil*W';
            for ll=1:EM_iter
                % On calcule les métriques de l'E.M. %
                for n=1:Lf
                    Rmini=10^5;
                    for m=1:16
                        R(m,n)=(1/(Lf*std(ii)))*(abs(d(n,ii)-[MV(m,1),MV(m,2)]* [ww(n,:),zeros(1,Lh);zeros(1,Lh),ww(n,:)]*htil(:,ii))^2 + ...
                            real([MV(m,1),MV(m,2)]*[dsigma(n,n),dsigma(n+Lf,n);dsigma(n,n+Lf),dsigma(n+Lf,n+Lf)]*[MV(m,1),MV(m,2)]'))-log(PP(m+16*(n-1),ii));
                        if R(m,n)<Rmini
                            min=m;
                            Rmini=R(m,n);
                        end;
                    end;
                    ddn(Nt*(n-1)+1,ii)=MV(min,1);
                    xx1(n,ii)=ddn(Nt*(n-1)+1,ii);
                    ddn(n*Nt,ii)=MV(min,2);
                    xx2(n,ii)=ddn(n*Nt,ii); 
                end;
                T=[diag(xx1(:,ii)),diag(xx2(:,ii))];
                % Calcul des matrices htilde et sigmatildeh %
                Q=inv((1/std(ii))*W'*T'*T*W+pinv(sigmah))*W'*T'*(1/std(ii));
                htil(:,ii)=Q*d(:,ii);
                sigmatil=sigmah-Q*T*W*sigmah;
                dsigma=W*sigmatil*W';
                ll=ll+1;
            end;
            % Calcul des LLR's %
            for n=1:Lf
                delta1(n,ii)=log(exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
                    log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n)));
                delta2(n,ii)=log(exp(-R(5,n))+exp(-R(6,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(13,n))+exp(-R(14,n))+exp(-R(15,n))+exp(-R(16,n)))-...
                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(3,n))+exp(-R(4,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(11,n))+exp(-R(12,n)));
                delta3(n,ii)=log(exp(-R(3,n))+exp(-R(4,n))+exp(-R(7,n))+exp(-R(8,n))+exp(-R(11,n))+exp(-R(12,n))+exp(-R(15,n))+exp(-R(16,n)))-...
                    log(exp(-R(1,n))+exp(-R(2,n))+exp(-R(5,n))+exp(-R(6,n))+exp(-R(9,n))+exp(-R(10,n))+exp(-R(13,n))+exp(-R(14,n)));
                delta4(n,ii)=log(exp(-R(2,n))+exp(-R(4,n))+exp(-R(6,n))+exp(-R(8,n))+exp(-R(10,n))+exp(-R(12,n))+exp(-R(14,n))+exp(-R(16,n)))-...
                    log(exp(-R(1,n))+exp(-R(3,n))+exp(-R(5,n))+exp(-R(7,n))+exp(-R(9,n))+exp(-R(11,n))+exp(-R(13,n))+exp(-R(15,n)));
                zeta((n-1)*4+1,ii)=delta1(n,ii);
                zeta((n-1)*4+2,ii)=delta2(n,ii);
                zeta((n-1)*4+3,ii)=delta3(n,ii);
                zeta((n-1)*4+4,ii)=delta4(n,ii);
            end;
        end;
        delta=reshape(zeta,block_size,1);
        % Calcul des informations extrinsèques %
        delta=delta-extrinsic;
        % Calcul des probabilités P0 et P1 %
        P0=1./(1+exp(delta));
        for i=1:block_size
            % On évite les problèmes de log 0 %
            if P0(i)>=0.99999
                P0(i)=0.99999;
            elseif P0(i)<=0.00001
                P0(i)=0.00001;
            end;
        end;
        P1=1-P0;
        P1;
        [z_hat,success,k,APP]=ldpc_decode(P0,P1,H);
        extrinsic=APP-delta;
        P0=1./(1+exp(extrinsic));
        for i=1:block_size
            % On évite les problèmes de log 0 %
            if P0(i)>=0.99999
                P0(i)=0.99999;
            elseif P0(i)<=0.00001
                P0(i)=0.00001;
            end;
        end;
        P1=1-P0;
        x_hat = z_hat(size(G,2)+1-size(G,1):size(G,2));
        x_hat = x_hat';
        o = xor(x_hat,gseq);
        numerror4 = sum(o);
        BER4(da,num_packet) = numerror4/(taille);    
    end;  
    TEB(da)=sum(BER(da,:),2)/Nbpacket;
    TEB1(da)=sum(BER1(da,:),2)/Nbpacket;
    TEB2(da)=sum(BER2(da,:),2)/Nbpacket;
    TEB3(da)=sum(BER3(da,:),2)/Nbpacket;
    TEB4(da)=sum(BER4(da,:),2)/Nbpacket;
end;

grid,x = 13:1.:14.;
grid on ;
grid,semilogy(x,TEB(:),'b', x, TEB1(:),'b*', x, TEB2(:),'r', x, TEB3(:),'go', x, TEB4(:),'b--');

return;














