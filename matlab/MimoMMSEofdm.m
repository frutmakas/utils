function  [TEB,SER] = MIMO1OFDM
% On considère une modulation MDP-4 transmise en OFDM %
% sur deux antennes d'émission en utilisant le schéma d'Alamouti %
% M=((s0, s1),(-s1* s0*)) %
n = 100;
L = 8;
% taille d'un paquet en nombre de symboles MDP-4 %
M = 2; 
Nbpacket = 100;
Lf = 8;
% taille de la FFT: Lf doit être un multiple du nombre d'antennes d'émission %
Nt = 2;
% nombre d'antennes de transmission %
training_length=5;
% frequence Doppler maximale dans le specte fréquentiel du canal %
fd = 50;
coding_rate = 1;
nb_step = 10;
Gn=[];
Gn1=[];
Gn2=[];
SNR = [1,2,3,4,5,6];
for xx=1:length(SNR)
    for jj=1:Nbpacket
        num_error(xx,jj)=0;  
    end;
end;
for dd=1:length(SNR)
    SER(dd,1)=0;
end;

for da=1:length(SNR)
    for num_packet=1:Nbpacket
        feet=n*M*Lf*Nt;
        seq(:,num_packet)=(rand(1,feet)>0.5)';
        mseq=seq(:,num_packet);
        
        for i=1:(n*Lf*Nt)
            sn=mseq(M*(i-1)+1:M*i);
            if(sn==[1 1]')
                vn(i)=exp(j*pi/4);
            elseif(sn==[0 1]')
                vn(i)=exp(j*(pi/4+pi/2));
            elseif(sn==[0 0]')
                vn(i)=exp(j*(pi/4+pi));
            else
                vn(i)=exp(j*(pi/4+3*pi/2));
            end;
        end;
        
        % création du codage spatio-temporel % 
        for jj=1:(n*Lf)
            wn(Nt*(jj-1)+1:Nt*jj)=vn(Nt*(jj-1)+1:Nt*jj);
            tn(Nt*(jj-1)+1)=-conj(vn(Nt*jj));
            tn(Nt*jj)=conj(vn(Nt*(jj-1)+1));
            s1(Nt*(jj-1)+1)=vn(Nt*(jj-1)+1);  
            s1(Nt*jj)=tn(Nt*(jj-1)+1);
            s2(Nt*(jj-1)+1)=vn(Nt*jj);
            s2(Nt*jj)=tn(Nt*jj);
        end;
        
        % creation de la modulation multiporteuse OFDM %
        for i=1:n*Nt
            w1(:,i)=conj((s1(Lf*(i-1)+1:Lf*i))');
            w2(:,i)=conj((s2(Lf*(i-1)+1:Lf*i))');
            dcn(:,i)=conj((wn(Lf*(i-1)+1:Lf*i))');
        end;
        
        % On a à créer deux canaux de propagation sur chacun des deux flux %
        % de sortie des deux antennes %
        
        % première antenne %
        x1=fadingtaps(fd,1,1,n*Nt,512);
        % deuxième antenne %
        x2=fadingtaps(fd,1,1,n*Nt,512);
        
        % On construit les répliques fréquentielles décalées de ce trajet multiple
        % On construit Lf colonnes chacune ayant Ns (nombre symboles) lignes 
        % Ainsi, on obtient une matrice temps fréquence (colonnes => fréquences et 
        % lignes => temps)
        
        for k=1:Lf
            for ii=1:n*Nt
                y1(k,ii)=x1(ii);
                y2(k,ii)=x2(ii);  
            end;
            ww(k,1)=exp(-j*2*pi*((k-1)/Lf));
            zy1(k,:)=y1(k,:).*ww(k,1);
            zy2(k,:)=y2(k,:).*ww(k,1);
        end;
        
        % Création du signal X1,n.H1,n %
        for i=1:n*Nt
            dd1(:,i)=diag(w1(:,i))*zy1(:,i);
            dd2(:,i)=diag(w2(:,i))*zy2(:,i);
            [d(:,i),std]=awgn(coding_rate,Lf,SNR(da),dd1(:,i)+dd2(:,i));
        end;
        
        for i=1:L
            Gn=[Gn;diag(w1(:,i)),diag(w2(:,i))];
            rn((i-1)*Lf+1:i*Lf,1)=d(:,i);
            Gn1 = [Gn1;diag(w1(:,i))];   
            Gn2 = [Gn2;diag(w2(:,i))]; 
        end;
        hest=pinv(Gn)*rn;
        
        % Test de l'EM et de l'algorithme SAGE %
        % On initialise : algorithme EM %
        
        hest1EM = inv(Gn1'*Gn1)*Gn1'*rn;
        hest2EM = inv(Gn2'*Gn2)*Gn2'*rn;
        
        for ii = 1:nb_step
            % On commence par calculer Zi,n (cf article Georghiadhes IEEE Trans. on %
            % Commun. janvier 2003) %
            
            Z1 = Gn1*hest1EM;
            Z2 = Gn2*hest2EM;
            
            Y1 = Z1 + 0.5*(rn-Z1-Z2);
            Y2 = Z2 + 0.5*(rn-Z1-Z2);
            
            hest1EM = inv(Gn1'*Gn1)*Gn1'*Y1;
            hest2EM = inv(Gn2'*Gn2)*Gn2'*Y2;
            ii = ii+1;
        end;
        
        % algorithme SAGE %
        hest1sage = inv(Gn1'*Gn1)*Gn1'*rn;
        hest2sage = inv(Gn2'*Gn2)*Gn2'*rn;
        Z1 = Gn1*hest1sage;
        Z2 = Gn2*hest2sage;
        
        for nn = 1:nb_step
            % On commence par calculer Zi,n (cf article Georghiadhes IEEE Trans. on %
            % Commun. janvier 2003) %
            
            Y1 = rn-Z2;
            hest1sage = inv(Gn1'*Gn1)*Gn1'*Y1;
            K1 = Gn1*hest1sage;
            
            Y2 = rn-Z1;
            hest2sage = inv(Gn2'*Gn2)*Gn2'*Y2;
            K2 = Gn2*hest2sage;
            Z1=K1;
            Z2=K2;
            ii = ii+1;
        end;
        
        Gn=[];
        Gn1=[];
        Gn2=[];
        
        % A partir de maintenant l'estimateur opère en mode DD: Decision Directed %
        ddn(:,1:L)=dcn(:,1:L);
        for ii=L+1:n*Nt
            for k=1:Lf/Nt
                H(1,:)=[hest(Nt*(k-1)+1),hest(Nt*(k-1)+1+Lf)];   
                H(2,:)=[conj(hest(Nt*k+Lf)),-conj(hest(Nt*k))];  
                A=conj([d((k-1)*Nt+1,ii),conj(d(k*Nt,ii))]');
                %x((k-1)*Nt+1:k*Nt)=pinv(H)*A;
                x((k-1)*Nt+1)=(H(1,2)*A(2)-H(2,2)*A(1))/(H(2,1)*H(1,2)-H(1,1)*H(2,2));
                x(k*Nt)=(H(2,1)*A(1)-H(1,1)*A(2))/(H(1,2)*H(2,1)-H(1,1)*H(2,2));
                ddn((k-1)*Nt+1:k*Nt,ii)=conj([cos(pi/4)*(sign(real(x((k-1)*Nt+1)))+j*sign(imag((x((k-1)*Nt+1))))),...
                        cos(pi/4)*(sign(real(x(k*Nt)))+j*sign(imag(x(k*Nt))))]');
                w1((k-1)*Nt+1:k*Nt,ii)=conj([ddn((k-1)*Nt+1,ii),-conj(ddn(k*Nt,ii))]');
                w2((k-1)*Nt+1:k*Nt,ii)=conj([ddn(k*Nt,ii),conj(ddn((k-1)*Nt+1,ii))]');
            end; 
            % Une fois que l'on a déterminé le vecteur des symboles décidés, on réestime le canal % 
            % On recrée d'abord le codage STBC %
            for pp=1:L
                rn((pp-1)*Lf+1:pp*Lf,1)=d(:,pp+ii-L);
                Gn=[Gn;diag(w1(:,pp+ii-L)),diag(w2(:,pp+ii-L))];
            end;
            hest=pinv(Gn)*rn;
            Gn=[];
        end;
        % On compte les erreurs dans le paquet transmis %
        for i=L+1:n*Nt 
            err(:,i)=fix(dcn(:,i)-ddn(:,i));
        end;
        for i=L+1:n*Nt 
            for kk=1:Lf
                if abs(err(kk,i))>0.1
                    num_error(da,num_packet)=num_error(da,num_packet)+1;
                end;
            end;
        end; 
        TEB(da,num_packet)=num_error(da,num_packet)/((n*Nt-L)*Lf);
    end;
    SER(da) = sum(TEB(da,:),2)/Nbpacket;
end;
grid,x = 1.:1.:6.;
grid on ;
grid,semilogy(x,SER(:),'b');
grid on ;
grid,title('SER versus Es/N0');

return;













