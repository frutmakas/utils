function [TEB,SER] = mimowang
% On considère une modulation MDP-4 transmise en OFDM %
% sur deux antennes d'émission en utilisant le schéma d'Alamouti %
% M=((s0, s1),(-s1* s0*)) %
n = 100;
L = 10;
% taille d'un paquet en nombre de symboles MDP-4 %
M = 2; 
Nbpacket = 3;
Lf = 8;
% taille de la FFT: Lf doit être un multiple du nombre d'antennes d'émission %
Nt = 2;
% nombre d'antennes de transmission %
Lh=3;
% étalement des trajets multiples maximum %
fd = 100;
% frequence Doppler maximale dans le specte fréquentiel du canal %
SNR = 20;
coding_rate = 1;
nb_step = 5;
Gn1=[];
Gn2=[];
W=[];

SNR = [15,16,17,18,19,20];
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
        % de sortie des deux antennes: chacun étant sélectif en fréquence: %
        % le premier comporte deux coefficients et le deuxième trois %
        
        % première antenne %
        x1=0.8*fadingtaps(fd,1,1,n*Nt,512);
        x2=0.6*fadingtaps(fd,1,1,n*Nt,512);
        
        % deuxième antenne %
        x3=0.8*fadingtaps(fd,1,1,n*Nt,512);
        x4=0.4*fadingtaps(fd,1,1,n*Nt,512);
        x5=0.4*fadingtaps(fd,1,1,n*Nt,512);
        
        % on construit alors la matrice des coefficients de IFFT %
        for k=1:Lf
            for i=1:Lh
                ww(k,i)=(1/sqrt(Lf))*exp(-j*(2*pi*(k-1)*(i-1))/Lf);
            end;
        end;
        
        W=[ww,zeros(Lf,Lh);zeros(Lf,Lh),ww];    
        
        % on construit la matrice des coefficients de Rayleigh du canal % 
        
        for time=1:n*Nt 
            H(:,time)=conj([x1(time),0,x2(time),x3(time),x4(time),x5(time)]');
        end;
        
        % construction du signa	émis %
        
        for jj=1:n*Nt 
            X(:,jj)=[diag(w1(:,jj)),diag(w2(:,jj))]*W*H(:,jj);
            [d(:,jj),std]=awgn(coding_rate,Lf,SNR(da),X(:,jj));
        end;
        
        for mm=1:L
            rn((mm-1)*Lf+1:mm*Lf,1)=d(:,mm);
            Gn1 = [Gn1;diag(w1(:,mm))];   
            Gn2 = [Gn2;diag(w2(:,mm))]; 
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
            nn = nn+1;
        end;
        Gn1=[];
        Gn2=[];
        
        % A partir de maintenant l'estimateur opère en mode DD: Decision Directed %
        ddn(:,1:L)=dcn(:,1:L);
        for ii=L+1:n*Nt
            for k=1:Lf/Nt
                HH(1,:)=[hest1sage(Nt*(k-1)+1),hest2sage(Nt*(k-1)+1)];   
                HH(2,:)=[conj(hest2sage(Nt*k)),-conj(hest1sage(Nt*k))];  
                A=conj([d((k-1)*Nt+1,ii),conj(d(k*Nt,ii))]');
                x((k-1)*Nt+1)=(HH(1,2)*A(2)-HH(2,2)*A(1))/(HH(2,1)*HH(1,2)-HH(1,1)*HH(2,2));
                x(k*Nt)=(HH(2,1)*A(1)-HH(1,1)*A(2))/(HH(1,2)*HH(2,1)-HH(1,1)*HH(2,2));
                ddn((k-1)*Nt+1:k*Nt,ii)=conj([cos(pi/4)*(sign(real(x((k-1)*Nt+1)))+j*sign(imag((x((k-1)*Nt+1))))),...
                        cos(pi/4)*(sign(real(x(k*Nt)))+j*sign(imag(x(k*Nt))))]');
                w1((k-1)*Nt+1:k*Nt,ii)=conj([ddn((k-1)*Nt+1,ii),-conj(ddn(k*Nt,ii))]');
                w2((k-1)*Nt+1:k*Nt,ii)=conj([ddn(k*Nt,ii),conj(ddn((k-1)*Nt+1,ii))]');
            end; 
            % Une fois que l'on a déterminé le vecteur des symboles décidés, on réestime le canal % 
            for pp=1:L
                rn((pp-1)*Lf+1:pp*Lf,1)=d(:,pp+ii-L);
                Gn1 = [Gn1;diag(w1(:,pp+ii-L))];   
                Gn2 = [Gn2;diag(w2(:,pp+ii-L))]; 
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
                Z1 = K1;
                Z2 = K2;
                nn = nn+1;
            end;
            Gn1=[];
            Gn2=[];
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
grid,x = 15.:1.:20.;
grid on ;
grid,semilogy(x,SER(:),'b');
grid on ;
grid,title('SER versus Es/N0');

return;














