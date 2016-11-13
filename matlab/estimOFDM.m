function  TEB = estimOFDM
% Dans cette version, le canal est estimé à l'aide de données non-codées.
% Les versions plus sophistiquées utilisent des données décidées en codage
% convolutif ou l'algorithme du urvivant
Nbpacket = 80;
SNR = [1,2,3,4,5,6,7,8];
num_error(1:length(SNR),1:Nbpacket)=0;
for da = 1:length(SNR)
    for num_packet=1:Nbpacket
% M nombre d'états de la modulation, ici une MDP4 %
n=100;   %nombre d'échantillons temporels
M=2;
Lf=5; %nombre d'échantillons fréquentiels = taille de la FFT 
%frequence Doppler maximale dans le specte fréquentiel du canal 
fd=100;
feet=n*M*Lf;
seq(:,num_packet)=(rand(1,feet)>0.5)';
mseq=seq(:,num_packet);
training_length=10;
% creation de la modulation MDP4 %
for i=1:n*Lf
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

% creation de la modulation multiporteuse OFDM %
for i=1:n
   w(:,i)=conj((vn(Lf*(i-1)+1:Lf*i))');
end;

% construction du canal de propagation bi-dimensionnel
% temps-espace, sélectif en temps et en fréquence 
x=fadingtaps(fd,1,1,n,512);

% On construit les répliques fréquentielles décalées de ce trajet multiple
% On construit Lf colonnes chacune ayant Ns (nombre symboles) lignes 
% Ainsi, on obtient une matrice temps fréquence (colonnes => fréquences et 
% lignes => temps)
for k=1:Lf
    for ii=1:n
        y(k,ii)=x(ii);
    end;
end;

for k=1:Lf
    ww(k,1)=exp((k-1)*2*j*pi/Lf);
    zy(k,:)=y(k,:).*ww(k,1);
end;

for i=1:n
    dd(:,i)=w(:,i).*zy(:,i);
    [d(:,i),std]=awgn(1,Lf,SNR(da),dd(:,i));
end;
sigma=std;
        
 % Calcul temporel des coefficients du canal de propagation, on suppose 
 % que les "training_length" premiers symboles sont connus, on utilise 
 % un estimateur des moindres carrés
 for k=1:training_length
     yob((k-1)*Lf+1:k*Lf,1)=d(:,k);
     A((k-1)*Lf+1:k*Lf,:)=diag(w(:,k));
 end;
 h=inv(A'*A+sigma*eye(Lf))*A'*yob;
 
 % Calcul de la fonction d'autocorrélation fréquentielle, on la calcule 
 % à l'aide de la fonction xcorr
 Sf=xcorr(h);
% Construction de la matrice de Toeplitz de corrélation fréquentielle du
% canal
Af=toeplitz(Sf(Lf:-1:1),Sf(Lf:1:2*Lf-1));
[U, S, V]=svd (Af);
% Calcul des coefficients de la matrice diagonale Phi(w)
% On utilise la distribution de Jakes pour le spectre Doppler
for m=1:Lf
a(m)=2*S(m,m)/(2*pi*fd*sigma);
if (a(m)>=1)
        b(m)=(pi/2)*a(m)-sqrt(a(m)^2-1)*(pi/2-asin(1/a(m)));
    else
        b(m)=(pi/2)*a(m)+sqrt(1-a(m)^2)*log((1+sqrt(1-a(m)^2))/a(m));
    end;
   gg(m)=exp(-2*fd*(b(m)+log(a(m)/2)));
   f(m)=1-gg(m);
end;

% Construction de la matrice de l'estimateur fréquentiel 
C=diag(f);
nn=training_length+1;
D=U*C*V'*h;
for i=1:Lf 
    CF(i)=D(i)/abs(D(i));
end;
for i=1:nn-1
    db(:,i)=w(:,i);
end;

% On opére maintenant en mode DD: Decision-Directed
dn(:,nn)=CF'.*d(:,nn);
for k=1:Lf
          if ((real(dn(k,nn))>0) & (imag(dn(k,nn))>0))
             db(k,nn)=exp(j*pi/4);
         elseif((real(dn(k,nn))>0) & (imag(dn(k,nn))<0))
                db(k,nn)=exp(j*(pi/4+3*pi/2));
            elseif((real(dn(k,nn))<0) & (imag(dn(k,nn))>0))
                db(k,nn)=exp(j*(pi/4+pi/2));
            else 
                db(k,nn)=exp(j*(pi/4+pi));
            end;
        end;
        
for ii=nn:n-1
    [C,CF]=estim_canal(db(:,ii-nn+2:1:ii),d(:,ii-nn+2:1:ii),nn-1,Lf,sigma,fd);
    dn(:,ii+1)=CF'.*d(:,ii+1);
for k=1:Lf
          if ((real(dn(k,ii+1))>0) & (imag(dn(k,ii+1))>0))
             db(k,ii+1)=exp(j*pi/4);
         elseif((real(dn(k,ii+1))>0) & (imag(dn(k,ii+1))<0))
                db(k,ii+1)=exp(j*(pi/4+3*pi/2));
            elseif((real(dn(k,ii+1))<0) & (imag(dn(k,ii+1))>0))
                db(k,ii+1)=exp(j*(pi/4+pi/2));
            else 
                db(k,ii+1)=exp(j*(pi/4+pi));
            end;
        end;
    end;
    
    % Calcul du taux d'erreur par symbole ou par élément binaire
for i=nn:n 
   err(:,i)=db(:,i)-w(:,i);
end;
for i=nn:n
    for kk=1:Lf
        if err(kk,i)~=0
            num_error(SNR(da),num_packet)=num_error(SNR(da),num_packet)+1;
        end;
    end;
end;
end;
end;

    %E(i)=D(i)-zy(i,training_length+1);
    %FF(i)=h(i)-zy(i,training_length+1);
    %E(i)=abs(E(i))^2;
    %FF(i)=abs(FF(i))^2;
    %end;
    %QQ=sum(E)/Lf;
    %QR=sum(FF)/Lf;
end
TEB(1:length(SNR))=sum(num_error,2)/(Nbpacket*Lf*n);
grid,x = 1.:1.:8.;
grid on ;
grid,semilogy(x,TEB(:),'b');
grid on ;
grid,title('BER versus Eb/N0')
    
return