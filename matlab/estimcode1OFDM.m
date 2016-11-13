function  TEB = estimcode1OFDM
%M nombre d'états de la modulation, ici une MDP4 %
n=20;
g=[1 0 0 1 1; 1 1 1 0 1];
coding_rate=0.5;
[p,K]=size(g);
m=K-1;
M=2;
%taille de la FFT %
Lf=5;
training_length=5;
%frequence Doppler maximale dans le specte fréquentiel du canal 
fd=100;
Nbpacket=1000;
SNR = [1,2,3,4,5,6,7,8];
estim_length = 5;
num_error(1:length(SNR),1:Nbpacket)=0;
for da=1:length(SNR)
    for num_packet=1:Nbpacket
feet=n*M*Lf;
seq(:,num_packet)=(rand(1,feet)>0.5)';
mseq=seq(:,num_packet);

% codage de la séquence d'information mseq avec un code 
% convolutif de polynomes 23-35
seq_cod=encode_block(g, mseq');

% creation de la modulation MDP4 %
for i=1:(n*Lf*(1/coding_rate)+K-1)
   sn=(seq_cod(M*(i-1)+1:M*i))';
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
    [d(:,i),std]=awgn(coding_rate,Lf,SNR(da),dd(:,i));
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
% On utilise une profondeur de décodage (length_depth)
% égale ici à 5*Lf=25 symboles MDP4
for pp=1:n-nn+1-estim_length
for k=1:estim_length
dn(:,k)=CF'.*d(:,nn-2+k+pp);
end;
ysd=reshape(dn,estim_length*Lf,1); 
for i=1:estim_length*Lf             
qn(M*(i-1)+1:M*i)=[real(ysd(i)) imag(ysd(i))]';
end;
qn(estim_length*Lf*M+1:estim_length*Lf*M+2*(m-1))=0.;
[init_state,qn_dec]=slidevipostate(g,qn,Lf);
sd(:,pp)=(qn_dec(1:Lf))';
bc=encode_blockstate(g,qn_dec(1:Lf),init_state);
db(:,nn-1+pp)=conj(modmdp4(bc(1:(1/coding_rate)*Lf)))';
% On peut réestimer le canal de propagation à l'aide de ces symboles
% décidés (démodulés puis remodulés)
[C,CF]=estim_canal(db(:,pp+1:1:nn-1+pp),d(:,pp+1:1:nn-1+pp),nn-1,Lf,sigma,fd);
end;

% Pour les Lf*estim_length derniers symboles, on se sert de la dernière estimation de CF 
for kk=1:estim_length
  dn(:,kk)=CF'.*d(:,n-estim_length+kk);  
end;
  ysd=reshape(dn,estim_length*Lf,1); 
for i=1:estim_length*Lf             
qn(M*(i-1)+1:M*i)=[real(ysd(i)) imag(ysd(i))]';
end;
qn(estim_length*Lf*M+1:estim_length*Lf*M+2*(m-1))=0.;
[init_state,qn_dec]=slidevipostate(g,qn,Lf);
for jj=1:estim_length
sd(:,n-nn+1-estim_length+jj)=(qn_dec(Lf*(jj-1)+1:Lf*jj))';
end;
sx=reshape(sd,Lf*(n-training_length),1);
vx=mseq(training_length*Lf+1:n*Lf);
dsf=xor(sx,vx);
num_error(SNR(da),num_packet)=sum(dsf)/(Lf*(n-training_length));
end;
end;

    %E(i)=D(i)-zy(i,training_length+1);
    %FF(i)=h(i)-zy(i,training_length+1);
    %E(i)=abs(E(i))^2;
    %FF(i)=abs(FF(i))^2;
    %end;
    %QQ=sum(E)/Lf;
    %QR=sum(FF)/Lf;
    
TEB(1:length(SNR))=sum(num_error,2)/Nbpacket;
grid,x = 1.:1.:8.;
grid on ;
grid,semilogy(x,TEB(:),'b');
grid on ;
grid,title('BER versus Eb/N0');
    
return