 
function [C,CF]=estim_canal(C,D,training_length,Lf,sigma,fd);

for k=1:training_length
     yob((k-1)*Lf+1:k*Lf,1)=C(:,k);
     A((k-1)*Lf+1:k*Lf,:)=diag(D(:,k));
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
D=U*C*V'*h;
for i=1:Lf 
    CF(i)=conj(D(i))/abs(D(i));
end;