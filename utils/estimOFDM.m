function z = estimOFDM
%M nombre d'états de la modulation, ici une MDP4 %
n=10;
M=2;
size=16;
%taille de la FFT %
mseq = (gene_sequence(n*M*size))';
% creation de la modulation MDP4 %
for i=1:n*size
   sn=mseq(M*(i-1)+1:M*i);
   if(sn==[0 0])
      vn(i)=exp(j*pi/4);
   elseif(sn==[0 1])
      vn(i)=exp(j*(pi/4+pi/2));
   elseif(sn==[1 1])
      vn(i)=exp(j*(pi/4+pi));
   else
      vn(i)=exp(j*(pi/4+3*pi/2));
   end;
end;
% creation de la modulation multiporteuse OFDM %
for i=1:n
   w=vn(size*(i-1)+1:size*i);
   z(i,:)=sqrt(size)*ifft(w,size);
end;
% construction du canal de propagation bi-dimensionnel
% temps-espace, sélectif en temps et en fréquence 

% premier trajet 
x = canalOFDM(0.001,1,10);
% On construit les répliques fréquentielles décalées de ce trajet multiple
% On construit size colonnes chacune ayant Ns (nombre symboles) lignes 
% Ainsi, on obtient une matrice temps fréquence (colonnes => fréquences et 
% lignes => temps

for j=1:size
   zz(:,j)=x'.*exp((j-1)*2*pi/size);
end;

% Creation du signal de sortie du canal, on calcule le produit scalaire de 
% z(i,:) et zz(:,j)


return


