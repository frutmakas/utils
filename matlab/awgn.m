function [z,d] = awgn(coding_rate,n,SNR,x)
if nargin ~=4
   error('Nombre d''arguments d''appel incorrect');
end;
varz = std(x).^2;
bruit = varz*0.5*10^(-SNR/10)*(1/coding_rate);
b = randn(1,n)+j*randn(1,n);
moy = mean(b);
bc = b - moy;
varb = std(bc).^2;
dbr = sqrt((bruit/varb))*bc;
z = x(:) + dbr(:);
d = std(dbr).^2;