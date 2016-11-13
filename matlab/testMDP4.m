function  TEB = testMDP4
SNR = [1,2,3,4,5];
num_error(1:length(SNR))=0;
for q=1:length(SNR)
n=100000;   %nombre d'échantillons temporels
M=2;
feet=n*M;
seq=(rand(1,feet)>0.5)';

% creation de la modulation MDP4 %
for i=1:n
   sn=seq(M*(i-1)+1:M*i);
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
    [d,std]=awgn(1,n,SNR(q),vn);
        
for k=1:n
          if ((real(d(k))>0) & (imag(d(k))>0))
             db(k)=exp(j*pi/4);
         elseif((real(d(k))>0) & (imag(d(k))<0))
                db(k)=exp(j*(pi/4+3*pi/2));
            elseif((real(d(k))<0) & (imag(d(k))>0))
                db(k)=exp(j*(pi/4+pi/2));
            else 
                db(k)=exp(j*(pi/4+pi));
            end;
        end;
        
for ii=1:n
   seq_dec(M*(ii-1)+1:M*ii)=[0.5*(sign(real(db(ii)))+1);0.5*(sign(imag(db(ii)))+1)]';
end;
dsf=xor(seq_dec,seq');
num_error(q)=sum(dsf)/feet;
end;

num_error
TEB(1:length(SNR))=num_error(1:length(SNR));
grid,x = 1.:1.:5.;
grid on ;
grid,semilogy(x,TEB(:),'b');
grid on ;
grid,title('BER versus Eb/N0');
    
return
    
    
    
    