function [O]=forros(H,x) 
% initialisation
M=40;
N=70;
for m=1:M
    z(m,:)=x;
end

%iteration start
for t=1:10
    t
    %step 1
    sigma_mn=(sign(z)+1)/2;
    for m=1:M
        sigma_m(m)=mod(sum(sigma_mn(m,:)),2);
    end
    z
    for n=1:N
        for m=1:M
            if (H(m,n)==1 )                L(m,n)=(-1)^(~(mod( sum( sigma_mn(m,[1:n-1,n+1:N]) ),2)))*min(abs(z(m,[1:n-1,n+1:N])));            end
        end
    end
    L
    %step2
    for m=1:M
        for n=1:N
            z(m,n)=x(n)+sum(L([1:m-1,m+1:M],n));
        end
    end
    zn=x(n)+sum(L(1:M,:));
    zn
    %step3
    c=(sign(zn)+1)/2;
    ctrl_chk=H*c';
    ctrl_chk
    if sum(ctrl_chk)==0
        break;
    end
end
t
O=c;
    
        
        
    
        
        
    
    
