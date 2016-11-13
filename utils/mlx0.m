% [X_kappa] = mlx0(y, W, h, sigma_h_chapeau)
function [X_kappa] = mlx0(y, W, h, sigma_h_chapeau)

psk4gray(1)=1;
psk4gray(2)=1i;
psk4gray(3)=-1i;
psk4gray(4)=-1;
    
for k=1:8
    [dummy, Lf2] = size(W);
    Lf=Lf2/2;
    for L=1:Lf
        wf(1, L) = exp(2*pi*i*(L-1)/8);
    end
    Wfp(1,1:Lf)=wf;
    Wfp(2,Lf+(1:Lf))=wf;
    qi_min=1e20; 
    k;
   for p1=1:4
        for p2=1:4
            clear xtest;
            xtest(1,1)=psk4gray(p1);
            xtest(1,2)=psk4gray(p2);
            xtest=xtest';
            temp = W*sigma_h_chapeau*W';
            sigma_h_chapeau_k = [[temp(k,k) temp(8+k, k)];[temp(k,k+8) temp(8+k, k+8)]];
            qi=abs(y(k)-xtest'*Wfp*h)^2; %+xtest'*sigma_h_chapeau_k*xtest
            if(qi_min > qi) 
                qi_min=qi;
                xv = xtest';
            end
        end
    end
    X_kappa(k,k)=xv(1);
    X_kappa(k,8+k)=xv(2);
end
    
            
    