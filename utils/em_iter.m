%[new_X_kappa] = em_iter(y, X_kappa, W, h, sigma_h_chapeau, psk_constellation)
function [new_X_kappa] = em_iter(y, X_kappa, W, h, sigma_h_chapeau, psk_constellation)

clear new_X_kappa;
for k=1:8
    qi_min =1e20;
    clear x_min subsigma temp;
    temp = W*sigma_h_chapeau*W';
    subsigma(1,1)=temp(k,k);
    subsigma(1,2)=temp(k+8,k);
    subsigma(2,1)=temp(k,k+8);
    subsigma(2,2)=temp(k+8,k+8);
    %calcul de Wf'
    clear Wfp;
    clear wf;
    [dummy, Lf2] = size(W);
    Lf=Lf2/2;
    for L=1:Lf
        wf(1, L) = exp(2*pi*i*(L-1)/8);
    end
    Wfp(1,1:Lf)=wf;
    Wfp(2,Lf+(1:Lf))=wf;
        
    for p1=1:4
        for p2=1:4
            x = [psk_constellation(p1) psk_constellation(p2)]';
            qi = abs(y(k)-x'*Wfp*h)^2 + x'*subsigma*x;
            if(qi<qi_min)
                clear x_min;
                x_min=x';
                qi_min=qi;
            end
        end
    end
    new_X_kappa(k,k)=x_min(1);
    new_X_kappa(k,k+8)=x_min(2);
end