% [x_em]=em(ant1, ant2, path_1_1, path_1_2, W, Xpilot, variance, sigmaH, sigmaH_cross)
function [x_em]=em(ant1, ant2, path_1_1, path_2_1, W, Xpilot, variance, sigmaH, sigmaH_cross)
[dummy, T] = size(ant1);
totaltime = T/8;
pcount = 0;
tapstime=1;
[dummy, Lf]=size(W);
Lf=Lf/2;
clear psk4gray;
psk4gray(1)=1;
psk4gray(2)=1i;
psk4gray(3)=-1i;
psk4gray(4)=-1;
x_em_time=1;
for t=1:8:T
    h(1:Lf,1)=path_1_1(tapstime,1:Lf).';
    h(2*(1:Lf),1)=path_2_1(tapstime,1:Lf).';
    if (pcount==0)
        zp = randn(8,1)*variance+rand(8,1)*i*variance;
        yp = Xpilot*W*h+zp;
        [h_est, sigma_h_chapeau] = calc_h(yp, Xpilot, W, variance, sigmaH, sigmaH_cross);
    end
    pcount = pcount+1;
    if(pcount==20) 
        pcount=0;
    end
    Xtx = ant2mat(ant1(t:t+7), ant2(t:t+7));
    zp = randn(8,1)*variance+rand(8,1)*i*variance;
    y = Xtx*W*h+zp;
    X_kappa = mlx0(y, W, h_est, sigma_h_chapeau);
    for iter=1:3
        new_X_kappa = em_iter(y, X_kappa, W, h_est, sigma_h_chapeau, psk4gray);
        clear X_kappa;
        X_kappa=new_X_kappa;
    end
    
    for k=1:8
        x_em(x_em_time) = X_kappa(k,k);
        x_em_time=x_em_time+1;
        x_em(x_em_time) = X_kappa(k,k+8);
        x_em_time=x_em_time+1;
    end
    
end
    
    
