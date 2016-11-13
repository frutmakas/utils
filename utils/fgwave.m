function fg=fgwave(nb_rays,taps,symbol_T)
symbol_T=symbol_T*1e-6;
Tmax=10e-6;
b_tau=1e-6;
for i=1:nb_rays
    tau(i)=-b_tau*log10(1-rand*(1-exp(-Tmax/b_tau)));
    for k=1:taps
        tap_factor=((k-1)*symbol_T-tau(i))/symbol_T;
        if (abs(tap_factor)<1)
            fg(i,k)=1-abs(tap_factor);
        else
            fg(i,k)=0;
        end
    end
end