function fd=freqd(nb_rays,max_doppler)
for i=1:nb_rays
fd(i)=max_doppler*cos(2*pi*rand);
end

