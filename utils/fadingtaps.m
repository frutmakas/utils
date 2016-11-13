%Generating Fading Taps, August 1997
%Tao Wu Version 1.1 - Date: August 29, 1997
%Version 2.1 - Date: Nov 1998
%
% INPUT:
% max_doppler: Maximum doppler frequency
% symbol_time: The symbol time in micro seconds
% filter_length: The filter tap number
% nb_samples: The number of samples to generate
% nb_rays: The number of simulated rays
%
% INTERIM:
% fd:frequency spectrum
% g:signal waveform
%
% OUTPUT:
% The output is a matrix of complex sample values for all the taps.
% The kth column is the kth tap, where k=1,2....
function taps=fadingtaps(max_doppler,symbol_time,filter_length,nb_samples,nb_rays)
fd=freqd(nb_rays,max_doppler);
g=fgwave(nb_rays,filter_length,symbol_time);
%tap2=fadingmain(symbol_time,filter_length,nb_samples,nb_rays,1000*rand,fd,g)';
taps=tapsgen(symbol_time, filter_length, nb_samples, nb_rays, fd, g,  1000*rand);
%m=1:filter_length;k=1:nb_samples;
%taps(k,m)=tap2(k,2*m-1)+i*tap2(k,2*m);

