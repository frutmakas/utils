
#define nsp_UsesSampleGen
#define nsp_UsesVector
#define nsp_UsesTransform
#define nsp_UsesConversion

#include "nsp.h"

#include "tools.h"
#include "rayleigh.h"
#include "interpol/interpol.h"

/*
function r = rayleigh( fd, fs, Ns )
% r = rayleigh(fd,fs,N)
%
% A Rayleigh fading simulator based on Clarke's Model
% Creates a Rayleigh random process with PSD determined
% by the vehicle's speed.
%
% INPUTS:
%   fd = doppler frequency
%        set fd = v*cos(theta)/lambda
%        v = velocity (meters per second)
%        lambda = carrier wavelength (meters)
%        theta = angle w/ respect to tangent (radians).
%   fs = sample frequency (Samples per second)
%   Ns = number of samples of the Rayleigh fading
%        process to produce
%
% OUTPUTS:
%   r  = row vector containing Ns samples of the Rayleigh
%        fading process 
*/

double *RayleighFading(double doppler_freq, double sample_freq, int nbSamples, int *output_size) {

/*	
N = 8;
while(N)
   if (N < 2*fd*Ns/fs) 
      N = 2*N;
   else
      break;
   end
end 
*/
	double dtemp = log2(2*doppler_freq*nbSamples/sample_freq);
	double N = (int) pow(2.0, dtemp<3.0 ?3.0 : ceil(dtemp));


/*	
% number of ifft points (for smoothing)
N_inv = ceil(N*fs/(2*fd));
*/

	double N_inv = ceil(N*sample_freq/(2*doppler_freq));
/*
% determine the frequency spacings of signal after FFT
delta_f = 2*fd/N;
*/
	double delta_f = 2*doppler_freq/N;

/*
% determine time spacing of output samples
delta_T_inv = 1/fs;
*/
	double delta_T_inv = 1/sample_freq;

/*
%%%%%%%%%%% Begin Random Input Generation %%%%%%%%%%%%
% fprintf( 'Generating Input\n');
% generate a pair of TIME DOMAIN gaussian i.i.d. r.v.'s

I_input_time = randn(1,N);
Q_input_time = randn(1,N);
*/

	double *I_input_time = nspdMalloc(N); // freed = OK 
	double *Q_input_time = nspdMalloc(N); // freed = OK 
	nspdbRandGaus(&rayleighDGausStatePtr, I_input_time, (int)N);
	nspdbRandGaus(&rayleighDGausStatePtr, Q_input_time, (int)N);

/*
% take FFT
I_input_freq = fft(I_input_time);
Q_input_freq = fft(Q_input_time);
*/
	int n_2 = (int) N/2;
	int n_21 = (int) ceil(n_2+1);
	DCplx *I_input_freq = nspzMalloc(n_21); // freed = OK 
	DCplx *Q_input_freq = nspzMalloc(n_21); // freed = OK 


	nspdRealFftNip(I_input_time, I_input_freq, (int)log2(N), NSP_Forw);
	
	nspdRealFftNip(Q_input_time, Q_input_freq, (int)log2(N), NSP_Forw);
	

/*
%%%  Generate Doppler Filter's Frequency Response %%%
% fprintf( 'Generating Doppler Filter Function\n');
% filter's DC component

  SEZ(1) = 1.5/(pi*fd);
*/
	double *SEZ = nspdMalloc(n_21); // freed = OK 
	SEZ[0] =  1.5/(NSP_PI*doppler_freq);

/*% 0 < f < fd 
for j=2:N/2
   f(j) = (j-1)*delta_f;
   SEZ(j) = 1.5/(pi*fd*sqrt(1-(f(j)/fd)^2));
   SEZ(N-j+2) = SEZ(j);
end
*/
	double *f=nspdMalloc(n_2); // freed = NO
	double doppler_pi=NSP_PI*doppler_freq; // optimisation
	double inv_doppler_freq = 1.0/doppler_freq; // optimisation
	f[0]=0; // <<<<<<<<<<<<< ????? 
	for(int i = 1; i< n_2-1; i++) {
		f[i]=i*delta_f;
		SEZ[i] = 1.5/( doppler_pi*sqrt( 1-pow( f[i]*inv_doppler_freq ,2.0 ) ) );
	}
/*
% use a polynomial fit to get the component at f = fd
% p = polyfit( f(N/2-3:N/2), SEZ(N/2-3:N/2), 3);
% SEZ(N/2+1) = polyval( p, f(N/2)+delta_f );
% k = N/2 - 1;
*/

	int k = 4; double dy=0.0;

	polint(&(f[n_2-1-k]), &(SEZ[n_2-1-k]), k, f[n_2-1]+delta_f, &(SEZ[n_2]), &dy); // peut etre pb indices !! 
	
/*	
	%%%%%%%% Perform Filtering Operation %%%%%%%%%%%%%%
	% fprintf( 'Computing Output\n' );
	% pass the input freq. components through the filter

	I_output_freq = I_input_freq .* sqrt(SEZ);
	Q_output_freq = Q_input_freq .* sqrt(SEZ);

*/
	nspdbSqrt1(SEZ, n_21);
	
	double *im_temp = nspdMalloc(n_21); // freed = OK 
	nspdbZero(im_temp, n_21);
	DCplx *DSEZ = nspzMalloc(n_21); // freed = OK 
	nspzb2RealToCplx(SEZ, im_temp, DSEZ, n_21);

	DCplx *I_output_freq = nspzMalloc(n_21); // freed = OK
	nspzbMpy3(I_input_freq, DSEZ, I_output_freq, n_21); 

	DCplx *Q_output_freq = nspzMalloc(n_21); // freed = OK
	nspzbMpy3(Q_input_freq, DSEZ, Q_output_freq, n_21);

	nspFree(im_temp);
	nspFree(SEZ);
	nspFree(DSEZ);
	nspFree(I_input_freq);
	nspFree(Q_input_freq);


// take inverse FFT

/*
I_temp = [I_output_freq(1:N/2) zeros(1,N_inv-N) I_output_freq(N/2+1:N)];
I_output_time = ifft(I_temp);
*/
	double * I_temp=nspdMalloc(N); // freed = OK
	nspdCcsFftNip(I_output_freq, I_temp, (int)log2(N), NSP_Inv);

	/*
Q_temp = [Q_output_freq(1:N/2) zeros(1,N_inv-N) Q_output_freq(N/2+1:N)];
Q_output_time = ifft(Q_temp);
*/
	double * Q_temp=nspdMalloc(N); // freed = OK
	nspdCcsFftNip(Q_output_freq, Q_temp, (int)log2(N), NSP_Inv);

	nspFree(I_output_freq);
	nspFree(Q_output_freq);


// % take magnitude squared of each component and add together
/*
	for j=1:N_inv
		r(j) = sqrt( (abs(I_output_time(j)))^2 + (abs(Q_output_time(j)))^2);
	end
*/
	nspdbSqr1(I_temp, N);
	nspdbSqr1(Q_temp, N);

	double *r = nspdMalloc(N); // freed = NO
	

	nspdbAdd3(I_temp, Q_temp, r, N);
	nspdbSqrt1(r, N);

// % normalize and compute rms level
/*	
	rms = sqrt( mean( r.*r ) );
	r = r(1:Ns)/rms;
*/
	double *r2 = nspdMalloc(N); // freed = OK
	nspdbCopy(r, r2, N);
	nspdbSqr1(r2, N);

	double rms_inv = 1/sqrt(nspdMean(r2,N));
	nspdbMpy1(rms_inv, r, N);
	nspFree(r2);
	nspFree(I_temp);
	nspFree(Q_temp);
	*output_size=N;

/*
	grid,x = 1:Ns;
	grid on;
	plot(x,r(:),'r');
	gtext('frequence doppler = 0.01')
	grid on;
*/
	return r;

}

