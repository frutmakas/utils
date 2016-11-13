function XX = canalOFDM( fd, fs, Ns )
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

N = 8;
while(N)
   if (N < 2*fd*Ns/fs) 
      N = 2*N;
   else
      break;
   end
end 

% number of ifft points (for smoothing)
N_inv = ceil(N*fs/(2*fd));

% determine the frequency spacings of signal after FFT
delta_f = 2*fd/N;

% determine time spacing of output samples
delta_T_inv = 1/fs;

%%%%%%%%%%% Begin Random Input Generation %%%%%%%%%%%%
% fprintf( 'Generating Input\n');
% generate a pair of TIME DOMAIN gaussian i.i.d. r.v.'s
I_input_time = randn(1,N);
Q_input_time = randn(1,N);
R_input_time = randn(1,N);
S_input_time = randn(1,N);
T_input_time = randn(1,N);
U_input_time = randn(1,N);
V_input_time = randn(1,N);
W_input_time = randn(1,N);
X_input_time = randn(1,N);
Y_input_time = randn(1,N);
Z_input_time = randn(1,N);



% take FFT
I_input_freq = fft(I_input_time);
Q_input_freq = fft(Q_input_time);
R_input_freq = fft(R_input_time);
S_input_freq = fft(S_input_time);
T_input_freq = fft(T_input_time);
U_input_freq = fft(U_input_time);
V_input_freq = fft(V_input_time);
W_input_freq = fft(W_input_time);
X_input_freq = fft(X_input_time);
Y_input_freq = fft(Y_input_time);
Z_input_freq = fft(Z_input_time);


%%%  Generate Doppler Filter's Frequency Response %%%
% fprintf( 'Generating Doppler Filter Function\n');
% filter's DC component
SEZ(1) = 1.5/(pi*fd);

% 0 < f < fd 
for j=2:N/2
   f(j) = (j-1)*delta_f;
   SEZ(j) = 1.5/(pi*fd*sqrt(1-(f(j)/fd)^2));
   SEZ(N-j+2) = SEZ(j);
end

% use a polynomial fit to get the component at f = fd
% p = polyfit( f(N/2-3:N/2), SEZ(N/2-3:N/2), 3);
% SEZ(N/2+1) = polyval( p, f(N/2)+delta_f );
% k = N/2 - 1;

k = 3;
p = polyfit( f(N/2-k:N/2), SEZ(N/2-k:N/2), k );
SEZ(N/2+1) = polyval( p, f(N/2)+delta_f );

%%%%%%%% Perform Filtering Operation %%%%%%%%%%%%%%
% fprintf( 'Computing Output\n' );
% pass the input freq. components through the filter
I_output_freq = I_input_freq .* sqrt(SEZ);
Q_output_freq = Q_input_freq .* sqrt(SEZ);
R_output_freq = R_input_freq .* sqrt(SEZ);
S_output_freq = S_input_freq .* sqrt(SEZ);
T_output_freq = T_input_freq .* sqrt(SEZ);
U_output_freq = U_input_freq .* sqrt(SEZ);
V_output_freq = V_input_freq .* sqrt(SEZ);
W_output_freq = W_input_freq .* sqrt(SEZ);
X_output_freq = X_input_freq .* sqrt(SEZ);
Y_output_freq = Y_input_freq .* sqrt(SEZ);
Z_output_freq = Z_input_freq .* sqrt(SEZ);


% take inverse FFT
I_temp = [I_output_freq(1:N/2) zeros(1,N_inv-N) I_output_freq(N/2+1:N)];
I_output_time = ifft(I_temp);

Q_temp = [Q_output_freq(1:N/2) zeros(1,N_inv-N) Q_output_freq(N/2+1:N)];
Q_output_time = ifft(Q_temp);

R_temp = [R_output_freq(1:N/2) zeros(1,N_inv-N) R_output_freq(N/2+1:N)];
R_output_time = ifft(R_temp);

S_temp = [S_output_freq(1:N/2) zeros(1,N_inv-N) S_output_freq(N/2+1:N)];
S_output_time = ifft(S_temp);

T_temp = [T_output_freq(1:N/2) zeros(1,N_inv-N) T_output_freq(N/2+1:N)];
T_output_time = ifft(T_temp);

U_temp = [U_output_freq(1:N/2) zeros(1,N_inv-N) U_output_freq(N/2+1:N)];
U_output_time = ifft(U_temp);

V_temp = [V_output_freq(1:N/2) zeros(1,N_inv-N) V_output_freq(N/2+1:N)];
V_output_time = ifft(V_temp);

W_temp = [W_output_freq(1:N/2) zeros(1,N_inv-N) W_output_freq(N/2+1:N)];
W_output_time = ifft(W_temp);

X_temp = [X_output_freq(1:N/2) zeros(1,N_inv-N) X_output_freq(N/2+1:N)];
X_output_time = ifft(X_temp);

Y_temp = [Y_output_freq(1:N/2) zeros(1,N_inv-N) Y_output_freq(N/2+1:N)];
Y_output_time = ifft(Y_temp);

Z_temp = [Z_output_freq(1:N/2) zeros(1,N_inv-N) Z_output_freq(N/2+1:N)];
Z_output_time = ifft(Z_temp);


% take magnitude squared of each component and add together
for j=1:Ns
   r(j) = sqrt((abs(I_output_time(j)))^2+(abs(Q_output_time(j)))^2+(abs(R_output_time(j)))^2+...
      (abs(S_output_time(j)))^2+(abs(T_output_time(j)))^2+(abs(U_output_time(j)))^2+(abs(V_output_time(j)))^2+...
      (abs(W_output_time(j)))^2+(abs(X_output_time(j)))^2+(abs(Y_output_time(j)))^2+(abs(Z_output_time(j)))^2);
end

% normalize and compute rms level
rms = sqrt(mean( r(1:Ns).*r(1:Ns)));
XX = (I_output_time(1:Ns)+Q_output_time(1:Ns)+R_output_time(1:Ns)+...
   S_output_time(1:Ns)+T_output_time(1:Ns)+U_output_time(1:Ns)+...
   V_output_time(1:Ns)+W_output_time(1:Ns)+X_output_time(1:Ns)+...
   Y_output_time(1:Ns)+Z_output_time(1:Ns))/rms;

return