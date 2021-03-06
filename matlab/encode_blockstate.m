function y = encode_blockstate(g, x, state)
% Usage: y = encode_block(g, x)
%
% This function takes as input an entire block of information bits 'x'
% (which are arranged in a row vector), and the coeficients of the 
% generator polynomials 'g', and returns as output an entire 
% convolutionally encoded codeword 'y'.  Tail bits are automatically 
% appended to force the encoder back to the all-zeros state.
%
% The encoder matrix 'g' has n rows and K columns, where the code rate
% is 1/n and the constraint length is K.  Each row corresponds to one
% of the n encoder polynomials and each column of the row is a 
% coefficient of the encoder polynomial.
% determine the constraint length (K), memory (m), and rate (1/n)
% and number of information bits.
[n,K] = size(g);
m = K - 1;
[temp,L_info] = size(x);
% initialize the state vector
state = bin_state(state-1,m);

% zero pad the codeword
L_total = L_info+m;
x(1,1:L_total) = [x(1:L_info) zeros(1,m)];
% generate the codeword
for i = 1:L_total
   input_bit = x(1,i);
   [output_bits, state] = encode_bit(g, input_bit, state);
   y(n*(i-1)+1:n*i) = output_bits;
end

