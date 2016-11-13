function int_state = int_state( state )
% Usage: int_state = int_state( state )
%
% Converts a row vector 'state' of m bits into a integer 
% 'int_state' (which is base 10).

[dummy, m] = size( state );

for i = 1:m
   vect(i) = 2^(m-i);
end

int_state = state*vect';


