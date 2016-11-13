function [init_state,x_hat] = slidevipostate( g, r, Lf)
% Usage: x_hat = viterbi( g, r )
%
% Soft decision Viterbi decoding algorithm.
%
% Takes an entire convolutionally encoded (and possibly corrupted)
% codeword with values between [-1 and 1] and returns an estimate
% of the information bits (not including tail bits).
% 15-06-98 CANCES Jean-Pierre
[n,K] = size(g);
m = K - 1;
max_state = 2^m;
[temp, rec_size] = size(r);
L_total = rec_size/n;
L_info = L_total - m;
length_depth=5*Lf;

% set infiinity to an arbitrarily large value
inf = 10^5;

% initialize trellis and path matrices
trellis = inf*ones(max_state,L_total);
trellis(1:max_state,1)= 0.0;
path = zeros(max_state,L_total);
new_path = path;

% initialize output and transition matrices
for state=1:max_state
   state_vector = bin_state( state-1, m );
   [out_0, state_0] = encode_bit(g, 0, state_vector);
   [out_1, state_1] = encode_bit(g, 1, state_vector);
   output(state,:) = [out_0 out_1];
   transition(state,:) = [(int_state(state_0)+1) (int_state(state_1)+1)];
end

% now determine trellis and path matrices for times 1 through length_depth
counter = 1;
for time=1:length_depth
   y_segment = r(1,counter:counter+n-1);
   counter = counter + n;
   for state=1:max_state
      for i=0:1
         hypothesis = 2*output( state, n*i+1:n*(i+1) )-1;
         next_state = transition( state, i+1 );

         % compute squared Euclidian distance
         square_dist = 0;
         for j = 1:n
            % y_segment(j) == 0 indicates an erasure
            if y_segment(j) 
               square_dist = square_dist + (hypothesis(j)-y_segment(j))^2;
            end
         end

         % update path metric
         path_metric = square_dist + trellis(state,time);
         if path_metric < trellis( next_state, time+1 )
            trellis( next_state, time+1 ) = path_metric;
            new_path( next_state, 1:time+1 ) = [path(state,1:time) i];
            antecedent(i+1,next_state)=state;
         end
      end
   end
   path = new_path;
   xx(time,:)=[antecedent(1,1:max_state/2),antecedent(2,max_state/2+1:max_state)];
end

%selection of the best surviving paths among the remaining num_states ones 
trellis_min=trellis(1,length_depth+1); jj=1.;
for i=1:max_state
    if trellis(i,length_depth+1)<trellis_min
        trellis_min=trellis(i,length_depth+1);
        jj=i;
    end;
end;
state_fin(length_depth+1)=jj;

%trace_back for determining initial state and symbol
x_hat(1)=path(jj,2);
etat(length_depth+1)=jj;
for ii=1:length_depth
   etat(length_depth+1-ii)=xx(length_depth+1-ii,etat(length_depth+2-ii));
end;
    init_state = etat(1);               
    
% now determine trellis and path matrices for times length_depth+1 through L_info 
for time=length_depth+1:L_info
y_segment = r(1,counter:counter+n-1);
   counter = counter + n;
   for state=1:max_state
      for i=0:1
         hypothesis = 2*output( state, n*i+1:n*(i+1) )-1;
         next_state = transition( state, i+1 );

         % compute squared Euclidian distance
         square_dist = 0;
         for j = 1:n
            % y_segment(j) == 0 indicates an erasure
            if y_segment(j) 
               square_dist = square_dist + (hypothesis(j)-y_segment(j))^2;
            end
         end

         % update path metric
         path_metric = square_dist + trellis(state,time);
         if path_metric < trellis( next_state, time+1 )
            trellis( next_state, time+1 ) = path_metric;
            new_path( next_state, 1:time+1 ) = [path(state,1:time) i];
            %antecedent(i+1,next_state)=state;
         end
      end
   end
   path = new_path;
   %selection of the best surviving paths among the remaining num_states ones 
   trellis_min=trellis(1,time+1); jj=1.;
   for i=1:max_state
    if trellis(i,time+1)<trellis_min
        trellis_min=trellis(i,time+1);
        jj=i;
    end;
end;
x_hat(time-length_depth+1) = path(jj,time-length_depth+2);
state_fin(time)=jj;
end;
x_hat(L_info-length_depth+2:L_info) = path(jj,L_info-length_depth+3:L_info+1);

return

