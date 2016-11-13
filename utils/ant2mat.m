% X]=ant2mat(ant1, ant2)
function [X]=ant2mat(ant1, ant2)

[m,n]=size(ant1);
for k=1:n
    X(k,k)=ant1(k);
    X(k,n+k)=ant2(k);
end
