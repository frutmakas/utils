function x=gene_sequence(n)
for i=1:10
    x(i)=1;
end;
for i=11:n
    x(i)=xor(x(i-7),x(i-3));
end;

