function d=dist(y,z)
[n,p]=size(y);
zz=y-z;
dd=zz'*zz;
d=sqrt(dd);
return;