function ar=red2d(a,n,nr) % resize matrix a using only 
% each nr-th point
ar=zeros(n/nr,n/nr);
for k1=1:n/nr
    for k2=1:n/nr
    l1=(k1-1)*nr+1; l2=(k2-1)*nr+1; ar(k1,k2)=a(l1,l2); 
    end
end
end