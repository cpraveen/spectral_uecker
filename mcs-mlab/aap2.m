function ph=aap2(uh,vh) %anti-aliased product, 2D 
% in: uh,vh from fft with n samples
n=length(uh); m=n*3/2;uhp=zeros(m,m);vhp=zeros(m,m);ph=zeros(n,n); 
for j=1:n/2
    uhp(j,1:n/2)=uh(j,1:n/2); uhp(j,m-n/2+1:m)=uh(j,n/2+1:n); 
    vhp(j,1:n/2)=vh(j,1:n/2); vhp(j,m-n/2+1:m)=vh(j,n/2+1:n); 
end
for j=1:n/2 
    uhp(m-n/2+j,1:n/2)=uh(n/2+j,1:n/2);vhp(m-n/2+j,1:n/2)=vh(n/2+j,1:n/2); 
    uhp(m-n/2+j,m-n/2+1:m)=uh(n/2+j,n/2+1:n); vhp(m-n/2+j,m-n/2+1:m)=vh(n/2+j,n/2+1:n); 
end
%uh, uhp
up=ifft2(uhp); vp=ifft2(vhp); w=up.*vp; wh=fft2(w);
for j=1:n/2
    ph(j,:)=1.5*[wh(j,1:n/2) wh(j,m-n/2+1:m)]; 
    ph(n/2+j,:)=1.5*[wh(m-n/2+j,1:n/2) wh(m-n/2+j,m-n/2+1:m)]; 
end;
end