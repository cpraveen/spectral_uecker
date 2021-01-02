function [g1f,g2f]=lowp2d(g1,g2,k) % used in ns2d.m 
% trivial lowpass in rectangle with width k/2
n=length(g1); nc=n/2+1;g1f=zeros(n,n);g2f=g1f;
g1=fftshift(fft2(g1)); g2=fftshift(fft2(g2));
g1f(nc-k:nc+k, nc-k:nc+k)=g1(nc-k:nc+k, nc-k:nc+k);
g2f(nc-k:nc+k, nc-k:nc+k)=g2(nc-k:nc+k, nc-k:nc+k);
g1f=ifft2(ifftshift(g1f));g2f=ifft2(ifftshift(g2f));
end