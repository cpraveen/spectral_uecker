function [u1,u2]=ns2dstep(u1,u2,n,h,t); %  used in ns2d.m
% perform dt=h time step, u1,u2 in x-space (in and out), no anti-aliasing
global X Y av bv mm w1 w2 w3 w4 g1h g2h; 
v1h=fft2(u1.*u1); v2h=fft2(u1.*u2);v3h=fft2(u2.*u2);
[g1,g2]=g(X,Y,t); g1h=fft2(g1); g2h=fft2(g2); % comment out for stat. forcing
for j=1:n % differentiate in Fourier 
    w1(j,:)=av.*v1h(j,:); w2(:,j)=bv.*v2h(:,j);
    w3(j,:)=av.*v2h(j,:); w4(:,j)=bv.*v3h(:,j);
end
r1=-w1-w2+g1h; r2=-w3-w4+g2h; 
u1=mm.*(fft2(u1)+h*r1);u2=mm.*(fft2(u2)+h*r2);
[u1,u2]=proj2d1(u1,u2); u1=ifft2(u1); u2=ifft2(u2);
end
