function [u1,u2]=ns2dstepaa(u1,u2,n,h,t); % used in ns2d.m
% perform dt=h time step, u1,u2 in x-space (in and out) 
global X Y av bv mm w1 w2 w3 w4 g1h g2h; 
u1h=fft2(u1); u2h=fft2(u2); 
v1h=aap2(u1h,u1h); v2h=aap2(u1h,u2h);v3h=aap2(u2h,u2h); % nonlinearity 
[g1,g2]=g(X,Y,t); g1h=fft2(g1); g2h=fft2(g2); % comment out for stat. forcing
for j=1:n                  % differentiate nonlin. in Fourier 
    w1(j,:)=av.*v1h(j,:); % loop over x=column index
    w2(:,j)=bv.*v2h(:,j); % loop over y=row index
    w3(j,:)=av.*v2h(j,:); w4(:,j)=bv.*v3h(:,j);
end
r1=-w1-w2+g1h; r2=-w3-w4+g2h; 
u1=mm.*(u1h+h*r1);u2=mm.*(u2h+h*r2); % integrate 
[u1,u2]=proj2d1(u1,u2);              % project
u1=ifft2(u1); u2=ifft2(u2);          % return to x-space 
end
