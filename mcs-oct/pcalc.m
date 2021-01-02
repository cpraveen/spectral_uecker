function p=pcalc(u1,u2,n,t);
global X Y av bv p4 p5 w1 w2 w3 w4; 
v1h=fft2(u1.*u1); v2h=fft2(u1.*u2); v3h=fft2(u2.*u2);
for j=1:n % diff in Fourier 
    w1(j,:)=av.*v1h(j,:); w2(:,j)=bv.*v2h(:,j);
    w3(j,:)=av.*v2h(j,:); w4(:,j)=bv.*v3h(:,j);
end
[g1,g2]=g(X,Y,t); g1h=fft2(g1); g2h=fft2(g2); 
ph=p4.*(-w1-w2+g1h)+p5.*(-w3-w4+g2h); p=ifft2(ph);
end