% checks of (anti)aliasing (Fig.4 in [U09MCS])
n=8; dx=2*pi/n; x=0:dx:2*pi-dx; xfine=0:dx/50:2*pi;
u=cos(3*x); clf;plot(x,u,'--ko',xfine,cos(3*xfine));
axis([0 2*pi -1.1 1.7]);
h=legend('u (discrete)','cos(3x)',2);
break;  % work around matlab cell mode by 'copy paste by hand' 
uh=fft(u);up=ifft(aap(uh,uh)); 
plot(x,u.*u,'--ko',x,up,'--ks'); hold on; plot(xfine,0.5+0.5*cos(6*xfine));
axis([0 2*pi -0.1 1.5]);
h=legend('u^2 (discrete)', 'u^2 (anti-alias)','0.5*(1+cos(6x))',2);
axis([0 2*pi -0.1 1.5]);


