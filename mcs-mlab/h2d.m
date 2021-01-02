% h2d.m, solving u_t=Delta u by FFT 
n=32; dx=2*pi/n; dy=dx; x=0:dx:2*pi-dx; y=x; [X Y]=meshgrid(x,y); 
eta=1; h=0.1; mu=ones(n,n); % holds F-multipliers 
for j1=1:n; k1=j1-n/2-1; % fill multiplier matrix mu as if fft(u) was centered 
    for j2=1:n; k2=j2-n/2-1;ks=k1*k1+k2*k2; 
        mu(j1,j2)=(1-(1-eta)*ks^2*h)/(1+eta*ks^2*h);
    end
end
mu=fftshift(mu); % now adapt mu to true fft
u0=1+sin(X)+sin(2*Y); umin=min(min(u0)); umax=max(max(u0));
clf;p1=surf(X,Y,u0,'EraseMode','background');
set(gca,'FontSize',16);axis([0 2*pi 0 2*pi umin umax]);
more=1;u=u0; uf=fft2(u); more=askmore(more);
while more==1 % integration loop
  uf=mu.*uf;  % actual integration, next line for plotting
  u=ifft2(uf); set(p1,'zdata',u);drawnow; more=askmore(more);
end

