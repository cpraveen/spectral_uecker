% ac2d.m, solving u_t=Delta u+u-u^3 by FFT 
n=128; dx=2*pi/n; dy=dx; x=0:dx:2*pi-dx; y=x; [X Y]=meshgrid(x,y); 
nu=0.001; h=0.1; mu=ones(n,n); % holds F-multipliers 
for j1=1:n; k1=j1-n/2-1; % fill multiplier matrix mu as if fft(u) was centered 
    for j2=1:n; k2=j2-n/2-1;ks=k1*k1+k2*k2; mu(j1,j2)=1/(1+nu*ks*h);end
end
mu=fftshift(mu); % now adapt mu to true fft
u0=zeros(n,n); nb=8; % generate IC 'blockwise random', blocksize nb
for j1=0:n/nb-1 
    for j2=0:n/nb-1
        u0(j1*nb+1:(j1+1)*nb, j2*nb+1:(j2+1)*nb)=0.5*(rand(1,1)-0.5)*ones(nb,nb); 
    end
end
figure(1);clf; colormap cool; p1=surf(X,Y,u0);
set(gca,'FontSize',20); shading interp;
axis([0 2*pi 0 2*pi -1.1 1.1]);view(-10,70);grid off; 
tits=['t=0 '];title(tits);
%pcolor(x,y,u0);colormap Gray; colorbar;
set(gca,'FontSize',16);colormap Gray; shading interp;
more=1;u=u0; t=0; uf=fft2(u); more=askmore(more);
while more==1 % integration loop
  uf=mu.*(uf+h*fft2(u-u.^3));  u=ifft2(uf); t=t+h;
  %clf; pcolor(x,y,u);shading interp; colorbar;
  clf; p1=surf(X,Y,u);
  set(gca,'FontSize',20);colormap cool; shading interp;
  axis([0 2*pi 0 2*pi -1.1 1.1]);view(-10,70);grid off; 
  tits=['t=', num2str(t)];title(tits);
  more=askmore(more);
end

