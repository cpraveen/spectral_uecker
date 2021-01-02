% ac.m solving u_t=u_xx+u-u^3 by semi-implicit time stepping via FFT 
n=128; dx=2*pi/n; x=0:dx:2*pi-dx; nu=0.001; h=0.1; 
ssteps=5; % ssteps=small steps, internal integration loop 
psteps=20;% psteps=plotting steps=number of time slices to plot 
tmax=h*(psteps-1)*ssteps;t=0; % holds F-multipliers 
nv = fftshift(-n/2:1:n/2-1);mu=1./(1+nv.^2*h*nu); % holds F-multipliers 
% generate IC 'blockwise random' 
nb=4; u0=[]; for k=1:n/nb; u0=[u0 0.5*(rand(1,1)-0.5)*ones(1,nb)]; end
tv=0:h*ssteps:h*(psteps-1)*ssteps; 
[X,T]=meshgrid(x,tv); [N,T]=meshgrid(-n/2:n/2-1,tv);
ua=zeros(psteps,n); ua(1,:)=u0; u=u0; uh=fft(u0);% ua holds u(x,t) for plotting 
va=zeros(psteps,n); va(1,:)=abs(fftshift(uh))/n; % for F-coefficients 
for i=2:psteps
    for j=1:ssteps 
    fh=fft(u-u.*u.*u); uh=mu.*(uh+h*fh);u=ifft(uh);
    end
    ua(i,:)=real(u);va(i,:)=abs(fftshift(uh))/n;t=t+ssteps*h;
end
figure(1);clf;surf(X,T,ua);set(gca,'FontSize',16);colormap cool;
axis([0 2*pi-dx 0 tmax -1.1 1.1]);view(4,40); grid off;
figure(2);clf;surf(N,T,va);set(gca,'FontSize',16);colormap cool;
axis([-n/2 n/2-1 0 tmax 0 0.8]);view(4,40); grid off;
%%
figure(1);clf;
waterfall(X,T,ua);colormap(1e-6*[1 1 1]);
axis([0 2*pi 0 tmax]);view(2,44);pbaspect([1 1 .33]);grid off
%%
% long-time behaviour
more=1;more=ask('longer time? (Y/n)',more);
while more==1
for i=1:100
    ff=fft(u-u.*u.*u); uf=fft(u);uf=mu.*(uf+h*ff);u=ifft(uf);
    vh=fft(v);v2h=aap(vh,vh);v3h=aap(vh,v2h);
    fh=vh-v3h; vh=mu.*(vh+h*fh);v=ifft(vh);
end
t=t+100*h;figure(2);clf;plot(x,real(u),x,real(v),'-o');
title(['t= ',num2str(t)]);more=ask('longer time? (Y/n)',more);
end

