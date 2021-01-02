% ks.m solving KS by semi-implicit time stepping via FFT 
n=128; dx=2*pi/n; x=0:dx:2*pi-dx; h=dx/20;  sf=1/20; 
% sf=scale-factor due to rescaling [0 2pi/sf]->[0 2pi]
u0=1./cosh(4*(x-2))+0.5./cosh(2*(x-4));
ist=0.1/h; tn=50; tmax=h*(tn-1)*ist;t=0; % ist=internal steps, tn=plot-steps
mu=zeros(1,n); kv=zeros(1,n); % holds F-multipliers 
for k=1:n/2-1; mp=1/(1+h*(sf^3*k^4-sf*k^2)); km=1i*k;
    mu(k+1)=mp; mu(n+1-k)=mp; kv(k+1)=km; kv(n+1-k)=-km;
end
mu(1)=1;mu(n/2+1)=1/(1+h*(sf^3*(n/2)^4-sf*(n/2)^2));
tv=0:h*ist:h*(tn-1)*ist; [X,T]=meshgrid(x,tv); [N,T]=meshgrid(-n/2:n/2-1,tv);
ua=zeros(tn,n); ua(1,:)=u0; u=u0; % ua holds u(x,t) for plotting 
uh=fft(u); uha=zeros(tn,n); uha(1,:)=abs(fftshift(uh))/n; % for plotting FT
for i=2:tn
    for j=1:ist 
    %fh=fft(-0.5*u.*u);
    fh=-0.5*aap(uh,uh); % anti-aliasing 
    uh=mu.*(uh+h*kv.*fh);u=ifft(uh);
    end
    ua(i,:)=real(u);uha(i,:)=abs(fftshift(uh))/n;t=t+ist*h;
end
%%
figure(1);clf;
waterfall(X,T,ua);colormap(1e-6*[1 1 1]);set(gca,'FontSize',16)
axis([0 2*pi 0 tmax]);view(2,44);pbaspect([1 1 .33]);grid off
%% the same over original scales
figure(1);clf;
waterfall(X/sf,T/sf,ua);colormap(1e-6*[1 1 1]);set(gca,'FontSize',16)
axis([0 2*pi/sf 0 tmax/sf]);view(2,44);pbaspect([1 1 .33]);grid off

%%
figure(1);clf;colormap Gray;pcolor(x,tv,ua); colorbar;
 shading interp;axis([0 2*pi 0 tmax]); 
%%
kplot=40;[N,T]=meshgrid(0:kplot,tv);
figure(2);clf;colormap Gray;
pcolor(N,T,uha(:,n/2+1:n/2+kplot+1)); colorbar;
axis([0 kplot 0 tmax]);shading interp;
%%
kplot=50;[N,T]=meshgrid(0:kplot,tv);
figure(2);clf;colormap Gray;
waterfall(N,T,uha(:,n/2+1:n/2+kplot+1));colormap(1e-6*[1 1 1]);
axis([0 kplot 0 tmax]);view(2,44);pbaspect([1 1 .33]);grid off
%%
figure(2);clf;surf(N,T,uha/n);set(gca,'FontSize',16);colormap cool;
axis([-n/2 n/2+1 0 tmax]);view(4,40)
% long-time behaviour
more=1;more=ask('longer time? (Y/n)',more);
while more==1
for i=1:100; fh=fft(-0.5*u.*u); uh=mu.*(uh+h*kv.*fh);u=ifft(uh); end
figure(3);clf;plot(x,real(u),x,abs(fftshift(uh))/n,'-o');
t=t+100*h;title(['t= ',num2str(t)]);more=ask('longer time? (Y/n)',more);
end
%%
s1=1; s2=1/20;
k=0:22; kf=0:0.1:22;
plot(k,s2*k.^2-s2^3*k.^4, '-ko'); set(gca,'FontSize',20);axis([0 22 -5 5]);
%%
s=1; k=0:5; 
plot(k,s*k.^2-s^3*k.^4, '-ko'); set(gca,'FontSize',20);axis([0 5 -600 10]);
