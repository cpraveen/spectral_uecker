% solving u_t=nu*u_xx+u-u_x u by semi-implicit time stepping via FFT 
n=32; dx=2*pi/n; x=0:dx:2*pi-dx; u0=1./cosh(4*(x-pi)); nu=0.01; h=0.01; 
ist=10; tn=40; tmax=h*(tn-1)*ist;t=0;
nv = fftshift(-n/2:1:n/2-1);mu=1./(1+nv.^2*h*nu); 
kv=zeros(1,n); % holds F-multipliers for spectral diff. 
for k=1:n/2-1; km=1i*k;kv(k+1)=km; kv(n+1-k)=-km; end
tv=0:h*ist:h*(tn-1)*ist; [X,T]=meshgrid(x,tv); [N,T]=meshgrid(-n/2:n/2-1,tv);
ua=zeros(tn,n); ua(1,:)=u0; u=u0; % ua holds u(x,t) for plotting 
uh=fft(u); uha=zeros(tn,n); uha(1,:)=abs(fftshift(uh)); % for plotting FT
for i=2:tn
    for j=1:ist 
    % fh=fft(-0.5*u.*u); % no anti-aliasing 
    fh=-0.5*aap(uh,uh); % anti-aliasing 
    uh=mu.*(uh+h*kv.*fh);u=ifft(uh);
    end
    ua(i,:)=real(u);uha(i,:)=abs(fftshift(uh));t=t+ist*h;
end
figure(1);clf;surf(X,T,ua);colormap cool;
axis([0 2*pi 0 tmax -0.1 1.1]);view(4,40); grid off;
tics('x',[0 3 6]);tics('y',[0 2 4]);tics('z',[0 1]);
figure(2);clf;surf(N,T,uha/n);colormap cool;
axis([-n/2 n/2-1 0 tmax]);view(4,40), grid off
tics('x',[-n/2+1 0 n/2-1]);tics('y',[0 2 4]);tics('z',[0 0.1 0.2]);
% long-time behaviour
more=1;more=ask('longer time? (Y/n)',more);
while more==1
for i=1:100; fh=fft(-0.5*u.*u); uh=mu.*(uh+h*kv.*fh);u=ifft(uh); end
figure(3);clf;plot(x,real(u),x,abs(fftshift(uh))/n,'-o');
t=t+100*h;title(['t= ',num2str(t)]);more=ask('longer time? (Y/n)',more);
end

