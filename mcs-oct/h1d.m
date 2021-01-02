% solving u_t=u_xx by FFT and impl. (eta>0) time stepping
% Note:  if c_k is the k-th DFT-coeff then 
%         fft(u)=(c_0 c_1...c_(n/2-1)  c_(-n/2)  c_(-n/2+1)...c_(-1))
% i.e.: position   1   2      n/2        n/2+1    n/2+2        n
n=50; dx=2*pi/n; x=0:dx:2*pi-dx; eta=0; h=0.3; tn=11; tmax=h*(tn-1);
nv=fftshift(-n/2:1:n/2-1); % wave numbers 
mu=(1-h*(1-eta)*nv.^2)./(1+eta*nv.^2*h); % F-multiplier vector 
u0=2+sin(x)+sin(2*x); t=0:h:h*(tn-1); [X,T]=meshgrid(x,t); 
ua=zeros(tn,n); ua(1,:)=u0; u=u0; uf=fft(u);
for i=2:tn
  uf=mu.*uf; u=ifft(uf);ua(i,:)=u; % back to x only for plotting
end
surf(X,T,real(ua)); view(10,50);axis([0 2*pi 0 tmax 0 4]); colormap cool;
grid off; %tics('x',[0 3 6]);tics('y',[0 0.5 1]);tics('z',[0 2 4]);
%print -F:20 h1d.eps


