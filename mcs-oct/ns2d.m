% solving Navier--Stokes in a 2D periodic box via FFT, 
% including visualization of pressure and power-spectra (if desired). 
% For some background see nse-fft.pdf 
% at www.staff.uni-oldenburg.de/hannes.uecker/hupublic.html
% or 'A short ad hoc introduction to spectral methods for parabolic PDE 
%     and the Navier--Stokes equations'
%    In Summer school MCS 09, Universität Oldenburg
% 
% Comment out undesired (plot)commands, e.g., plot of p, plot of g etc.
% Definition of forcing in g.m
% Time-stepping in ns2dstepaa.m, using an anti-aliased product aap2. 
% Auxiliary functions: vof.m, asknu.m
% Some global variables defined for convenience.
% 
% For questions/comments mail to hannes.uecker@uni-oldenburg.de
% -----------------------------------------------------------------------
rey=10; nl=max(ceil(log2(rey^0.75))+1,6); n=2^nl; n=32;
rey=asknu('Rey',rey);n=asknu('n',n);dx=2*pi/n, dy=dx; 
np=2; np=asknu('Use each np-th point in quiver plots', np);  
tp=1;tp=asknu('dT (time-interval between plots)',tp); t=0; 
h=dx/5; h=asknu('h',h); intsteps=ceil(tp/h); h=tp/intsteps; 
exsteps=5; exsteps=asknu(['number of steps of length ' num2str(tp) ' before input'],exsteps); 
global X Y av bv mm w1 w2 w3 w4 p1 p2 p3 p4 p5 g1h g2h; 
x=0:dx:2*pi-dx; y=x; [X Y]=meshgrid(x,y); more=1;
kv=[0:(n/2-1) (-n/2):-1];[KX KY]=meshgrid(kv,kv); 
kmod=sqrt(KX.^2 + KY.^2)+1e-10; % matrix of |k| (for  plotting power spectra)
xr=red2d(X,n,np); yr=red2d(Y,n,np); t=0; % for quiver 
av=zeros(1,n); % av=dx, hence belongs to k2, hence row-vector
bv=zeros(n,1); % bv=dy, hence belongs to k1, hence column-vect. 
mm=zeros(n,n); w1=zeros(n,n);w2=zeros(n,n);w3=zeros(n,n);w4=zeros(n,n); 
p1=zeros(n,n); p2=p1; p3=p1;p4=p1; p5=p1; p=zeros(n,n);
makemult(n,rey,h); % generate matrices
[g1,g2]=g(X,Y,t); g1h=fft2(g1); g2h=fft2(g2); % store g1h,g2h here once 
% if forcing is stationary, and comment out call of g in ns2dstep 
% ------------------ plot forcing for checking ---------------------
g1r=red2d(g1,n,np); g2r=red2d(g2,n,np); figure(1);
clf;quiver(xr,yr,g1r,g2r,2,'r'); 
axis([0 2*pi 0 2*pi]); ask('forcing (return)',more);
% ------------ show projection of rhs and pressure   ---------------
%[gp1,gp2,ph]=proj2d2(g1h,g2h);gp1=ifft2(gp1); gp2=ifft2(gp2); p=ifft2(ph); 
%p=p+g1h(1,1).*X/n+g2h(1,1).*Y/n; % set linear pressure for zero mean flow
%gp1=red2d(gp1,n,np);gp2=red2d(gp2,n,np); 
%hold on;quiver(xr,yr,real(gp1),real(gp2),3,'k');
%more=ask('rhs post-proj (return to continue)',more);
%figure(2);clf;surf(X,Y,real(p));  axis([0 2*pi 0 2*pi]);
%set(gca,'FontSize',16);colormap gray;view(-10,60);
%more=ask('pressure part of rhs (return to continue)',more);
% -------------------------------- actual start of program 
u1=0.*X; u2=0.*X; t=0; % IC
% project IC and calculate initial pressure (only needed if u1,u2 not zero) 
%u1h=fft2(u1); u2h=fft2(u2);[u1h,u2h]=proj2d1(u1h,u2h); 
%u1=ifft2(u1h); u2=ifft2(u2h);
%quiver(X,Y,u1,u2); ask('IC (return)',more); % visualize IC 
while more==1 % main loop ------------------------------------------
    cfl=[];   % to give user feedback on CFL 
    for j1=1:exsteps % steps till next user input 
         for j2=1:intsteps % inner integration steps  
      %  [u1,u2]=ns2dstep(u1,u2,n,h,t);t=t+h; % (no anti-alias)
         [u1,u2]=ns2dstepaa(u1,u2,n,h,t);t=t+h; % anti-alias version 
        end % inner int. loop done, now plot 
% ---- quiver reduced flow field
    xr=red2d(X,n,np); yr=red2d(Y,n,np);  
    u1r=red2d(u1,n,np); u2r=red2d(u2,n,np); 
    figure(1);clf;quiver(xr,yr,real(u1r),real(u2r),0.8); 
    axis([0 2*pi 0 2*pi]);
    ua=u1.*u1+u2.*u2;umax=max(max(ua)); umax=sqrt(umax); 
    cfl=[cfl real(umax*h/dx)]; 
    title(['t= ',num2str(t),',   Rey=',num2str(rey),', n=',num2str(n),', h=', num2str(h),', max|u|=', num2str(real(umax))]);
% ----- calculate and plot vorticity  
    u1h=fft2(u1); u2h=fft2(u2); voh=vof(u1h,u2h); vo=ifft2(voh); % diff in Fourier 
    figure(2); clf; pcolor(x,y,real(vo));colormap gray; shading interp; colorbar;
    cmax=max(max(abs(vo))); 
    title(['t= ',num2str(t),',   Rey=',num2str(rey),', n=',num2str(n),', h=',...
        num2str(h),', max(|curl|)=', num2str(cmax)]);
% --  plot vorticity-hat, useful for 'complicated' (turb) flows
    vohs=(abs(voh)+1e-10).^2/n^4;  
    figure(3); clf;surf(fftshift(KX),fftshift(KY),log10(fftshift(vohs)));
    axis([-n/2 n/2 -n/2 n/2 -30 0]); view(-10,60); 
    tits=['t= ',num2str(t),',   Rey=',num2str(rey),', n=',num2str(n),', log(|w_k|^2)']; 
    title(tits); colormap gray;shading interp;colorbar;
% -- scatter-plot of vorticity power-spectrum to check scaling behaviour 
    figure(4); clf; loglog(kmod,abs(vohs)+1e-20,'b *','MarkerSize',10);
    title(tits); axis([1 n/2 1e-8 1]); 
    hold on;loglog(2:10,(2:10).^(-2),'k');text(10,.05,'K^{-2}');
% -- calc and plot pressure, usually not needed 
   %p=pcalc(u1,u2,n,t);figure(5); clf; surf(X,Y,real(p)); 
   %set(gca,'FontSize',16);axis([-pi pi -pi pi]);view(-10,60);colormap Gray;
   %title(['t= ',num2str(t),',   Rey=',num2str(rey),',  p']);  
    end;  % outer int. loop, ------------------------------------------- 
% return control to user to ask for Rey etc 
    rey=asknu('Rey',rey); tp=asknu('dT (plot-interval)',tp); intsteps=ceil(10*tp/dx); 
    'Make sure that CFL=umax*h/dx<1/2. So far:', cfl
    h=tp/intsteps;  h=asknu('h',h); intsteps=ceil(tp/h); makemult(n,rey,h); 
    np=asknu('Use each np-th point in quiver plots', np);  
    exsteps=asknu(['new number of steps of length ',num2str(tp)],exsteps);
end


        