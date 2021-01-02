function [g1,g2]=g(X,Y,t) %a selection of right hand sides. 
Xt=X-pi; Yt=Y-pi; % makes notation easier 
% --------- Uncomment the one you like
%dec=exp(-2*X.*X-2*Y.*Y); g1=-dec.*X; g2=-dec.*Y; % a sink
%dec=exp(-2*X.*X-2*Y.*Y); g1=dec.*X; g2=-dec.*Y; % a saddle
%dec=exp(-2*X.*X-Y.*Y); g1=dec.*X; g2=dec.*Y; % a source
%dec=exp(-2*X.*X-2*Y.*Y); g1=dec.*Y*sin(t); g2=-dec.*X*sin(t); % a vortex
%dec=exp(-X.*X-Y.*Y); g1=dec.*(sin(X+Y)+2*sin(0.2*t)); g2=-dec.*Y; % shear like
dec=exp(-4*Xt.^2-4*Yt.^2); g1=dec*cos(0.0*t).*(2+tanh(Yt)); g2=0.*Yt;% a kick in the middle
% --------- for "turbulence" simulations 
%dec1=exp(-2*(Xt-1).^4-2*(Yt-1).^4);dec2=exp(-2*(Xt+1).^2-(Yt+0.5).^4); 
%g1=dec1.*(Yt-1)*cos(t); 
%g2=-dec1.*(Xt-1)*cos(t)+dec2.*(Xt+1)*cos(t/5); 
%[g1,g2]=lowp2d(g1,g2,6); % low-pass if desired
end


        