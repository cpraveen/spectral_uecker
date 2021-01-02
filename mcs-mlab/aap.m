function ph=aap(uh,vh) %anti-aliased product, 
% in: uh,vh from fft with n samples
n=length(uh);m=n*3/2;
uhp=[uh(1:n/2) zeros(1,(m-n)) uh(n/2+1:n)]; % pad uhat with zeros 
vhp=[vh(1:n/2) zeros(1,(m-n)) vh(n/2+1:n)]; % pad vhat with zeros 
up=ifft(uhp); vp=ifft(vhp); w=up.*vp; wh=fft(w);
ph=1.5*[wh(1:n/2) wh(m-n/2+1:m)]; % extract F-coefficients 
end