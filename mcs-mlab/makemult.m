function makemult(n,rey,h)  
% generate multiplier matrices for ns2d.m; 
% av,bv for diff in Fourier (dx resp dy)
% mm for time stepping; p1,p2,p3 needed in projection step, p4,p5 in pcalc
% matrices are generated as if FFT was centered, i.e., 
% k=(* -n/2+1 ... -1 | 0 1 ... n/2-1) 
% ('0' at position n/2+1, * is unused for even n) 
% and then swapped by fftshift to confirm to matlab fft 
% such that fftshift(k)=(0 1 ... n/2-1 | * -n/2+1 ... -1) 
'calculating multipliers..'
global av bv mm p1 p2 p3 p4 p5;
for j1=1:n % fill av,bv,mm 
    k1=j1-n/2-1; % attention: j1=row-index corresponds to k2=x-Four.-coeff. 
    av(1,j1)=-1i*k1;bv(j1,1)=-1i*k1;
    for j2=1:n
        k2=j2-n/2-1;ks=k1*k1+k2*k2; mm(j1,j2)=1/(1+h*ks/rey);
    end
end
for j1=1:n % fill p1,p2,p3 
    for j2=1:n
        k1=j2-n/2-1;k2=j1-n/2-1; 
        % attention: j1=row-index corresponds to k2=y-Four.-coeff. 
        k1s=k1*k1; k2s=k2*k2; ks=k1s+k2s;
        if (ks==0) p1(j1,j2)=0; p2(j1,j2)=0; p3(j1,j2)=0; 
            p4(j1,j2)=0; p5(j1,j2)=0; 
        else mk=1/ks;
            p1(j1,j2)=mk*k2s; p2(j1,j2)=-mk*k1*k2; p3(j1,j2)=mk*k1s;
            p4(j1,j2)=-mk*1i*k1; p5(j1,j2)=-mk*1i*k2; % for p-calc 
        end
    end
end
av=fftshift(av); bv=fftshift(bv); mm=fftshift(mm);p1=fftshift(p1); 
p2=fftshift(p2); p3=fftshift(p3); p4=fftshift(p4); p5=fftshift(p5); 
end