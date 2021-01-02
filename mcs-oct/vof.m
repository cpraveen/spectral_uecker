function vo=vof(f1,f2); % calculate vort. in Fourier
global av bv; n=length(f1); v1=zeros(n,n);v2=zeros(n,n);
for j=1:n % diff in Fourier 
    v1(j,:)=av.*f2(j,:); % loop over x=column index
    v2(:,j)=bv.*f1(:,j);
end
vo=v1-v2;
end