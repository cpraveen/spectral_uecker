% u'=u, explicit vs implicit Euler
n=50; h=0.1; u=zeros(n,1);v=zeros(n,1);w=zeros(n,1);
t=0:h:(n-1)*h;u(1)=1;v(1)=1;w(1)=1;
for k=1:n-1;
    u(k+1)=u(k)+h*u(k); v(k+1)=v(k)/(1-h);
end
clf;plot(t,exp(t),t,u,'-ko',t,v,'-k*');
break 
f=@(t,x) x; [t,x] = ode45(f,[0 5],1); 
clf;plot(t,exp(t),t,x,'-ks');