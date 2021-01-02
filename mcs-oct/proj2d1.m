function [u1,u2]=proj2d1(f1,f2);
global p1 p2 p3; 
u1=p1.*f1+p2.*f2; u2=p2.*f1+p3.*f2; % done for zero mean flow
%u1(1,1)=f1(1,1); u2(1,1)=f2(1,1); % use this for periodic pressure
end