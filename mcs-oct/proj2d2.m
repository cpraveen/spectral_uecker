function [u1,u2,p]=proj2d2(f1,f2);
global p1 p2 p3 p4 p5; 
u1=p1.*f1+p2.*f2; u2=p2.*f1+p3.*f2;
p=p4.*f1+p5.*f2;   % done for zero mean flow (comment out next line)
u1(1,1)=f1(1,1); u2(1,1)=f2(1,1); %for periodic pressure set mean flow
end