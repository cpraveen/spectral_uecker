function huwf(X,Y,Z,c,lw) % emulate matlab-waterfall-plot in octave, 
% c=color, e.g., 'k', lw=linewidth, e.g.=2
n=length(X(:,1)); 
plot3(X(1,:),Y(1,:),Z(1,:),'linewidth',lw,c); hold on;
for i=2:n;  plot3(X(i,:),Y(i,:),Z(i,:),'linewidth',lw,c); end
hold off;
end 