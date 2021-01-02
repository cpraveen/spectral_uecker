function huwf(X,Y,Z)
% emulate matlab-waterfall-plot in octave
% X,Y,Z=meshgrids
n=length(X(:,1)); 
plot3(X(1,:),Y(1,:),Z(1,:)); hold on;
for i=2:n;  plot3(X(i,:),Y(i,:),Z(i,:)); end
hold off;
end 