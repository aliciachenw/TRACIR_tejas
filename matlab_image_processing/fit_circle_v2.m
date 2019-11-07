function circle = fit_circle_v2(XY)
% a more robust fitting algorithm based on the least-square approximation
% reference: https://www.mathworks.com/matlabcentral/fileexchange/15060-fitcircle-m
B = [XY(:,1).^2+XY(:,2).^2, XY(:,1),XY(:,2),ones(size(XY,1),1)];
[U,S,V] = svd(B);
u = V(:,4);
a = u(1);
b = u(2:3);
c = u(4);
centroid = -b./(2*a);
circle.xc = centroid(1); 
circle.yc = centroid(2); 
circle.rad = sqrt((norm(b)/(2*a))^2-c/a);
end 