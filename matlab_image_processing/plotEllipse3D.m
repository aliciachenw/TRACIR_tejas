function plotEllipse3D(fig1,ellipse,R,Rcal,sx,sy,pcal,p,image_frame_origin,image_num)
center = [ellipse.yc,ellipse.xc] - double(image_frame_origin);
center_3d=((R*(Rcal*[diag([sx,sy])*center';0]+pcal))+p)';

theta=0:0.1:2*pi;
points = zeros(3, length(theta));
for t = 1:length(theta)
    x = ellipse.xc+ellipse.a*cos(theta(t))*cos(ellipse.alpha)-ellipse.b*sin(theta(t))*sin(ellipse.alpha);
	y = ellipse.yc+ellipse.a*cos(theta(t))*sin(ellipse.alpha)+ellipse.b*sin(theta(t))*cos(ellipse.alpha);
    pt = [y,x] - double(image_frame_origin);
    points(:,t) = ((R*(Rcal*[diag([sx,sy])*pt';0]+pcal))+p);
end

%% to plot circles
figure(fig1);
% xlim([-90 -60 ]); ylim([600 620 ]); zlim([340 440])
% axis([-90 -60 600 620 340 440])
view(-25,7)
daspect([1 1 1]);
%  axis equal
% points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
plot3(points(1,:),points(2,:),points(3,:),'r-');
%% to plot centers only
hold on;
% scatter3(center(1),center(2),center(3),'r*');
textscatter3(center_3d(1),center_3d(2),center_3d(3),string(image_num));
end