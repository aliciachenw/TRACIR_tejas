function [ellipse, edge_points] = circle_detection_wanwen_v2(binary_image,params, type)
%%% ----------------------------
% input:
%
% binary_img: input binary image
% params: input parameters (include us image ranges, start points, minimal
% radius, maximal radius, window's half length
% type: 'circle' or 'ellipse' ('ellipse' fitting doesn't work for now)
%
% output:
% ellipse: detect circle / ellipse (circle struct: has xc, yc, rad)
% edge_points: the points that are used for fitting
%%%
img_out = binary_image;
[img_width, img_height] = size(binary_image);

start_point = params.start_point;
min_rad = params.min_rad;
max_rad = params.max_rad;
cut_xmin = params.cut_xmin;
cut_ymin = params.cut_ymin;
cut_width = params.cut_width;
cut_height = params.cut_height;
half_window = params.half_window;

% mask the non-ultrasound part in the image
for i = 1:img_width
    for j = 1:img_height
        if (i>=cut_xmin && i<=cut_xmin+cut_height) && (j>=cut_ymin && j<=cut_ymin+cut_width)
        else
            img_out(i,j)=0;
        end
    end
end
%imtool(img_out);
%pause;

% ray
step_theta = 0.01;
theta_range = 0:step_theta:2*pi;
rad_range = linspace(min_rad, max_rad, 10);
edge_points = []; % selected points

for i = 1:length(theta_range)
    candidate_bright = 0;
    candidate_points = [];
    for j = 1:length(rad_range)
        a = round(start_point(1)+rad_range(j)*cos(theta_range(i)));
        b = round(start_point(2)+rad_range(j)*sin(theta_range(i)));
        if img_out(a,b) ==255 && mean(mean(img_out(a-half_window:a+half_window,b-half_window:b+half_window))) >255*0.3 ...
            && candidate_bright < mean(mean(img_out(a-half_window:a+half_window,b-half_window:b+half_window)))
            candidate_points = [a,b];
            candidate_bright = mean(mean(img_out(a-half_window:a+half_window,b-half_window:b+half_window)));
        end
    end
    if ~isempty(candidate_points)
        edge_points = [edge_points;candidate_points];
    end
end

if strcmp(type,'circle')
    if size(edge_points,1)>3
        ellipse = fit_circle_v2(edge_points);
    else
        ellipse = [];
    end
    %{
    %% output the results
    if ~isempty(ellipse)
        img_out = binary_image;
        img_out(round(ellipse.xc),round(ellipse.yc)) = 255; % centroid
        for pt = 1:size(edge_points,1)
            img_out(round(edge_points(pt,1)),round(edge_points(pt,2))) = 100; % detected points in img for fitting
        end
        for theta = 0:0.01:2*pi % fitted circle
            a = round(ellipse.xc+ellipse.rad*cos(theta));
            b = round(ellipse.yc+ellipse.rad*sin(theta));
            if a > 0 && a<size(img_out,1) && b>0 && b<size(img_out,2)
                img_out(a,b) = 255;
            end
        end
        figure;
        montage({binary_image, img_out});
        pause;
    end
    %}
else 
    if size(edge_points,1)>6
        ellipse = fit_ellipse_v3(edge_points); % ellipse fitting (doesn't work now)
    else
        ellipse = [];
    end

end

