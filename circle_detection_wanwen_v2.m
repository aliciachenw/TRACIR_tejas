function [ellipses, edge_points] = circle_detection_wanwen_v2(binary_image, img_base,params, type)
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

pre_start_points = params.start_point;
min_rad = params.min_rad;
max_rad = params.max_rad;
cut_xmin = params.cut_xmin;
cut_ymin = params.cut_ymin;
cut_width = params.cut_width;
cut_height = params.cut_height;
half_window = params.half_window;
max_I = max(max(img_base));
min_I = min(min(img_base));
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


% cluster
img_label = bwlabel(img_out);
            
% ray
step_theta = 0.01;
theta_range = 0:step_theta:2*pi;
rad_range = linspace(min_rad, max_rad, 50);
dr = rad_range(2) - rad_range(1);
edge_points = []; % selected points

% random start points

gridx = cut_xmin:(2 * min_rad):cut_xmin+cut_height;
gridy = cut_ymin:(2 * min_rad):cut_ymin+cut_width;
grid_labels = zeros(length(gridx)-1,length(gridy)-1);

start_points = [];
for i = 1:size(pre_start_points,1)
    flag = 0;
    for xx = 1:length(gridx)-1
        for yy = 1:length(gridy)-1
            if pre_start_points(i,1)> gridx(xx) && pre_start_points(i,1)<gridx(xx+1) &&...
                    pre_start_points(i,2)>gridy(yy) && pre_start_points(i,2)<gridy(yy+1)
                if grid_labels(xx,yy) == 0
                    grid_labels(xx,yy) = 1;
                    flag = 1;
                    if isempty(start_points)
                        start_points = [start_points;pre_start_points(i,:)];
                    else
                        cf = 0;
                        for st = 1:size(start_points,1)
                            if norm(start_points(st,:)-pre_start_points(i,:)) < 2 * min_rad
                                cf = 1;
                                break
                            end
                        end
                        if cf == 0
                            start_points = [start_points;pre_start_points(i,:)];
                        end
                    end
                    
                    break
                end
            end
        end
        if flag == 1
            break;
        end
    end
end

for xx = 1:length(gridx)-1
    for yy = 1:length(gridy)-1
        if grid_labels(xx,yy) == 0
            cf = 0;
            pt = double([(gridx(xx)+gridx(xx+1))/2,(gridy(yy)+gridy(yy+1))/2]);
            for st = 1:size(start_points,1)
                
                if norm(start_points(st,:)-pt) < 2 * min_rad
                	cf = 1;
                    break
                end
            end
            if cf == 0
            	start_points = [start_points;pt];
            end
        end
    end
end

ct = 1;

for st = 1:size(start_points,1)
    edge_points = [];  
    for i = 1:length(theta_range)
        candidate_bright = 0;
        candidate_points = [];
        for j = 1:length(rad_range)
            a = round(start_points(st,1)+rad_range(j)*cos(theta_range(i)));
            b = round(start_points(st,2)+rad_range(j)*sin(theta_range(i)));
            if img_out(a,b) ==255 && mean(mean(img_out(a-half_window:a+half_window,b-half_window:b+half_window))) >255*0.3 ...
                && candidate_bright < mean(mean(img_out(a-half_window:a+half_window,b-half_window:b+half_window))) 
            %{
                in_a = round(a - 2 * dr * cos(theta_range(i)));
                in_b = round(b - 2 * dr * sin(theta_range(i)));
                out_a = round(a + 2 * dr * cos(theta_range(i)));
                out_b = round(b + 2 * dr * sin(theta_range(i)));
                if img_base(in_a,in_b) < img_base(out_a,out_b)
                    candidate_points = [a,b];
                    candidate_bright = mean(mean(img_out(a-half_window:a+half_window,b-half_window:b+half_window)));
                end
                %}
                candidate_points = [a,b];
                candidate_bright = mean(mean(img_out(a-half_window:a+half_window,b-half_window:b+half_window)));
            end
        end
        if ~isempty(candidate_points)
            edge_points = [edge_points;candidate_points];
        end
    end
    if strcmp(type,'circle')
        if size(edge_points,1)>6
            ellipse = fit_circle_v2(edge_points);
            if ~isempty(ellipse)
                flag = 1;
                if ct ~= 1
                    for el = 1:length(ellipses)
                        if norm([ellipse.xc - ellipses(el).xc,ellipse.yc - ellipses(el).yc]) < min_rad
                            flag = 0;
                            break;
                        end
                    end
                end
            end
            if flag == 1 && ellipse.rad < max_rad && ellipse.rad > min_rad && ...
                    ellipse.xc > cut_xmin - max_rad && ellipse.xc < cut_xmin + cut_height + max_rad && ...
                    ellipse.yc > cut_ymin - max_rad && ellipse.yc < cut_ymin + cut_width + max_rad
                ellipses(ct) = ellipse;
                ct = ct + 1;
                for i = 1:length(theta_range)
                	theta = theta_range(i);
                	x = ellipse.xc+ellipse.rad*cos(theta);
                	y = ellipse.yc+ellipse.rad*sin(theta);
                	img_out(round(x-dr):round(x+dr),round(y-dr):round(y+dr)) = 0;
                end
            end
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
        if size(edge_points,1)>=6
            ellipse = fit_ellipse_v3(edge_points);
            if ~isempty(ellipse)
                flag = 1;
                if ct ~= 1
                    for el = 1:length(ellipses)
                        if norm([ellipse.xc - ellipses(el).xc,ellipse.yc - ellipses(el).yc]) < min_rad
                            flag = 0;
                            break;
                        end
                    end
                end
                if flag == 1 && ellipse.a < 1.5 * max_rad && ellipse.a > 0.5 * min_rad && ...
                        ellipse.b < 1.5 * max_rad && ellipse.b > 0.5 * min_rad && ...
                        max([ellipse.a/ellipse.b, ellipse.b/ellipse.a]) < 2 && ...
                        ellipse.xc > cut_xmin - max_rad && ellipse.xc < cut_xmin + cut_height + max_rad && ...
                        ellipse.yc > cut_ymin - max_rad && ellipse.yc < cut_ymin + cut_width + max_rad
                    % if mean(mean(img_base(ellipse.xc - ellipse.a:ellipse.xc+ellipse.a,ellipse.yc - ellipse.b:ellipse.yc+ellipse.b))) < (max_I + min_I) / 2
                        
                        ellipses(ct) = ellipse;
                        ct = ct + 1;
                        %{
                        for ed = 1:size(edge_points,1)
                            a = edge_points(ed,1);
                            b = edge_points(ed,2);
                            label = img_label(a,b);
                            img_out(img_label == label) = 0;
                        end
                        %}
                        for i = 1:length(theta_range)
                            theta = theta_range(i);
                            x = ellipse.xc+ellipse.a*cos(theta)*cos(ellipse.alpha)-ellipse.b*sin(theta)*sin(ellipse.alpha);
                            y = ellipse.yc+ellipse.a*cos(theta)*sin(ellipse.alpha)+ellipse.b*sin(theta)*cos(ellipse.alpha);
                            img_out(round(x-dr):round(x+dr),round(y-dr):round(y+dr)) = 0;
                        end
                    % end
                end
            end
        end
    end  
end

if ct == 1
    ellipses = [];
end

end


