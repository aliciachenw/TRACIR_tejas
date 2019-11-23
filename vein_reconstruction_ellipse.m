%US_Scan_flat
%US_image
tic()
clear;
clc;
clf; 
close all;
count = 1;

%% define wanwen parameters
params.pix_shift_x = 1;
params.pix_shift_y = 1;
params.thresh_bin = 4; %4       
params.gauss_sigm = 3;

params.min_rad = 30;
params.max_rad = 50;
params.half_window = 1;

start_frame = 40;
%% read files
data_file = 'datasets/data_14Sep_2';

load('calibration.mat');
load('usprobe_pose.mat');
addpath(genpath('rvctools/'));
addpath('.\TRACIR_wanwen');
load_probe_position(data_file);
imageList = dir(strcat(data_file,'/*.jpg'));
image_frame_origin = uint32([158,37]);
% image_frame_origin = uint32([0,0]);


max_circ_per_img = 3;
center3d_list = zeros(size(imageList,1),3*max_circ_per_img ); 

resluts_save_dir = strcat(data_file,'/results');
mkdir(resluts_save_dir)
fig_reconstruction = figure;
daspect([1 1 1]);
 axis([-90 -60 600 620 340 440])
init_flag = 0;
 % fig_centroid_3d = figure;
%% run code
for i=start_frame: size(imageList,1)
    I_original = imread(strcat(data_file,'/',imageList(i).name) );
    
    if i == start_frame
        disp("select imaging area")
        [J,rect1] = imcrop(I_original);
        rect1 = uint32(rect1);
        rect_coord = [rect1(2),rect1(2)+rect1(4), rect1(1),rect1(1)+rect1(3)];
        offset = [rect_coord(3),rect_coord(1)] ;
        params.cut_xmin = rect1(2)+10;
        params.cut_ymin = rect1(1)+10;
        params.cut_height = rect1(4)-10;
        params.cut_width = rect1(3)-10;
%         fig_centroid_3d = figure;
    end
    
    if init_flag == 0
        disp("select initial vessle area")
        [J,rect2] = imcrop(I_original);
        start_point = [rect2(2)+rect2(4)/2,rect2(1)+rect2(3)/2];
        params.start_point = start_point;
        init_flag = 1;
    end
    
    img_out = shift_filter_tejas(I_original,params);
%     I2= I_original(rect_coord(1):rect_coord(2), rect_coord(3):rect_coord(4));
    disp(i)
    
    q=UnitQuaternion(us_pose(i,4),us_pose(i,5:7));
    p=us_pose(i,1:3)'*1000;
    R=q.R;
    %[centers,radii] = imfindcircles(I2,[30 70],'Sensitivity',0.915) %0.915
    %     [BWsdil,centers,radii] = robust_circle_v1(I2);
    % [circle, edge_points] = circle_detection_wanwen_v2(img_out,params,'circle');
    [ellipses, edge_points] = circle_detection_wanwen_v2(img_out,I_original, params, 'ellipse');
    centers = [];
    if ~isempty(ellipses) 
        stpoints = [];
        for el = 1:length(ellipses)
            ellipse = ellipses(el);
            centers = [ellipse.yc, ellipse.xc];
            stpoints = [stpoints;ellipse.xc,ellipse.yc];

            plot_on = 1;
            if plot_on == 1
                fig2 = figure;
                figure(fig2);
                [m,n] = size(I_original);
                if 0<ellipse.xc && ellipse.xc<m && ellipse.yc>0 && ellipse.yc<n
                    I_original(round(ellipse.xc), round(ellipse.yc), 1) = 255;
                    I_original(round(ellipse.xc), round(ellipse.yc), 2) = 255;
                    I_original(round(ellipse.xc), round(ellipse.yc), 3) = 0;
                end
                for theta = 0:0.01:2*pi
                    x = ellipse.xc+ellipse.a*cos(theta)*cos(ellipse.alpha)-ellipse.b*sin(theta)*sin(ellipse.alpha);
                    y = ellipse.yc+ellipse.a*cos(theta)*sin(ellipse.alpha)+ellipse.b*sin(theta)*cos(ellipse.alpha);
                    a = round(x);
                    b = round(y);
                    if a > 0 && a<m && b>0 && b<n
                        I_original(a, b, 1) = 255;
                        I_original(a, b, 2) = 0;
                        I_original(a, b, 3) = 0;
                    end
                end
                imshow(I_original);
                hold on;

                pause()
                close (fig2) ;
             end
             centers = centers-double(image_frame_origin);
             contri(count) = i;
             count = count+1;
        end
        params.start_point = stpoints;
    end
    
    if(length(ellipses)==0)
        continue;
    end
    for j=1:length(ellipses)
        plotEllipse3D(fig_reconstruction,ellipses(j),R,Rcal,sx,sy,pcal,p,image_frame_origin,i);
%         pause()
%         plotCenters3D(fig_centroid_3d,center_3d,i);
        hold on;
        if j > max_circ_per_img
            break
        end
            
    end

    %% saving results image by image
%     savename = strcat('/reconstruct_result_',string(i),'.png');
%     full_path = strcat(resluts_save_dir,savename);
%     saveas(fig_reconstruction,full_path);
end
figure (1)
daspect([1 1 1]);
axis vis3d;
title('Ultrasound vein reconstruction using circle hough transform');
xlabel('x(mm)');
ylabel('y(mm)');
zlabel('z(mm)');

toc()