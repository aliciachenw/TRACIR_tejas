function [ img_out] = shift_filter_tejas(img_base,params)
%SHIFT_FILTER_TEJAS Summary of this function goes here
 

% for data_14Sep_2 

pix_shift_x=params.pix_shift_x ;
pix_shift_y=params.pix_shift_y ;
thresh_bin=params.thresh_bin;
gauss_sigm=params.gauss_sigm;


%% filtering 
img_base_gauss =imgaussfilt(img_base,gauss_sigm);
% img_base_gauss =img_base;

%% shift image and take average
img_left=uint8(zeros(size(img_base_gauss)));
img_right=uint8(zeros(size(img_base_gauss)));
img_top=uint8(zeros(size(img_base_gauss)));
img_bottom=uint8(zeros(size(img_base_gauss)));


img_bottom(pix_shift_y+1 :end,:,:) = img_base_gauss(1:end-pix_shift_y, :, :);
img_top(1:end-pix_shift_y ,:,:) = img_base_gauss(pix_shift_y+1 :end, :, :);
img_right(:,pix_shift_x +1:end,:) = img_base_gauss(:,1:end-pix_shift_x, :);
img_left(:,1:end-pix_shift_x,:) = img_base_gauss(:,pix_shift_x+1:end, :);
 
I1 = img_base_gauss-img_left;  
I2 = img_base_gauss-img_right;
I3 = img_base_gauss-img_top;
I4 = img_base_gauss-img_bottom;

avg_img =(I1+I2+I3+I4)/4; 
a=rgb2gray(avg_img );
 

img_out=uint8(zeros(size(a)));
[idx,idy] = find(a>thresh_bin); 
for i = 1:size(idx)
    img_out(idx(i),idy(i)) = 255;
end 

end

