function [ max_map ] = calculate_max_map_R( im , wide )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[v , h] = size(im);
n = 3;
wide2 = wide;
bg = zeros(v,h);
%% Calculated thd_map (basic method)
im_local = sort(im(:));
n_loc = numel(im_local);
bg(1:v, 1:h) = im_local(round(n_loc/2));
sd(1:v, 1:h) = im_local(round(n_loc/2)) - im_local(round(n_loc/4));
thd_map = bg + n * sd;

%% Calculate local maxima
gsdil = imdilate(im,strel('disk',wide));
gsdiff = im-gsdil;
locmaxind = find(gsdiff==0);

max_map = zeros(v,h);
max_map(locmaxind) = im(locmaxind);
max_map(max_map < thd_map) = 0;
max_map = max_map .* (im - bg);