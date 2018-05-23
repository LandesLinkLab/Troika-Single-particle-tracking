function [ max_map ] = calculate_max_map_R( im , wide , pctoverthd)
%UNTITLED Summary of this function goes here
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
    
%% calculate max_map
% the idea of this part also originates from Arnauld Serge's work.
center = im(1+wide:v-wide, 1+wide:h-wide); %remove image borders
max_map = zeros(v, h); %size of full image
pos_check = max_map;
thd_map_og = thd_map;
if numel(thd_map) > 1
    thd_map = thd_map(1+wide:v-wide, 1+wide:h-wide);
end
% each center is compared with its local threshold
max_map(1+wide:v-wide, 1+wide:h-wide) = center > thd_map;
% figure; imagesc(max_map); axis image off
% the two loops below are intended to select the local maximums that meet 
% two conditions: 
% 1, the selected neighbors are not brighter than the center;
% 2, the selected neighbors are also brighter than the local threshold.
%RB - variable wide is used as both the check for how local a maxima is, as
%well as ensuring all values within the circle are over the threshold.

% Reformat loop over specified values rather than checking each round
[tmpjs,tmpis] = meshgrid(-wide:wide,-wide:wide);
ijradval = tmpis.^2 + tmpjs.^2 <= wide2^2;
jvals = tmpjs(ijradval);
ivals = tmpis(ijradval);
% Variable to count the number of neighbors over thd
% neighborcount = zeros(size(im));
nhood = zeros(2*wide+1);
nhood(ijradval) = 1;
neighborcount = imfilter(double(im > thd_map_og),nhood);
% figure(1); imagesc(neighborcount); axis image off
% figure(462); imagesc(nhood); axis image off
% thdval = thd_map(1)

for ttt = 1:numel(jvals)
    i = ivals(ttt);
    j = jvals(ttt);
    % Assign to matrix pox_check within a boarder defined by wide.
    % Shift the image around and compare with center (the original
    % image) to determine locality of maxima, and compare with
    % thd_map to ensure all points are over the local threshold
    pos_check(1+wide:v-wide, 1+wide:h-wide) = ...
        im(1+wide+i:v-wide+i, 1+wide+j:h-wide+j) <= center;
%     neighborcount(1+wide:v-wide, 1+wide:h-wide) = ...
%         neighborcount(1+wide:v-wide, 1+wide:h-wide) + ...
%         double(im(1+wide+i:v-wide+i, 1+wide+j:h-wide+j) > thd_map);
    max_map = max_map .* pos_check;
%         figure(24326); imagesc(max_map); axis image off
    pos_check = zeros(v, h);
end
% max(neighborcount(:))
% min(neighborcount(:))
max_map = max_map .* (im - bg);
max_map(neighborcount < pctoverthd*numel(jvals)) = 0;
