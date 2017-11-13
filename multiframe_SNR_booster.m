function [ frames_out ] = multiframe_SNR_booster( frames )
%multiframe_SNR_booster Performs SNR_booster on each frame
%   See SNR_booster.m

T = size(frames,3);
frames_out = zeros(size(frames));
for t = 1:T
    frames_out(:,:,t) = SNR_booster(frames(:,:,t));
end
