function im_out = SNR_booster(frame)
%% to increase SNR in each frame by convolution
% this emthod is highly inspired by Arnauld Serge's code for particle
% tracking, for more information about their method, look up: "Dynamic
% multiple-target tracing to probe spatiotemporal cartography of cell
% membranes" and their online code.
% There is a easier function in matlab: conv2, which can also perform 2-D
% convolution of two images. It is even faster. However, that function
% suffering from the boundary problem, that's why we don't use it here. Here
% we use a Fourier transform method. Petersen 1993, Eqn. 12 - Fourier method
%
% UPDATED 11/28/17 by Rashad Baiyasi to function for T frames

% input: frame, or a single image as a 2-D matrix.
% output: im_out, a 2-D matrix after the convolution
[v, h ,T] = size(frame);
ex_w = 3;% we always use a 3X3 matrix to convolute with the raw image.
% to use Fourier transform for convolution, the mask need to has the same
% size as the image.
mask = zeros(v, h);
im_out = zeros(v,h,T);
mask(ceil(v/2):ceil(v/2)+ex_w-1,ceil(h/2):ceil(h/2)+ex_w-1) = ones(ex_w, ex_w);
% Convolution
ft_mask = fft2(mask);
for t = 1:T
    ft_im = fft2(frame(:,:,t));
    im_out(:,:,t) = real(fftshift(ifft2(ft_mask .* ft_im))) / ex_w^2;
end