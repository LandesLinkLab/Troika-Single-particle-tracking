im = SNR_booster(Data1(:,:,1));
wide = 2;

tic
max_map = calculate_max_map(im,wide);
toc

tic
max_map_R = calculate_max_map_R(im,wide);
toc

figure(15274)
subplot(1,2,1)
imagesc(max_map);
axis image off

subplot(1,2,2)
imagesc(max_map_R);
axis image off