wide = 2;
[Y,X,T] = size(Data1);
Data2 = SNR_booster(Data1);

max_maps = zeros(size(Data1));
max_maps_R = max_maps;

wb = waitbar(0,'Calculating max maps...');
locs = cell(T,1);
locsR = cell(T,1);
tic
for t= 1:T
    tmp = calculate_max_map(Data2(:,:,t),wide);
    max_maps(:,:,t) = tmp;
    [R,C] = ind2sub([Y,X],find(tmp));
    locs{t} = [C,R];
%     max_maps_R(:,:,t) = calculate_max_map_R(Data1(:,:,t),wide);
    waitbar(t/T,wb)
end
toc
close(wb)

wb = waitbar(0,'Calculating max maps...');
tic
for t= 1:T
%     max_maps(:,:,t) = calculate_max_map(Data1(:,:,t),wide);
    tmp = calculate_max_map_R(Data2(:,:,t),wide,.9);
    max_maps_R(:,:,t) = tmp;
    [R,C] = ind2sub([Y,X],find(tmp));
    locsR{t} = [C,R];
    waitbar(t/T,wb)
end
toc
close(wb)

equalcheck = isequal(max_maps,max_maps_R)
figure(1)
plot(squeeze(sum(sum((max_maps-max_maps_R)))));

LocVisGUI(Data1,locs,locsR)