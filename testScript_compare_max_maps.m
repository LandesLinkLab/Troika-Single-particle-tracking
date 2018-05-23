wide = 2;
T = size(Data1,3);

max_maps = zeros(size(Data1));
max_maps_R = max_maps;

wb = waitbar(0,'Calculating max maps...');
tic
for t= 1:T
    max_maps(:,:,t) = calculate_max_map(Data1(:,:,t),wide);
%     max_maps_R(:,:,t) = calculate_max_map_R(Data1(:,:,t),wide);
    waitbar(t/T,wb)
end
toc
close(wb)

wb = waitbar(0,'Calculating max maps...');
tic
for t= 1:T
%     max_maps(:,:,t) = calculate_max_map(Data1(:,:,t),wide);
    max_maps_R(:,:,t) = calculate_max_map_R(Data1(:,:,t),wide);
    waitbar(t/T,wb)
end
toc
close(wb)

equalcheck = isequal(max_maps,max_maps_R)