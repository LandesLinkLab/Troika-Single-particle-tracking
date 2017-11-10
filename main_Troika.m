clear position trjR
search_r = 6; % input the search radius
%% load the file
% this part require the file to be a .mat file, which contain a 3
% dimensional structure. The three dimensions are y, x, and t.
[filename, filepath] = uigetfile;
frames = importdata([filepath,filename]);
if isstruct(frames)
        frames = frames.Data1;
end
%% increase SNR and then identify particles for each frame
for t = 1 : size(frames,3)
    im = SNR_booster(frames(:,:,t));
    position(t).p = particle_identify(im);
end
%% mapping
disp('mapping')
% initial the trajectory matrix: first diminsion is time, second is x, y,
% and Gaussian width, third dimension is the particle number.
for p = 1 : numel(position(1).p(:,1))
    trjR(1, :, p) = position(1).p(p, :);
end
map_new = [1: p; 1: p]';% match the map between two frames to particle numbers in trjR
termied = [];% record the terminated particles
for time = 1 : numel(position)-1
    time;
    map1_2 = mapping_frames(position(time).p, position(time+1).p, search_r);
    [map_new, termied] = linking_map(map_new, map1_2, position, time, termied, search_r);
    trjR = map2trj(trjR, map_new, position(time + 1).p, time + 1);
end