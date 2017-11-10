clear position trjR
search_r = 6; % input the search radius
%% load the file
% this part require the file to be a .mat file, which contain a 3
% dimensional structure. The three dimensions are y, x, and t.
path='E:\Widefielddata\111716\mat\';
savepath='E:\Widefielddata\111716\spt_r_6\';
cd(path)
[filename, filepath] = uigetfile('*.mat','multiselect','on');

if ~iscell(filename) % if only one file was selected, put it into a cell array for compatability
    temp=filename;
    clear filename
    filename=cell(1);
    filename{1}=temp;
    clear temp
end

for ifile = 1:length(filename)
disp(['Loading ' filename{ifile}])
frames = importdata([filepath,filename{ifile}]);

if isstruct(frames)
        frames = frames.Data1;
end
%% increase SNR and then identify particles for each frame
eventflag=0;
for t = 1 : size(frames,3)
    disp(['Frame ' num2str(t)])
    im = SNR_booster(frames(:,:,t));
    position(t).p = particle_identify(im);
    if ~isempty(position(t).p)
        if eventflag == 0
            eframe=t;
        end
        eventflag=1;
    end
end
if eventflag == 0
    disp(['No Trajectories Identified in ' filename{ifile}])
    continue
end
%% mapping
disp('mapping')
% initial the trajectory matrix: first diminsion is time, second is x, y,
% and Gaussian width, third dimension is the particle number.
for p = 1 : numel(position(eframe).p(:,1))
    trjR(eframe, :, p) = position(eframe).p(p, :);
end
map_new = [1: p; 1: p]';% match the map between two frames to particle numbers in trjR
termied = [];% record the terminated particles
for time = 1 : numel(position)-1
    time;
    map1_2 = mapping_frames_2(position(time).p, position(time+1).p, search_r);
    [map_new, termied] = linking_map(map_new, map1_2, position, time, termied, search_r);
    trjR = map2trj(trjR, map_new, position(time + 1).p, time + 1);
end
save([savepath filename{ifile}(1:end-4) '_spt.mat'],'trjR')
clear position trjR
end