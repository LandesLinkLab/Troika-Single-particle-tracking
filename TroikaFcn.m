function [ trjR , position , frames ] = TroikaFcn( frames , varargin )
%TroikaFcn Function form of main_Troika
%   INPUT:  Data1 - raw data, array of size [Y,X,T]
%           varargin = {search_r,local_thd,Gauss_width,wide2,num_std}
%               search_r - search radius for linking trajectories
%               local_thd - to decide if local thresholding is needed
%               Gauss_width - estimated Gaussian standard deviation. This 
%                   parameter will decide the size of the fitting region.
%               wide2 - the wide threshold
%               num_std - How many standard deviations to add up as 
%                   threshold
%           
%   OUTPUT: [ trjR , position , frames ]
%% handling input variables
%def_args: search_r,local_thd,Gauss_width,wide2,num_std
def_args = { 1.5 , true , 3 , 2 , 3};

if nargin > 1
    tmpind = find(~cellfun(@isempty,varargin));
    def_args(tmpind) = varargin(tmpind);
end
[search_r, local_thd, Gauss_width, wide2, num_std] = def_args{:};
clearvars varargin def_args argin %clear up memory
    
%% particle identification
position = KKframe2particles(frames,local_thd,Gauss_width,wide2,num_std);

%% mapping
disp('mapping')
% initial the trajectory matrix: first diminsion is time, second is x, y,
% and Gaussian width, third dimension is the particle number.
% trjR = [];
%ADDED BY RASHAD 3/31/2017 - also changed '1's to 'k's in the for loop
%below
    for k = 1:numel(position)
        if numel(position(k).p(:,1))>0
            disp(['First localizations in frame ',num2str(k)]);
            break %if the k has particles in it, break out and keep that k
        end
    end
    if k == numel(position)
        disp('No particle localizations found, returning empty trjR')
        trjR = [];
        return;
    end
for p = 1 : numel(position(k).p(:,1))
    trjR(k, :, p) = position(k).p(p, :);
end
map_new = [1: p; 1: p]';% match the map between two frames to particle numbers in trjR
termied = [];% record the terminated particles
wb = waitbar(0,'Building trajectories');
tfact = 1/(numel(position)-1);
for time = 1 : numel(position)-1
    time;
    map1_2 = mapping_frames(position(time).p, position(time+1).p, search_r);
    [map_new, termied] = linking_map_R(... %Using Rashad fix to linking_map
        map_new, map1_2, position, time, termied, search_r);
    trjR = map2trj(trjR, map_new, position(time + 1).p, time + 1);
    waitbar(time*tfact,wb);
end
close(wb);
end

