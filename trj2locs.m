function [ clocs ] = trj2locs( trjR )
%trj2locs Hacks appart trjR to return cell array of localizations

[~,D,~] = size(trjR);
% trjBin = squeeze(trjR(:,1,:) ~= 0);
% numlocs = sum(trjBin,2); % Number of locs in each frame
numlocs = sum(squeeze(trjR(:,1,:) ~= 0),2); % Number of locs in each frame

% Rearrange trjR so we can get a full list of localizations, looping
% through frames, trajectories, and dimensions
trjR = permute(trjR,[3,1,2]);
trjR = trjR(trjR~=0);
newN = numel(trjR)/D;
trjR = reshape(newN,D);

clocs = mat2cell(trjR,numlocs,D);
