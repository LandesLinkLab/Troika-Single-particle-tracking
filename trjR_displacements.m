function [ disps  , t0s ] = trjR_displacements( trjR )
%Calculates the displacement from the starting point of each point in a 
%trajectory.
%   INPUT
%       trjR  - Troika output array with size [T,D,N] for T frames, D
%               dimensions, and N particles.
%   OUTPUT
%       disps - Cell array of the displacements from the origin for each
%               position along a particle's trajectory
%       t0s   - Time coordinate of the first frame for each trajectory. 
%               Variable 'disps' starts listing displacements in frame 
%               t0 + 1.

% Initialize variables
[~,~,N] = size(trjR);
t0s = zeros(1,N);
disps = cell(1,N);
 
% Main loop to extract trajectories
for n = 1:N
    xs = trjR(:,1,n); ys = trjR(:,2,n);
    ts = find(xs); % Get all on positions
    if isempty(ts) % If there are no particles, skip it
        t0s(n) = NaN;
        disps{n} = zeros(0,1);
        continue
    end
    t0 = ts(1); % Start frame
    t0s(n) = t0;
    tf = ts(end); % End frame
    if isequal(tf,t0) % If it's only one frame long, skip it
        disps{n} = zeros(0,1);
        continue
    end
    % Otherwise, calculate all the displacements
    xdisp = nan((tf-t0),1);
    xdisp(ts(2:end)-t0) = (xs(ts(2:end))-xs(t0)).^2;
    ydisp = nan((tf-t0),1);
    ydisp(ts(2:end)-t0) = (ys(ts(2:end))-ys(t0)).^2;
    disps{n} = sqrt(xdisp+ydisp);
end