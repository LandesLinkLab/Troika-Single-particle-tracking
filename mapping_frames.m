function map1_2 = mapping_frames(frame_position1, frame_position2, search_r)
%{
stopped using from 081314, Bo
if isempty(frame_position2) || isempty(frame_position1)
    % in case there is problem with the image
    map1_2 = [];
    return
end % if 
%}
%% modified by Bo, 081314
if isempty(frame_position2) && isempty(frame_position1)
    map1_2 = [];
    return
elseif isempty(frame_position2)
    nop1 = numel(frame_position1(:, 1));
    map1_2 = zeros(nop1, 2);
    map1_2(:,1) = 1:nop1;
    return
elseif isempty(frame_position1)
    nop2 = numel(frame_position2(:, 1));
    map1_2 = zeros(nop2, 2);
    map1_2(:,2) = 1:nop2;
    return
end
%%
nop1 = numel(frame_position1(:, 1));
nop2 = numel(frame_position2(:, 1));
map_2 = [];% record all the connections, solid and dashed, see Fig 8 of the paper
map1_2 = zeros(nop1, 2);% record the solid connections
left_1 = zeros(nop1, 1);% record if the particle has been mapped
left_2 = zeros(nop2, 1);% record how many times the particle has been mapped
%% first search
% As described in Fig 8 (b) and (c), for each partilce in frame1, we search
% for its neighbors within the search distance in frame2.
for i = 1 : nop1
    candi = find(abs(frame_position2(:, 1) - frame_position1(i, 1)) < search_r);
    candi = candi(abs(frame_position2(candi, 2) - frame_position1(i, 2)) < search_r);
    if numel(candi) == 1
        % easy case, only one neighbor found
        map_2(end+1,:) = [i, candi];
        map1_2(i, :) = [i, candi];
        left_1(i) = 1;
        left_2(candi) = left_2(candi) + 1;
    elseif numel(candi) > 1
        % the case that more than one neighors were found, 
        % find the nearest one
        dis = ((frame_position2(candi, 1) - frame_position1(i, 1)).^2 + ...
            (frame_position2(candi, 2) - frame_position1(i, 2)).^2).^0.5;
        dis = [dis, candi];
        dis = sortrows(dis, 1);
        map_2 = [map_2; i*ones(numel(candi),1), candi];
        map1_2(i, :) = [i, dis(1, 2)];% only record the shortest one
        left_1(i) = 1;
        left_2(dis(1, 2)) = left_2(dis(1, 2)) + 1;
    end % if numel
end % for i
indx = find(left_2 > 1);% search for local graph with conflicting
group1 = struct([]);% group1.f1 is the particle in frame1, and the group1.f2
% are particles linked to the group1.f1 particle
group2 = struct([]);% group2.f1 is the particle in frame1, and the group2.f2
% are particles linked to the group2.f1 particle
for i = 1 : numel(indx)
    group2(i).f1 = indx(i);
    group2(i).f2 = map_2((map_2(:, 2) == indx(i)), 1);
    for k = 1 : numel(group2(i).f2)
        group1(end+1).f1 = group2(i).f2(k);
        group1(end).f2 = map_2((map_2(:, 1) == group2(i).f2(k)), 2);
    end % for k
end % for i
%% Solve the local global optimum
while ~(isempty(group1) && isempty(group2))
    [f1, f2, group1, group2, edge] = graph_group(group1, group2);
    edge = sortrows(edge);
    ind = [];
    for i = 2 : numel(edge(:, 1))
        if sum(edge(i, :) == edge(i - 1, :)) == 2
            ind(end+1) = i;
        end
    end % for i
    edge(ind, :) = [];% remove the repeat links
    dis = ((frame_position2(edge(:,2), 1) - frame_position1(edge(:,1), 1)).^2 + ...
        (frame_position2(edge(:,2), 2) - frame_position1(edge(:,1), 2)).^2).^0.5;
    [dis1, dis2, map1, map2, map_plus] = short_dis(f1, f2, edge, dis);
    if dis2 > dis1 
        map_plus = [map_plus; map1];
        map1 = [];
    end
    if ~isempty(map_plus) % add particles with no links
        for i = 1 : numel(f1)
            if isempty(find(map_plus(:, 1) == f1(i),1))
                map_plus(end+1, :) = [f1(i), 0];
            end
        end % for i
    end % if 
    if ~isempty(map_plus)
        for i = 1 : numel(f2)
            if isempty(find(map_plus(:, 2) == f2(i),1))
                map_plus(end+1, :) = [0, f2(i)];
            end
        end % for i
    end % if
    % record the links in the output variable: map1_2
    for i = 1 : size(map_plus, 1)
        if map_plus(i, 1) ~= 0
            map1_2(map_plus(i, 1), :) = map_plus(i, :);
            left_1(map_plus(i, 1), :) = 1;
        else
            map1_2 = [map1_2; map_plus(i, :)];
        end
    end % for i
end % while
%% Final step
% increase the search distance by 1.5 to search for the neighbors for the
% lefted particles. It works in the very similar way as the first step.
single = find(left_1 == 0);
for i = 1 : numel(single)
    candi = find(abs(frame_position2(:, 1) - frame_position1(single(i), 1)) < 1.4*search_r);
    candi = candi(abs(frame_position2(candi, 2) - frame_position1(single(i), 2)) < 1.4*search_r);
    ind = [];
    for j = 1 : numel(candi)
        if ~isempty(find(map1_2(:, 2) == candi(j),1))
            ind(end + 1) = j;
        end
    end % for j
    candi(ind) = [];
    if numel(candi) == 1
        map1_2(single(i), :) = [single(i), candi];
    end
end % for i

for i = 1 : nop2
    if isempty(find(map1_2(:,2) == i,1))
        map1_2(end+1, :) = [0, i];
    end
end % for i
%}
end % function

%% 
function [local1, local2, varargout] = graph_group(group1, group2, varargin)
% graph_group use the particles and links recorded in group1 and group2 to
% search for the isolated local graph with mapping conflicting. This
% function use recursion, so the structure of group1 and group2 are
% specially designed.
if (isempty(group1) && isempty(group2)) && isempty(varargin)
    local1 = [];
    local2 = [];
    return
end
if isempty(varargin)
    if isempty(group1)
        local1 = group2(1).f2';
        local2 = group2(1).f1;
        n = numel(group2(1).f2);
        edge = [group2(1).f2, group2(1).f1*ones(n, 1)];
        group2(1) = [];
        varargout{1} = group1;
        varargout{2} = group2;
        varargout{3} = edge;
        return
    end
    local1 = group1(1).f1;
    local2 = group1(1).f2';
    n = numel(group1(1).f2);
    edge = [group1(1).f1*ones(n, 1), group1(1).f2];
    nuevo = group1(1).f2;
    group1(1) = [];
    % recursion
    [local2, local1, group2, group1, edge] = graph_group(group2, group1, nuevo, local2, local1, flipdim(edge,2));
    varargout{1} = group1;
    varargout{2} = group2;
    varargout{3} = edge;
else
    local1 = varargin{2};
    local2 = varargin{3};
    edge = varargin{4};
    nuevo1 = [];
    if isempty(varargin{1})
        varargout{1} = group1;
        varargout{2} = group2;
        varargout{3} = flipdim(edge,2);
        return
    else
        nuevo = varargin{1};
        to_del = [];
        for i = 1 :numel(group1)
            if ~isempty(find(nuevo == group1(i).f1,1))
                for j = 1 : numel(group1(i).f2)
                    if isempty(find(local2 == group1(i).f2(j),1))
                        local2(end + 1) = group1(i).f2(j);
                        nuevo1(end + 1) = group1(i).f2(j);
                        edge = [edge; group1(i).f1, group1(i).f2(j)];
                    else
                        edge = [edge; group1(i).f1, group1(i).f2(j)];
                    end
                end
                to_del(end + 1) = i;
            end
        end
        group1(to_del) = [];
        nuevo = nuevo1;
        % recursion
        [local2, local1, group2, group1, edge] = graph_group(group2, group1, nuevo, local2, local1, flipdim(edge,2));
        varargout{1} = group1;
        varargout{2} = group2;
        varargout{3} = flipdim(edge,2);
    end
end
end

%%
function [dis1 dis2, map1, map2, map_share] = short_dis(f1, f2, edge, dis)
% search for the shortest distance. The calculation could be huge if too
% many particles are involved.
% input: f1, f2 are particles in frame1 and 2. edge has all the links, dis
% record all the distance of the links.
% output: dis1 and 2 are the total distance corresponding to map1 and 2.
% And they should be the first and second shortest total distance.
% map_share record the links for both dis1 and dis2, and map1, map2 only
% record the different links.
n1 = numel(f1);
n2 = numel(f2);
dis2 = (n1 + n2) * 10;% search_r = 20
dis1 = dis2;
map1 = [];
map2 = [];
[~, Ie, ~] = unique(edge(:, 1));
ini = [0; Ie];
nums = diff(ini);
i = 1;
while i < prod(nums + 1)
    dis0 = 0;
    map = [];
    ind = i;
    new_i = 0;
    for j = 1 : numel(nums)
        if j < numel(nums)
            k = ceil(ind / prod(nums(j+1:end)+1));
            ind = ind - (k-1) * prod(nums(j+1:end)+1);
        else
            k = ind;
        end
        while k <= nums(j)
            if isempty(map) || (sum(map(:,1) == edge(ini(j)+k, 1)) == 0 && sum(map(:,2) == edge(ini(j)+k, 2)) == 0)
                dis0 = dis0 + dis(ini(j)+k);
                map(end + 1, :) = edge(ini(j)+k, :);
                break
            else
                k = k + 1;
                ind = 1;
            end
        end
        if j < numel(nums)
            new_i = new_i + (k-1) * prod(nums(j+1:end)+1);
        else
            new_i = new_i + k;
        end
    end
    if i == new_i
        i = i + 1;
    else
        i = new_i;
    end
    map = sortrows(map);
    if numel(map) < n1 + n2
        dis0 = dis0 + (n1 + n2 - numel(map)) * 10; %
    end
    if ~isequal(map,map1) && ~isequal(map,map2)
        if dis0 < dis1
            dis2 = dis1;
            dis1 = dis0;
            map2 = map1;
            map1 = map;
        elseif dis0 < dis2 && dis0 > dis1
            dis2 = dis0;
            map2 = map;
        end
    end
end
ind1 = [];
ind2 = [];
map_share = [];
for i = 1 : size(map1,1)
    for j = 1 : size(map2,1)
        if sum(map1(i,:)==map2(j,:)) == 2
            map_share = [map_share; map1(i, :)];
            ind1(end+1) = i;
            ind2(end+1) = j;
        end
    end
end
map1(ind1,:) = [];
map2(ind2,:) = [];
end