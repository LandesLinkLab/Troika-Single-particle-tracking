function [map_new, termied] = linking_map(map_pr, map1_2, position, time, termi, search_r)
% this function transfer the map between two consecutive frames to the map
% corresponding to the trajectory matrix
map_new = map_pr;
termied = [];
nueve = [];
makeup = [0,0];
if isempty(map1_2)
    map_new(:,2) = 0;
else
    e = find(map1_2(:,1)~=0,1,'last');
    nueve = map1_2(e+1:end,2);
    nueve(nueve==0)=[];
    if time >= 2 && ~isempty(termi)
        makeup = mapping_frames(position(time-1).p(termi(:,2),:), position(time+1).p(nueve,:), search_r);
        for k = 1:size(makeup,1)
            if makeup(k,1) ~= 0 && makeup(k,2) ~= 0
                makeup(k,:) = [termi(k,1), nueve(makeup(k,2))];
            elseif makeup(k,1) ~= 0
                makeup(k,1) = termi(k,1);%make sure the labels are correct
            end % for k
        end % for time
    end % if time
for i = 1 : numel(map_pr(:,2))
    if map_pr(i, 2) ~= 0
        map_new(i, 2) = map1_2(map_pr(i, 2), 2);
        if map1_2(map_pr(i, 2), 2) == 0
            termied(end+1,:) = [i, map_pr(i, 2)];
        end
    elseif ~isempty(makeup)
        if ~isempty(makeup(makeup(:,1)==i))
            map_new(i, 2) = makeup(makeup(:,1)==i,2);
        end
    else
        map_new(i, 2) = 0;
    end %if map_pr
end % for i
for j = max(map_pr(:,2)) + 1 : numel(map1_2(:,2))
    i = i + 1;
    map_new(i, 2) = map1_2(j, 2);
    map_new(i, 1) = i;
end % for
end % if isempty
end % function