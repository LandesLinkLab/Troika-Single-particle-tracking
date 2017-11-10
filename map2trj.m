function trjR = map2trj(R, map_new, frame_position, time)
% correct add the particles in the new frame to the trajectory matrix
for i = 1 : numel(map_new(:,1))
    if map_new(i, 2) ~= 0
        R(time, :, i) = frame_position(map_new(i, 2),:);
    elseif map_new(i, 2) == 0
        R(time, :, i) = 0;
    end
end
trjR = R;
end