% Remove empty particles first
trjBin = squeeze(trjR(:,1,:) ~= 0); % Extract the 'on' states
trj2 = trjR(:,:,sum(trjBin,1)>0);

% Then do the rest
trjBin = squeeze(trj2(:,1,:) ~= 0);
[T,~,N] = size(trj2);


disps = cell(1,N);
meandisps = zeros(1,N);
mediandisps = meandisps;
trajlen = meandisps;
maxdisps = meandisps;
for n = 1:N
    xs = trj2(:,1,n); ys = trj2(:,2,n);
    ts = find(xs);
    xvect = xs(ts); yvect = ys(ts);        
    tmpdisps = sqrt((xvect(2:end)-xvect(1)).^2 + (yvect(2:end)-yvect(1)).^2);
    disps{n} = tmpdisps;
    meandisps(n) = mean(tmpdisps);
    mediandisps(n) = median(tmpdisps);
    if numel(xvect) < 2
        maxdisps(n) = NaN;
    else
        maxdisps(n) = max(tmpdisps);
    end
    trajlen(n) = numel(xvect)-1;
end

%% figures
lnwdth = 1;
fntsz = 18;

figure(1)
hh = histogram(maxdisps,'LineWidth',lnwdth);
setFont(fntsz)
title('Maximum displacements for each particle')
xlabel('Displacement (px)')
ylabel('Frequency')
set(gca,'LineWidth',lnwdth)

GMModel = fitgmdist(maxdisps',2);
hold on
viewGMModel(GMModel,1)
hold off