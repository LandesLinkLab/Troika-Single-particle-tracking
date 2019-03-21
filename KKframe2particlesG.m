function positions = KKframe2particlesG( frames, varargin )
%KKframe2particlesG Applies the Gaussian fitting to SNR boosted frames
%and returns localizations as a structure
%   INPUT:  frames - input frames.
%           varargin - { local_thd , Gauss_width , wide2 , num_std }
%               local_thd - Default = true. Determines whether to use local
%                   thresholding. For a good image, local thresholding is 
%                   not necessary.
%               Gauss_width - Default = 3. The estimated Gaussian standard 
%                   deviation. This parameter will decide the size of the 
%                   fitting region.
%               wide2 - Default = 2. The wide threshold. Note - wide and
%                   wide2 are different. Together they can decide how many
%                   neighbors to include. More detail is discussed in Fig 6
%                   in the paper.
%               num_std - Default = 3. How many standard deviations to add
%                   up as a threshold.
%% particle identification
T = size(frames,3);
invT = 1/T;
% Initialize variable positions using the same fcn
positions.p = particle_identify([]);
positions = repmat(positions,1,T);
wb = waitbar(0,'Localizing particles');
for t = 1 : size(frames,3)
    im = SNR_booster(frames(:,:,t));
    positions(t).p = particle_identify_Gaussfit(im,varargin{:});
    waitbar(t*invT,wb);
end
close(wb);