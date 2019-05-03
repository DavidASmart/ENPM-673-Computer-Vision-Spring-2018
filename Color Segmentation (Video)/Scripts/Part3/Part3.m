%% PART 3.
% Buoy  Detection (40) 

%% ------------------------------------------------------------------------

% clear up memory and workspace
close all
clear
clc

% set up folder paths
current = pwd;

% read in EM parameters from 
ParamsFolder = '../../Output/Part2';

% output segmented images to
plot_path = '../../Output/Part3/Frames';

%% ------------------------------------------------------------------------
% get params

cd(ParamsFolder); 
load('EM.mat')
cd(current); % return


%% ------------------------------------------------------------------------
% segment images

for Frame = 1:200 % for all images (training & testing)
    cd(current); % return
    clc
    fprintf('\n Frame: %i',Frame)
    detectBuoy(Frame, muY, muR, muG, covarY, covarR, covarG, plot_path)
end

%% ------------------------------------------------------------------------
% generate video

% read in frames from
inpath = plot_path;

% output video to
outpath = '../../Output/Part3/Video';

cd(current); % return
frames2video(inpath,outpath);

%% ------------------------------------------------------------------------

cd(current); % return
