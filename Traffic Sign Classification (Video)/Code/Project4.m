%% *************************************************************************
% University of Maryland
% ENPM 673 Robotics Perception
% Project 4: Traffic Sign Recognition
% David A Smart
% Due on: May 16, 2018
% *************************************************************************

% clear up memory
% close all
% clear
clc

% debugging?
debugging = 0;

% load in svm classifiers
load('svmB.mat');
load('svmR.mat');

% set up folders and files
scriptsfolder = pwd; % current folder
inputfolder = '../Input/Video/';
ouputfolder = '../Output/';
imagefolder = '../Input/';

cd(inputfolder);
fileinfo = dir('*.jpg'); % all jpegs

% loop through images
for f = 1:length(fileinfo) % 2799
    
    % --------------------------------------------------------------
    % Class   - aprox. frame no. (visualy determined)
    % --------------------------------------------------------------
    % "00001" - 2124,(2505),2715
    % "00014" - 1964,2065
    % "00017" - 932,1022,1366
    % "00019" - 2170,2305
    % "00021" - 770
    %
    % "00035" - 982,1070,1432
    % "00038" - 865,1081,(1445),2475,2715
    % "00045" - 86,199,730,(1102),1910,2350,2395,2410,(2620),2715
    % --------------------------------------------------------------

    % print out
    clc; fprintf('\n frame: %i', f);
    
    cd(scriptsfolder); cd(inputfolder); % switch to video input folder
    Im = imread(fileinfo(f).name); % read in image
    cd(scriptsfolder); % switch back to scripts folder
    
    %% Detection
    
    % Thresholding in HSV Color Space
    fprintf('\n HSV Detection...');
    [Blue_BB, Red_BB,~] = HSV_Detect(Im, debugging);
    
    % Maximally Stable Extremal Region (MSER) Detection
    fprintf('\n MSER Detection...');
    [Blue_BB, Red_BB,~] = MSER_Detect(Im, Blue_BB, Red_BB, debugging);
    
    %% Classification
    
    fprintf('\n Classification...');
    Im_C = HOG_SVM_Classification(Im, Blue_BB, Red_BB, classifierB, classifierR, imagefolder);
    
    %% OUTPUT 
    
    % plot image with classifications
    figure(3); imshow(Im_C);
    
    if debugging
        P42(f) = getframe(figure(2));
    end
    
    % save figure for making video
    P43(f) = getframe(figure(3));
    
    
end

%% make video
cd(scriptsfolder); cd(ouputfolder); % change to ouput folder

if debugging
    vW = VideoWriter('Detection_MSER.avi');
    vW.FrameRate = 30; % fps
    open(vW);
    for f = 1:length(P42)
        writeVideo(vW,P42(f).cdata);
    end
    close(vW);
    clear vW
end

vW = VideoWriter('Classification.avi');
vW.FrameRate = 30; % fps
open(vW);
for f = 1:length(P43)
    writeVideo(vW,P43(f).cdata);
end
close(vW);
clear vW

cd(scriptsfolder) % return to scripts folder

%% ************************************************************************
