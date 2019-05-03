%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% David Smart
% ENPM 673 - Perception
% University of Maryland, College Park
% Project 2: Visual Odometry
% 3/25/2018 (due: 3/28/2018)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem Statement
% You are given camera frames from a driving car.
% scripts to extract intrinsic parameters
% The output should be a plot for the trajectory of the camera/car.
%
% You should implement the functions to estimate the Essential Matrix
% and Rotation/Translation Matrices
%
% DO NOT use Matlab’s Computer Vision Toolbox or any third party code...
% ...Except for comparing your results against that.
%
% Your submission SHOULD be a ZIP folder called [YourDirectoryID] proj2.zip
% 1. You will have a parent directory P2 Submission.
% 2. Under P2 Submission/VisualOdometry you will have three sub-folders: code, input and output.
% 3. You SHOULD also submit a report (Report.pdf) under P2 Submission/VisualOdometry folder.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clear workspace
close all
clear;
clc

fprintf('\n Setting Stuff Up \n');
%% Set Up Video
V = VideoWriter('VisualOdometry','MPEG-4');
V.FrameRate = 20;
open(V);

%% Initialize Camera Poses
% my code
% R_my = eye(3);
% t_my = [0; 0 ; 0];
% pos_my2 = [R_my,t_my;0,0,0,1];
% pos_my = pos_my2(1:3,4);
R_1 = eye(3);
t_1 = [0;0;0];
pos_my = [0,0];
% computer vision toolbox
% R_CVT = eye(3);
% t_CVT = [0; 0; 0];
% pos_CVT2 = [R_CVT,t_CVT;0,0,0,1];
% pos_CVT = pos_CVT2(1:3,4);
R_2 = eye(3);
t_2 = [0;0;0];
pos_CVT = [0,0];

%% Setup Image File Locations
ImageFolder = '../Input/Oxford_dataset/stereo/centre';
ImageNames = dir(fullfile(ImageFolder,'*.png')); % gather image file names

%% Extract the Camera Parameters using "ReadCameraModel"
ModelFolder = '../Input/Oxford_dataset/model/';
[fx, fy, cx, cy, G_camera_image, LUT] = ReadCameraModel(ImageFolder, ModelFolder);
% create intrinsic camera matrix
K = [fx, 0, cx; 
    0, fy, cy; 
    0, 0, 1];
cameraParams = cameraParameters('IntrinsicMatrix',K);

fprintf('\n Starting Visual Odometry Loop \n');
%% Main Loop **************************************************************
for k = 15:length(ImageNames) % for all images
    
    %% Read Image Files, Convert to Correct Format, Find Features
    [RGB, BW, points2] = Read_Convert_Detect(LUT,ImageFolder,ImageNames,k);
    
    % after the first image...
    if k > 15
        % -----------------------------------------------------------------
        %% **Use my code**
        [R_my,t_my] = myVisualOdometry(BW_old, BW, points1, points2, K);
        
        % Update Camera Pose
%         pos_my2 = [R_my,t_my;0,0,0,1]*pos_my2; % pos & orientation
%         pos_my = [pos_my, pos_my2(1:3,4)]; % just pos
        t_1 = t_1 + R_1*t_my;
        R_1 =  R_1*R_my;
        pos_my = [pos_my;[t_1(1), t_1(3)]];
        
        % -----------------------------------------------------------------
        %% **use Matlab's Computer Vision Toolbox for comparison**
        
        [R_CVT,t_CVT] = CVT_Compare(BW_old, BW, points1, points2, cameraParams);
        
        % Update Camera Pose
%         pos_CVT2 = [R_CVT,t_CVT';0,0,0,1]*pos_CVT2; % pos & orientation
%         pos_CVT = [pos_CVT, pos_CVT2(1:3,4)]; % just pos
        t_2 = t_2 + R_2*t_CVT';
        R_2 =  R_2*R_CVT;
        pos_CVT = [pos_CVT;[t_2(1), t_2(3)]];        

        % -----------------------------------------------------------------
        %% PLOT
        pos_plot_compare(RGB,pos_my,pos_CVT);
        % -----------------------------------------------------------------
        %% Compare
        dx = t_1(1) - t_2(1);
        dz = t_1(3) - t_2(3);
        drift = sqrt(dx^2 + dz^2);
        fprintf('\n frame no. = %i \t', k);
        fprintf('drift = %3.2f \n', round(drift,2));
        % -----------------------------------------------------------------
        %% Make Video
        frame = getframe(gcf);
        writeVideo(V,frame);
        % frame_No = k;
        % -----------------------------------------------------------------
    end
    
    %% update next iteration
    points1 = points2;
    BW_old = BW;
    
end

%% done :)
close(V)

fprintf('\n DONE! \n')
pause(); % wait for user to close

save('VO.mat','pos_my','pos_CVT')

%% clear workspace
close all
clear;
clc
%% ************************************************************************
