%% saveframes.m
%% PART 0.
% 1. 
% Extract and Save the images from the video sequence provided. 
% Divide into training and testing datasets
% Also make coppies where each of the bouys in the training dataset is
% cropped out.

% clean up workspace & memory
close all
clear
clc

% current location is P3_submission/ColorSeg/Scripts/Part0
ScriptsPart0Folder = pwd;
% want to read video frames from P3_submission/ColorSeg/Input
InputFolder = '../../Input'; % set to input folder
% want to save to P3_submission/ColorSeg/Images/TrainingSet/Frames/
TrainingSetFolder = '../../Images/TrainingSet/Frames';
% want to save to P3_submissionColorSeg/Images/TestSet/Frames/
TestSetFolder = '../../Images/TestSet/Frames';
% want to save cropped images to ColorSeg/Images/TrainingSet/CroppedBuoys/
CroppedTrainingSetFolder = '../../Images/TrainingSet/CroppedBuoys';

cd(InputFolder); % change to input folder
vO = VideoReader('detectbuoy.avi'); % Video Frames
k = 1; % initialize index counter
frame = zeros(vO.Height,vO.Width,3,vO.Duration*vO.FrameRate); % initialize frame info
framecropY = zeros(vO.Height,vO.Width,3,vO.Duration*vO.FrameRate); % initialize cropped frame info
framecropR = zeros(vO.Height,vO.Width,3,vO.Duration*vO.FrameRate); % initialize cropped frame info
framecropG = zeros(vO.Height,vO.Width,3,vO.Duration*vO.FrameRate); % initialize cropped frame info

while vO.CurrentTime < vO.Duration % for the entire movie
    
    if k > 1
        cd(ScriptsPart0Folder); % first return to the original folder
        cd(InputFolder); % change to input folder
    end
    
    % original frame
    frame(:,:,:,k) = im2double(readFrame(vO));
    
    if mod(vO.CurrentTime,vO.Duration/20) == 0 % for the 1/10 of the data randomly selected throughout the video
        
        figure(1);
        set(gcf, 'Position', [1536*0.25 864*0.25 1536*0.5 864*0.5]);
        imshow(frame(:,:,:,k)); % display
        pause(0.1);
        
        cd(ScriptsPart0Folder); % first return to the original folder
        cd(TrainingSetFolder); % change to training dataset folder
        imwrite(frame(:,:,:,k),strcat(num2str(k),'.jpg')) % save as <FrameNo>.jpg
        
        % crop to individual bouys
        cd(ScriptsPart0Folder); % first return to the original folder
        cd(CroppedTrainingSetFolder); % change to cropped training dataset folder
        
        clc; fprintf('\n YELLOW \n'); % yellow bouy (left-most)
        roiY = imellipse; % make roi
        pY = wait(roiY);
        maskY = createMask(roiY); % make mask
        maskY(:,:,2) = maskY; maskY(:,:,3) = maskY(:,:,1);
        framecropYtemp = frame(:,:,:,k);
        framecropYtemp(maskY == 0) = 0;
        framecropY(:,:,:,k) = framecropYtemp;
        close(1); figure(2); 
        imshow(framecropY(:,:,:,k)); % display
        set(gcf, 'Position', [1536*0.1 864*0.15 1536*0.8 864*0.75]);
        pause(0.1);
        imwrite(framecropY(:,:,:,k),strcat('Y_',num2str(k),'.jpg')) % save as <Color_FrameNo>.jpg
        close(2);
        
        clc; fprintf('\n RED \n'); % red bouy (middle)
        figure(1);
        set(gcf, 'Position', [1536*0.25 864*0.25 1536*0.5 864*0.5]);
        imshow(frame(:,:,:,k)); % display
        pause(0.1);
        roiR = imellipse; % make roi
        pR = wait(roiR);
        maskR = createMask(roiR); % make mask
        maskR(:,:,2) = maskR; maskR(:,:,3) = maskR(:,:,1);
        framecropRtemp = frame(:,:,:,k);
        framecropRtemp(maskR == 0) = 0;
        framecropR(:,:,:,k) = framecropRtemp;
        close(1); figure(2); 
        imshow(framecropR(:,:,:,k)); % display
        set(gcf, 'Position', [1536*0.1 864*0.15 1536*0.8 864*0.75]);
        pause(0.1);
        imwrite(framecropR(:,:,:,k),strcat('R_',num2str(k),'.jpg')) % save as <Color_FrameNo>.jpg
        close(2);
        
        clc; fprintf('\n GREEN \n'); % green bouy (right-most)
        figure(1);
        set(gcf, 'Position', [1536*0.25 864*0.25 1536*0.5 864*0.5]);
        imshow(frame(:,:,:,k)); % display
        pause(0.1);
        roiG = imellipse; % make roi
        pG = wait(roiG);
        maskG = createMask(roiG); % make mask
        maskG(:,:,2) = maskG; maskG(:,:,3) = maskG(:,:,1);
        framecropGtemp = frame(:,:,:,k);
        framecropGtemp(maskG == 0) = 0;
        framecropG(:,:,:,k) = framecropGtemp;
        close(1); figure(2);
        imshow(framecropG(:,:,:,k)); % display
        set(gcf, 'Position', [1536*0.1 864*0.15 1536*0.8 864*0.75]);
        pause(0.1);
        imwrite(framecropG(:,:,:,k),strcat('G_',num2str(k),'.jpg')) % save as <Color_FrameNo>.jpg
        close(2);
        
    else % for the second half, save as testing data
        cd(ScriptsPart0Folder); % first return to the original folder
        cd(TestSetFolder); % change to test dataset folder
        
        figure(1);
        set(gcf, 'Position', [1536*0.25 864*0.25 1536*0.5 864*0.5]);
        imshow(frame(:,:,:,k)); % display
        pause(0.1);
        
        imwrite(frame(:,:,:,k),strcat(num2str(k),'.jpg')) % save as <FrameNo>.jpg
        
    end

    k = k+1; % increment index
end

cd(ScriptsPart0Folder); % return

