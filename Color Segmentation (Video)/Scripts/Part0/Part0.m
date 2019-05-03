%% Part0.m
%% PART 0.
% 3. 
% Model the color distributions of the buoys using 1-D Gaussians. 
% Segment the buoys based on their color (red green and yellow).

% clean up workspace & memory
close all
clear
clc

% current folder is ...
ScriptsPart0Folder = pwd;
% need to read images from both...
TrainingSetFolder = '../../Images/TrainingSet/Frames';
TestSetFolder = '../../Images/TestSet/Frames';
% need to output to ...
OutputFolder = '../../Output/Part0/';

for Frame = 1:200 % for all images
    cd(ScriptsPart0Folder);
    segment1D(Frame); % do 1D segmentation
end


%% make video -------------------------------------------------------------
cd(ScriptsPart0Folder);cd(OutputFolder);
% read all the frames in order
for Frame = 1:200
    I = imread(strcat('seg_',num2str(Frame),'.jpg'));
    imshow(I);
    set(gcf,'Position',[1536*0.05 864*0.1 1536*0.8 864*0.8]);
    F(Frame) = getframe(gcf); % save figure
end

% Write to File
vW = VideoWriter('Part0_subpart3_IDseg.avi');
vW.FrameRate = 5;
open(vW)
for Frame = 1:200
    writeVideo(vW,F(Frame).cdata)
end
close(vW)
%% ------------------------------------------------------------------------

cd(ScriptsPart0Folder); % return