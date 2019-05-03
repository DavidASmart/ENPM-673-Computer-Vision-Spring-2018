%% averageHistogram.m
%% PART 0.
% 2.
% For each buoy, compute and visualize the average histogram for each color 
% channel of the cropped training dataset images.

% clean up workspace & memory
close all
clear
clc

% current location is P3_submission/ColorSeg/Scripts/Part0
ScriptsPart0Folder = pwd;
% want to read cropped images from ColorSeg/Images/TrainingSet/CroppedBuoys/
CroppedTrainingSetFolder = '../../Images/TrainingSet/CroppedBuoys';
% want to saves the histogram plots in ColorSeg/Output/Part0/
HistogramFolder = '../../Output/Part0/';

cd(CroppedTrainingSetFolder) % change to input folder
ImageInfo = dir ('*.jpg'); % gather image file names from folder

% initialize histogram information
GR = zeros(256,1);GG = zeros(256,1);GB = zeros(256,1);GC = 0;
RR = zeros(256,1);RG = zeros(256,1);RB = zeros(256,1);RC = 0;
YR = zeros(256,1);YG = zeros(256,1);YB = zeros(256,1);YC = 0;

for k = 1:length(ImageInfo) % for all images
    
    I = imread(ImageInfo(k).name); % read in image files
    Iname = ImageInfo(k).name; % image file name
    
    % **filtering**
    I = imgaussfilt(I);
    % I = medfilt2(I);
    
    
    BW = imbinarize(rgb2gray(I)); % distinguish between colored & none-colored pixels
    PixelList = regionprops(BW,'PixelList'); % locations of colored & none-colored pixels
    if ~isempty(PixelList) % if colored pixels exist
       
        PixelList = PixelList.PixelList; % extract locations of colored pixels

        R = imhist(I(PixelList(:,2),PixelList(:,1),1)); % red component
        G = imhist(I(PixelList(:,2),PixelList(:,1),2)); % green component
        B = imhist(I(PixelList(:,2),PixelList(:,1),3)); % blue component

        % sum over all images
        if Iname(1) == 'G' % if green bouy
            GR = GR + R;
            GG = GG + G;
            GB = GB + B;
            GC = GC + 1;
        elseif Iname(1) == 'R' % if red bouy
            RR = RR + R;
            RG = RG + G;
            RB = RB + B;
            RC = RC + 1;
        elseif Iname(1) == 'Y' % if yellow bouy
            YR = YR + R;
            YG = YG + G;
            YB = YB + B;
            YC = YC + 1;
        end

    end
end

% average historgrams
GR = GR/GC; GG = GG/GC; GB = GB/GC;
RR = RR/RC; RG = RG/RC; RB = RB/RC;
YR = YR/YC; YG = YG/YC; YB = YB/YC;

% save as <Color>_hist.jpg
cd(ScriptsPart0Folder) % first switch to original folder
cd(HistogramFolder) % then switch to output folder

figure(1); plot(GR,'r'); hold on, plot(GG,'g'); plot(GB,'b'); 
legend(' Red channel','Green channel','Blue channel'); hold off;
title('GREEN Bouy');
saveas(gcf,'G_hist.jpg')

figure(2); plot(RR,'r'); hold on, plot(RG,'g'); plot(RB,'b'); 
legend(' Red channel','Green channel','Blue channel'); hold off;
title('RED Bouy');
saveas(gcf,'R_hist.jpg')

figure(3); plot(YR,'r'); hold on, plot(YG,'g'); plot(YB,'b'); 
legend(' Red channel','Green channel','Blue channel'); hold off;
title('YELLOW Bouy');
saveas(gcf,'Y_hist.jpg')

% also save numerical values...
save('G_hist.mat','GR','GG','GB');
save('R_hist.mat','RR','RG','RB');
save('Y_hist.mat','YR','YG','YB');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(ScriptsPart0Folder); % return