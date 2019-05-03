%% PART 1.
% Gaussian  Mixture  Models  and  Maximum  Likelihood  Algorithm (30) 

close all
clear
clc

current = pwd;
plot_path = '../../Output/Part1/';

cd('../../Images/TrainingSet/Frames')
ImageInfo = dir ('*.jpg'); % gather image file names from folder
Im = imread(ImageInfo(1).name); % read in first image file in folder
cd(current)

% display image
figure;imshow(Im);

% **filtering**
Im2 = imgaussfilt(Im);
% Im2 = medfilt2(Im);
I = double(Im2);

% reformat image into array
dataI = zeros(5, size(I,2)*size(I,1));
k = 1;
for j = 1:size(I,1) % all y
    for i = 1:size(I,2) % all x
        dataI(:,k) = [i;size(I,1)-j;I(j,i,1)/255;I(j,i,2)/255;I(j,i,3)/255]; % x,y,R,G,B
        k = k+1;
    end
end


% do clustering
% [idxI, CI] = K_MEANS(50, dataI, plot_path);
[muI, covarI] = EM(10, dataI, plot_path);