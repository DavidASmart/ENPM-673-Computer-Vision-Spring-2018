%% ************************************************************************* 
% University of Maryland
% ENPM 673 Robotics Perception  
% Project 4: Traffic Sign Recognition
% David A Smart
% Due on: May 16, 2018 
% ************************************************************************* 
%
function [trainingFeatures,trainingLabels,testFeatures,testLabels] = Get_Data(Trainingfolder,Testingfolder)
% my previous function did not format the data correctly...
% I am now following the example for digit classification

    % get data
    trainingSet = imageDatastore(Trainingfolder,   'IncludeSubfolders', true, 'LabelSource', 'foldernames');
    testSet     = imageDatastore(Testingfolder, 'IncludeSubfolders', true, 'LabelSource', 'foldernames');
    
    %% Training Set
    numImages = numel(trainingSet.Files);
    cellSize = [2 2];

    for i = 1:numImages
        img = readimage(trainingSet, i); % read in image
        %img = imadjust(img,stretchlim(img),[]); % contrast normalization
        img = imresize(img,[64 64]); % resize to standard 64x64
        %img = rgb2gray(img);
        %img = imbinarize(img);
        trainingFeatures(i, :) = extractHOGFeatures(img, 'CellSize', cellSize);  
    end
    
    % Get labels for each image.
    trainingLabels = trainingSet.Labels;
    
    %% Testing Set
    numImages = numel(testSet.Files);

    for i = 1:numImages
        img = readimage(testSet, i); % read in image
        %img = imadjust(img,stretchlim(img),[]); % contrast normalization
        img = imresize(img,[64 64]); % resize to standard 64x64
        %img = rgb2gray(img);
        %img = imbinarize(img);
        testFeatures(i, :) = extractHOGFeatures(img, 'CellSize', cellSize);  
    end
    
    % Get labels for each image.
    testLabels = testSet.Labels;
    
end
