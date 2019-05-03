%% ************************************************************************* 
% University of Maryland
% ENPM 673 Robotics Perception  
% Project 4: Traffic Sign Recognition
% David A Smart
% Due on: May 16, 2018 
% ************************************************************************* 
%
%% Traffic Sign Detection and Traffic Sign Classification. 
% This is one way to perceive the recognition pipeline. 
% In the Detection stage we aim to extract possible candidates/regions 
% which contain some traffic sign (we do not care what the sign might be). 
% In the Classification stage, we go over each Region of Interest (RoI)
% and identify what sign it represents 
% i.e., given that we know a set of traffic signs we are now classifying 
% what this specific RoI represents. 
%  
% Note: You are given images from a driving car, 
% training/testing images for various signs 
% and the output should be a video submission - DATASET. 
%  

%% 2.2 Traffic Sign Classification - 50 Pts 
% You are given sample images for diferent signs, 
% you can resize the images to a standard size say (64 × 64) 
% and extract HOG features.

% clear up memory
close all
clear
clc

% need to get data?
gather = 1;

% set up folders and files
scriptsfolder = pwd; % current folder
TrainingfolderB = '../Input/Training/Blue';
TrainingfolderR = '../Input/Training/Red';
TestingfolderB = '../Input/Test/Blue';
TestingfolderR = '../Input/Test/Red';

% get data
clc; fprintf('\n Getting Data ... \n')
if gather
    % BLUE
    [trainingFeaturesB,trainingLabelsB,testFeaturesB,testLabelsB] = Get_Data(TrainingfolderB,TestingfolderB);
    cd(scriptsfolder);
    % save data so it can be simply loaded later
    save('testFeaturesB.mat','testFeaturesB','trainingLabelsB');
    save('trainingFeaturesB.mat','trainingFeaturesB','testLabelsB');
    % RED
    [trainingFeaturesR,trainingLabelsR,testFeaturesR,testLabelsR] = Get_Data(TrainingfolderR,TestingfolderR);
    cd(scriptsfolder);
    % save data so it can be simply loaded later
    save('testFeaturesR.mat','testFeaturesR','trainingLabelsR');
    save('trainingFeaturesR.mat','trainingFeaturesR','testLabelsR');
else
    load('trainingFeaturesB.mat');
    load('testFeaturesB.mat');
    load('trainingFeaturesR.mat');
    load('testFeaturesR.mat');
end

%% Train a multi-class SVM for various signs. 
clc; fprintf('\n Training SVM ... \n');
classifierB = fitcecoc(trainingFeaturesB, trainingLabelsB);
classifierR = fitcecoc(trainingFeaturesR, trainingLabelsR);

save('svmB.mat','classifierB');
save('svmR.mat','classifierR');

%% Test the classifer performance using test data
clc; fprintf('\n Testing SVM ... \n');
predictedLabelsB = predict(classifierB, testFeaturesB);
confMatB = confusionmat(testLabelsB, predictedLabelsB);
plotconfusion(testLabelsB, predictedLabelsB);
saveas(gcf,'confusionB.jpg');
AccuracyB = mean(predictedLabelsB == testLabelsB)*100;
fprintf('\n Accuracy of BLUE SVM = %3.1f %% \n',AccuracyB)

predictedLabelsR = predict(classifierR, testFeaturesR);
confMatR = confusionmat(testLabelsR, predictedLabelsR);
plotconfusion(testLabelsR, predictedLabelsR);
saveas(gcf,'confusionR.jpg');
AccuracyR = mean(predictedLabelsR == testLabelsR)*100;
fprintf('\n Accuracy of RED SVM = %3.1f %% \n',AccuracyR)

pause(0.1);
close all
clear
clc
