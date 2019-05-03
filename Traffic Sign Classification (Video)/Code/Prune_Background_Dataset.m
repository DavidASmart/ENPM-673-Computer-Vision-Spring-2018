%% PRUNE BACKGROUND DATASET
% reduce the size of the datasets to be more like those given for the signs

close all
clear
clc

%% set up
scriptsfolder = pwd;

TrainingfolderB = '../Input/Training/Blue/Background/';
TrainingfolderR = '../Input/Training/Red/Background/';
TestfolderB = '../Input/Test/Blue/Background/';
TestfolderR = '../Input/Test/Red/Background/';

n = 350; % ~ number of images per folder

%% Training-Blue
cd(TrainingfolderB);
fileinfo = dir('*.jpg'); % all jpegs

% loop through images
for f = 1:length(fileinfo)
    if rand > n/length(fileinfo) % reduce to 200 items
        delete(fileinfo(f).name)
    end
end

%% Training-Red
cd(scriptsfolder);cd(TrainingfolderR);
fileinfo = dir('*.jpg'); % all jpegs

% loop through images
for f = 1:length(fileinfo)
    if rand > n/length(fileinfo) % reduce to 200 items
        delete(fileinfo(f).name)
    end
end

%% Test-Blue
cd(scriptsfolder);cd(TestfolderB);
fileinfo = dir('*.jpg'); % all jpegs

% loop through images
for f = 1:length(fileinfo)
    if rand > n/length(fileinfo) % reduce to 200 items
        delete(fileinfo(f).name)
    end
end

%% Test-R
cd(scriptsfolder);cd(TestfolderR);
fileinfo = dir('*.jpg'); % all jpegs

% loop through images
for f = 1:length(fileinfo)
    if rand > n/length(fileinfo)  % reduce to 200 items
        delete(fileinfo(f).name)
    end
end


%% clean up
cd(scriptsfolder);

close all
clear
clc