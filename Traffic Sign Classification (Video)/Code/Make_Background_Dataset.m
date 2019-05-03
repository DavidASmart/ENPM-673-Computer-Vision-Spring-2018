%% MAKE BACKGROUND DATSET
% required for good SVM classification

close all
clear
clc


% set up folders
scriptsfolder = pwd;
inputfolder = '../Input/Video/';
TrainingfolderB = '../Input/Training/Blue/Background/';
TrainingfolderR = '../Input/Training/Red/Background/';
TestfolderB = '../Input/Test/Blue/Background/';
TestfolderR = '../Input/Test/Red/Background/';

cd(inputfolder);
fileinfo = dir('*.jpg'); % all jpegs

% loop through images
inc = 5;
for f = 1:inc:length(fileinfo) % 2799
    
    % print out
    clc; fprintf('\n frame: %i', f);
    
    cd(scriptsfolder); cd(inputfolder); % switch to video input folder
    Im = imread(fileinfo(f).name); % read in image
    cd(scriptsfolder);
    
    % find regions which the detector will pick out
    [Blue_BB, Red_BB,~] = HSV_Detect(Im, 0);
    [Blue_BB, Red_BB,~] = MSER_Detect(Im, Blue_BB, Red_BB, 0);
    
    % crop and save each accoring to color...
    % i will go through and remove the signs (true detections) after
    % BLUE ----------------------------------------------------------------
    if ~isempty(Blue_BB)
        for i = 1:size(Blue_BB,1) % for each blue region
            
            % get the bouding box cordinates
            B_bb = round(Blue_BB(i,:)); % [x,y,w,h]
            
            % crop what was detected
            J = imcrop(Im,B_bb);
            
            % switch to one of the dataset folders
            if  mod(f,2*inc) == 1
                cd(scriptsfolder); cd(TrainingfolderB);
            else % mod(f,20) == 0
                cd(scriptsfolder); cd(TestfolderB);
            end
            
            % save the cropped image as part of the dataset
            imwrite(J,strcat(num2str(f),'_',num2str(i),'.jpg'));
            
        end
    end
    % RED -----------------------------------------------------------------
    if ~isempty(Red_BB)
        for i = 1:size(Red_BB,1) % for each blue region
            
            % get the bouding box cordinates
            R_bb = round(Red_BB(i,:)); % [x,y,w,h]
            
            % crop what was detected
            J = imcrop(Im,R_bb);
            
            % switch to one of the dataset folders
            if  mod(f,2*inc) == 1
                cd(scriptsfolder); cd(TrainingfolderR);
            else % mod(f,20) == 0
                cd(scriptsfolder); cd(TestfolderR);
            end
            
            % save the cropped image as part of the dataset
            imwrite(J,strcat(num2str(f),'_',num2str(i),'.jpg'));
        end
    end
    
end

cd(scriptsfolder);

close all
clear
clc