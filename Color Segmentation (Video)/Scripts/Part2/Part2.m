%% PART 2.
% Color Model Learning (20) 
% Now, we extend the clustering concept to use it for our goal of color segmentation. 
%
% 1.
% Compute the average histogram for each color channel of the training set images. 
% This should help decide on the number of Gaussians to fit the color histogram.
%
% 2.
% Use the Expectation Maximization function from Part 1 to compute "mu & covar"
%
% save as EM_<Color>.jpg


% clean up workspace & memory
close all
clear
clc

% current location is ...
ScriptsPart2Folder = pwd;
% need the EM function from Part 1...
ScriptsPart1Folder = '../Part1/';
% want to saves the histogram plots and EM plots in...
plot_path = '../../Output/Part2/';

%% 1. ---------------------------------------------------------------------
% Compute the average histogram for each color channel of the training set images. 

[Y_data, R_data, G_data] = averageHistogram2;

cd(plot_path)
load('Y_data');
load('R_data');
load('G_data');
cd(ScriptsPart2Folder);

%% 2. ---------------------------------------------------------------------
% Use the Expectation Maximization function from Part 1 to compute "mu & covar"

% randomly select datapoints ... to reduce size
while size(Y_data,2) > 10000
    idY = round(rand(1,size(Y_data,2))); 
    Y_data = Y_data(:,idY == 1);
end

while size(R_data,2) > 10000
    idR = round(rand(1,size(R_data,2))); 
    R_data = R_data(:,idR == 1);
end

while size(G_data,2) > 10000
    idG = round(rand(1,size(G_data,2))); 
    G_data = G_data(:,idG == 1);
end

% convert to double data-type
Y_data = double(Y_data);
R_data = double(R_data);
G_data = double(G_data);

% remove black pixels which should not have been included in the cropped
% images to begin with...not sure where those came from
for i = size(Y_data,2):-1:1
    if sqrt(Y_data(1,i)^2 +  Y_data(2,i)^2 + Y_data(3,i)^2) < sqrt(3*(75^2))
        Y_data(:,i) = [];
    end
end
for i = size(R_data,2):-1:1
    if sqrt(R_data(1,i)^2 +  R_data(2,i)^2 + R_data(3,i)^2) < sqrt(3*(75^2))
        R_data(:,i) = [];
    end
end
for i = size(G_data,2):-1:1
    if sqrt(G_data(1,i)^2 +  G_data(2,i)^2 + G_data(3,i)^2) < sqrt(3*(75^2))
        G_data(:,i) = [];
    end
end

cd(plot_path); % change to folder specified
% plot data for refrence
figure; title('Yellow Buoy RGB Space'); hold on
for i = 1:size(Y_data,2) % for all data
    % Plot data in RGB space as its color 
    plot3(Y_data(1,i),Y_data(2,i),Y_data(3,i),'.','Color',[Y_data(1,i),Y_data(2,i),Y_data(3,i)]./255); 
end
view(45,45); xlabel('R');ylabel('G');zlabel('B');
saveas(gcf,'YellowBouy_RGB3Ddata.jpg'); % save

% plot data for refrence
figure; title('Red Buoy RGB Space'); hold on
for i = 1:size(R_data,2) % for all data
    % Plot data in RGB space as its color 
    plot3(R_data(1,i),R_data(2,i),R_data(3,i),'.','Color',[R_data(1,i),R_data(2,i),R_data(3,i)]./255); 
end
view(45,45); xlabel('R');ylabel('G');zlabel('B');
saveas(gcf,'RedBouy_RGB3Ddata.jpg'); % save

% plot data for refrence
figure; title('Green Buoy RGB Space'); hold on
for i = 1:size(G_data,2) % for all data
    % Plot data in RGB space as its color 
    plot3(G_data(1,i),G_data(2,i),G_data(3,i),'.','Color',[G_data(1,i),G_data(2,i),G_data(3,i)]./255); 
end
view(45,45); xlabel('R');ylabel('G');zlabel('B');
saveas(gcf,'GreenBouy_RGB3Ddata.jpg'); % save


cd(ScriptsPart2Folder); % return to the scripts folder


% Yellow Bouy ----------------------------------------------------------------

% make model
NY = 2; % number of gaussians
cd(ScriptsPart2Folder); % first switch to original folder
cd(ScriptsPart1Folder); % then switch to other folder
[muY, covarY] = EM(NY, double(Y_data), ''); % EM

% plot models
x = 1:256;
figure; hold on;
for i = 1:size(Y_data,2) % for all data
    for j = 1:NY % for each gaussian
        % Plot data in RBG space according to probability
        pYm(j) = mvnpdf(muY(:,j),muY(:,j),covarY{j});
        pY(j) = mvnpdf(Y_data(:,i),muY(:,j),covarY{j});
        p(j) = pY(j)/pYm(j);
    end
    c = [1 1 1] + max(p).*[0 0 -1]; % interpolate between white and yellow
    plot3(Y_data(1,i),Y_data(2,i),Y_data(3,i),'.','Color',c); 
end
view(-45,45); xlabel('R');ylabel('G');zlabel('B');
title('EM-Y'); hold off;
cd(plot_path) % switch to output folder
saveas(gcf,'EM_Y.jpg'); % save as EM_<Color>.jpg

% Red Bouy ----------------------------------------------------------------

% make model
NR = 2; % number of gaussians
cd(ScriptsPart2Folder); % first switch to original folder
cd(ScriptsPart1Folder); % then switch to other folder
[muR, covarR] = EM(NR, double(R_data), ''); % EM

% plot models
figure; hold on;
for i = 1:size(R_data,2) % for all data
    for j = 1:NR % for each gaussian
        % Plot data in RBG space according to probability
        pRm(j) = mvnpdf(muR(:,j),muR(:,j),covarR{j});
        pR(j) = mvnpdf(R_data(:,i),muR(:,j),covarR{j});
        p(j) = pR(j)/pRm(j);
    end
    c = [1 1 1] + max(p).*[0 -1 -1]; % interpolate between white and red
    plot3(R_data(1,i),R_data(2,i),R_data(3,i),'.','Color',c); 
end
view(-45,45); xlabel('R');ylabel('G');zlabel('B');
title('EM-R'); hold off;
cd(plot_path) % switch to output folder
saveas(gcf,'EM_R.jpg'); % save as EM_<Color>.jpg

% Green Bouy ----------------------------------------------------------------

% make model
NG = 2; % number of gaussians
cd(ScriptsPart2Folder); % first switch to original folder
cd(ScriptsPart1Folder); % then switch to other folder
[muG, covarG] = EM(NG, double(G_data), ''); % EM

% plot models
x = 1:256;
figure; hold on;
for i = 1:size(G_data,2) % for all data
    for j = 1:NG % for each gaussian
        % Plot data in RBG space according to probability
        pGm(j) = mvnpdf(muG(:,j),muG(:,j),covarG{j});
        pG(j) = mvnpdf(G_data(:,i),muG(:,j),covarG{j});
        p(j) = pG(j)/pGm(j);
    end
    c = [1 1 1] + max(p).*[-1 0 -1]; % interpolate between white and green
    plot3(G_data(1,i),G_data(2,i),G_data(3,i),'.','Color',c);
end
view(-45,45); xlabel('R');ylabel('G');zlabel('B');
title('EM-G'); hold off;
cd(plot_path) % switch to output folder
saveas(gcf,'EM_G.jpg'); % save as EM_<Color>.jpg


% ALL ------------------------------------------------------------------------

figure; hold on;
for i = 1:size(Y_data,2) % for all data
    for j = 1:NY % for each gaussian
        % Plot data in RBG space according to probability
        pYm(j) = mvnpdf(muY(:,j),muY(:,j),covarY{j});
        pY(j) = mvnpdf(Y_data(:,i),muY(:,j),covarY{j});
        p(j) = pY(j)/pYm(j);
    end
    c = [1 1 1] + max(p).*[0 0 -1]; % interpolate between white and yellow
    plot3(Y_data(1,i),Y_data(2,i),Y_data(3,i),'.','Color',c); 
end
for i = 1:size(R_data,2) % for all data
    for j = 1:NR % for each gaussian
        % Plot data in RBG space according to probability
        pRm(j) = mvnpdf(muR(:,j),muR(:,j),covarR{j});
        pR(j) = mvnpdf(R_data(:,i),muR(:,j),covarR{j});
        p(j) = pR(j)/pRm(j);
    end
    c = [1 1 1] + max(p).*[0 -1 -1]; % interpolate between white and red
    plot3(R_data(1,i),R_data(2,i),R_data(3,i),'.','Color',c); 
end
for i = 1:size(G_data,2) % for all data
    for j = 1:NG % for each gaussian
        % Plot data in RBG space according to probability
        pGm(j) = mvnpdf(muG(:,j),muG(:,j),covarG{j});
        pG(j) = mvnpdf(G_data(:,i),muG(:,j),covarG{j});
        p(j) = pG(j)/pGm(j);
    end
    c = [1 1 1] + max(p).*[-1 0 -1]; % interpolate between white and green
    plot3(G_data(1,i),G_data(2,i),G_data(3,i),'.','Color',c);
end
view(-45,45); xlabel('R');ylabel('G');zlabel('B');
title('EM-Y,R,&G'); hold off;
% save
saveas(gcf,'EM_RGB.jpg'); % save as EM_<Color>.jpg 

% also save numerical values...
save('EM.mat','muY','muR','muG','covarY','covarR','covarG');

%% ------------------------------------------------------------------------
cd(ScriptsPart2Folder); % switch to original folder