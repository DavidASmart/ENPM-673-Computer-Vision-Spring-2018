%% PART 1.
% Gaussian  Mixture  Models  and  Maximum  Likelihood  Algorithm (30) 

close all
clear
clc

plot_path = '../../Output/Part1/';

%% 1.
% Generate 100 data samples (each) from 3 1-D Gaussians, with different 
% means and variances. 

% set up disribution parameters
mu = rand(3,1)*100;
sigma = rand(3,1)*25/2;

% make gaussian distribution models
gm1 = gmdistribution(mu(1),sigma(1));
gm2 = gmdistribution(mu(2),sigma(2));
gm3 = gmdistribution(mu(3),sigma(3));

% generate data
n = 100;
data(1:n) = random(gm1,n);
data(n+1:2*n) = random(gm2,n);
data(2*n+1:3*n) = random(gm3,n);

% plot data for refrence
figure; title('Original 1D-Data Clusters'); hold on
plot(data(1:n),zeros(1,n),'.'); % plot data in each culster as a distinct color
plot(data(n+1:2*n),zeros(1,n),'.'); % plot data in each culster as a distinct color
plot(data(2*n+1:3*n),zeros(1,n),'.'); % plot data in each culster as a distinct color
plot(mu(:),0,'ko','MarkerFaceColor','k'); % plot cluster centers

% save
currentfolder = pwd; % save the current folder so that it can be returned to after saving figure
cd(plot_path); % change to folder specified
saveas(gcf,'Original1D3N.jpg');
cd(currentfolder); % return to the scripts folder

%% 2.
% Implement the Expectation Maximization (EM) algorithm.
% Refer to http://www.cs.ubbcluj.ro/~csatol/gep_tan/Bishop-CUED-2006.pdf
% 
% Submit the script ColorSeg/Scripts/Part1/EM.m 

%% 3.
% Recover the model parameters (means and variances) for 3 Gaussians 
% from the data created in step 1 using your implementation of Expectation 
% Maximization created in step 2. 

% [idx1d3n, C1d3n] = K_MEANS(3, data, plot_path);
[mu1d3n, covar1d3n] = EM(3,data,plot_path);

%% 4.
% Now model the data from step one use 4 Gaussians rather than 3.
%
% Save the plot of the computed Gaussians as
% ColorSeg/Output/Part1/EM1D4N.jpg

% [idx1d4n, C1d4n] = K_MEANS(4, data, plot_path);
[mu1d4n, covar1d4n] = EM(4,data,plot_path);

% Discuss what you observe.
% >> see P1_submission/Reports/ColorSeg.pdf

clc
fprintf('\n\n ORIGINAL Distribution Parameters:');
fprintf('\n means: ');
fprintf('\n\t m1: %f ',mu(1));
fprintf('\n\t m2: %f ',mu(2));
fprintf('\n\t m3: %f ',mu(3));
fprintf('\n covariance: ');
fprintf('\n\t c1: %f ',sigma(1));
fprintf('\n\t c2: %f ',sigma(2));
fprintf('\n\t c3: %f ',sigma(3));

fprintf('\n\n 3-Gaussian Estimated Distribution Parameters:');
fprintf('\n means: ');
fprintf('\n\t m1: %f ',mu1d3n(1));
fprintf('\n\t m2: %f ',mu1d3n(2));
fprintf('\n\t m3: %f ',mu1d3n(3));
fprintf('\n covariance: ');
fprintf('\n\t c1: %f ',covar1d3n{1});
fprintf('\n\t c2: %f ',covar1d3n{2});
fprintf('\n\t c3: %f ',covar1d3n{3});

fprintf('\n\n 4-Gaussian Estimated Distribution Parameters:');
fprintf('\n means: ');
fprintf('\n\t m1: %f ',mu1d4n(1));
fprintf('\n\t m2: %f ',mu1d4n(2));
fprintf('\n\t m3: %f ',mu1d4n(3));
fprintf('\n\t m4: %f ',mu1d4n(4));
fprintf('\n covariance: ');
fprintf('\n\t c1: %f ',covar1d4n{1});
fprintf('\n\t c2: %f ',covar1d4n{2});
fprintf('\n\t c3: %f ',covar1d4n{3});
fprintf('\n\t c4: %f ',covar1d4n{4});
fprintf('\n\n');
%% 5.
% Confirm your implementation for generalized D-Dimensional Gaussians.

% 2D Data -----------------------------------------------------------------
% set up disribution parameters
mu = rand(3,2)*50;
sigma = rand(3,2)*25;

% make gaussian distribution models
gm1 = gmdistribution(mu(1,:),sigma(1,:));
gm2 = gmdistribution(mu(2,:),sigma(2,:));
gm3 = gmdistribution(mu(3,:),sigma(3,:));

% generate data
n = 100;
data2d(:,1:n) = random(gm1,n)';
data2d(:,n+1:2*n) = random(gm2,n)';
data2d(:,2*n+1:3*n) = random(gm3,n)';

% plot data for refrence
figure; title('Original 2D-Data Clusters'); hold on
plot(data2d(1,1:n),data2d(2,1:n),'.'); % plot data in each culster as a distinct color
plot(data2d(1,n+1:2*n),data2d(2,n+1:2*n),'.'); % plot data in each culster as a distinct color
plot(data2d(1,2*n+1:3*n),data2d(2,2*n+1:3*n),'.'); % plot data in each culster as a distinct color
plot(mu(:,1),mu(:,2),'ko','MarkerFaceColor','k'); % plot cluster centers

% save
currentfolder = pwd; % save the current folder so that it can be returned to after saving figure
cd(plot_path); % change to folder specified
saveas(gcf,'Original2D3N.jpg');
cd(currentfolder); % return to the scripts folder

% [idx2D, C2D] = K_MEANS(3, data2d, plot_path);
[mu2D, covar2D] = EM(3, data2d, plot_path);

% 3D Data -----------------------------------------------------------------

% set up disribution parameters
mu = rand(3,3)*25;
sigma = rand(3,3)*25/2;

% make gaussian distribution models
gm1 = gmdistribution(mu(1,:),sigma(1,:));
gm2 = gmdistribution(mu(2,:),sigma(2,:));
gm3 = gmdistribution(mu(3,:),sigma(3,:));

% generate data
n = 100;
data3d(:,1:n) = random(gm1,n)';
data3d(:,n+1:2*n) = random(gm2,n)';
data3d(:,2*n+1:3*n) = random(gm3,n)';

% plot data for refrence
figure; title('Original 3D-Data Clusters'); hold on
plot3(data3d(1,1:n),data3d(2,1:n),data3d(3,1:n),'.'); % plot data in each culster as a distinct color
plot3(data3d(1,n+1:2*n),data3d(2,n+1:2*n),data3d(3,n+1:2*n),'.'); % plot data in each culster as a distinct color
plot3(data3d(1,2*n+1:3*n),data3d(2,2*n+1:3*n),data3d(3,2*n+1:3*n),'.'); % plot data in each culster as a distinct color
plot3(mu(:,1),mu(:,2),mu(:,3),'ko','MarkerFaceColor','k'); % plot cluster centers
view(45,45);

% save
currentfolder = pwd; % save the current folder so that it can be returned to after saving figure
cd(plot_path); % change to folder specified
saveas(gcf,'Original3D3N.jpg');
cd(currentfolder); % return to the scripts folder

% [idx3D, C3D] = K_MEANS(3, data3d, plot_path);
[mu3D, covar3D] = EM(3, data3d, plot_path);

% 5D Data! -----------------------------------------------------------------

% set up disribution parameters
mu(:,1:2) = rand(3,2)*50; % for position values
sigma(:,1:2) = rand(3,2)*25/2; % for position values
mu(:,3:5) = rand(3,3); % for color values
sigma(:,3:5) = rand(3,3)*0.1; % for color values

% make gaussian distribution models
gm1 = gmdistribution(mu(1,:),sigma(1,:));
gm2 = gmdistribution(mu(2,:),sigma(2,:));
gm3 = gmdistribution(mu(3,:),sigma(3,:));

% generate data
n = 100;
data5d(:,1:n) = random(gm1,n)';
data5d(:,n+1:2*n) = random(gm2,n)';
data5d(:,2*n+1:3*n) = random(gm3,n)';

% force color values to fit in range by hard threshold
for i = 1:3*n
    for j = 3:5
        if data5d(j,i) > 1
            data5d(j,i) = 1;
        elseif data5d(j,i) < 0
            data5d(j,i) = 0;
        end
    end
end

% plot data for refrence
figure; title('Original 5D-Data Clusters'); hold on
for i = 1:3*n
    plot(data5d(1,i),data5d(2,i),'s','MarkerFaceColor',[data5d(3,i) data5d(4,i) data5d(5,i)]); % plot data in each culster as true color
end
plot(mu(:,1),mu(:,2),'ko','MarkerFaceColor','k'); % plot cluster centers

% save
currentfolder = pwd; % save the current folder so that it can be returned to after saving figure
cd(plot_path); % change to folder specified
saveas(gcf,'Original5D3N.jpg');
cd(currentfolder); % return to the scripts folder

% [idx5D, C5D] = K_MEANS(3, data5d, plot_path);
[mu5D, covar5D] = EM(3, data5d, plot_path);
  
